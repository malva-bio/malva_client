import json
import time
import requests
import msgpack
import numpy as np
import pandas as pd
from typing import Dict, List, Any, Optional, Union, Tuple
from pathlib import Path
import logging
from urllib.parse import urljoin, urlparse, urlunparse
import h5py
import anndata as ad

from malva_client.exceptions import MalvaAPIError, AuthenticationError, SearchError, QuotaExceededError
from malva_client.models import (
    CellExpressionMatrixResult,
    SearchResult,
    CoverageResult,
    CoexpressionResult,
    UMAPCoordinates,
)
from malva_client.config import Config


def _expr_data_from_columnar(expr_col: dict) -> dict:
    """Reconstruct row-major expression_data from compact_v2 columnar format."""
    si = expr_col.get('si', [])
    n = len(si)
    if n == 0:
        return {
            'data': [],
            'samples': expr_col.get('samples', []),
            'cell_types': expr_col.get('cell_types', []),
            'columns': expr_col.get('columns', []),
            'celltype_sample_counts': expr_col.get('celltype_sample_counts', {}),
        }
    ci = expr_col.get('ci', [])
    v0 = expr_col.get('v0', [])
    v1 = expr_col.get('v1', [])
    v2 = expr_col.get('v2', [])
    v3 = expr_col.get('v3', [])
    v4 = expr_col.get('v4', [])
    return {
        'data': [[si[j], ci[j], v0[j], v1[j], v2[j], v3[j], v4[j]] for j in range(n)],
        'samples': expr_col.get('samples', []),
        'cell_types': expr_col.get('cell_types', []),
        'columns': expr_col.get('columns', []),
        'celltype_sample_counts': expr_col.get('celltype_sample_counts', {}),
    }


class GeneResult:
    """Lazy per-gene view over compact_v2 data; per-cell arrays expanded on access."""

    def __init__(self, gdata: dict, global_cells, global_samples):
        self._d = gdata
        self._gc = global_cells
        self._gs = global_samples

    @property
    def ncells(self) -> int:
        return len(self._d.get('_ci') or self._d.get('cell') or [])

    @property
    def expression_data(self) -> dict:
        ed = self._d.get('expression_data', {})
        if isinstance(ed, dict) and ed.get('_fmt') == 'col':
            return _expr_data_from_columnar(ed)
        return ed

    @property
    def cells(self) -> np.ndarray:
        ci = self._d.get('_ci')
        if ci is not None and self._gc is not None:
            return np.asarray(self._gc, dtype=np.uint32)[np.asarray(ci, dtype=np.int32)]
        return np.asarray(self._d.get('cell', []), dtype=np.uint64)

    @property
    def samples(self) -> np.ndarray:
        ci = self._d.get('_ci')
        if ci is not None and self._gs is not None:
            return np.asarray(self._gs, dtype=np.uint64)[np.asarray(ci, dtype=np.int32)]
        return np.asarray(self._d.get('sample', []), dtype=np.uint64)

    @property
    def expression(self) -> np.ndarray:
        return np.asarray(self._d.get('expression', []), dtype=np.float32)

    # Allow dict-style access for backward-compatible code that does result['cell']
    def __getitem__(self, key: str):
        if key == 'cell':
            return self.cells.tolist()
        if key == 'sample':
            return self.samples.tolist()
        if key == 'expression':
            return self.expression.tolist()
        if key == 'expression_data':
            return self.expression_data
        if key == 'ncells':
            return self.ncells
        return self._d[key]

    def get(self, key: str, default=None):
        try:
            return self[key]
        except KeyError:
            return default


class SearchResults:
    """Lazy wrapper for compact_v2 search results.  Instantiation is O(1)."""

    def __init__(self, data: dict):
        self._data = data          # full top-level response dict (job_id, status, …)
        r = data.get('results', data)
        self._compact = isinstance(r, dict) and r.get('_format') == 'compact_v2'
        self._gc = r.get('_global_cells') if self._compact else None
        self._gs = r.get('_global_samples') if self._compact else None
        self._r = r
        self.sample_metadata = r.get('_sample_metadata', {})

    # Allow dict-style field access for top-level metadata (status, job_id, …)
    def get(self, key, default=None):
        return self._data.get(key, default)

    @property
    def genes(self) -> list:
        return [k for k, v in self._r.items() if not k.startswith('_') and isinstance(v, dict)]

    def __getitem__(self, key):
        # Gene access (primary use case)
        if key in self._r and not key.startswith('_') and isinstance(self._r.get(key), dict):
            return GeneResult(self._r[key], self._gc, self._gs)
        # Top-level metadata access (job_id, status, …)
        return self._data[key]

    def __contains__(self, gene: str) -> bool:
        return gene in self._r and not gene.startswith('_')

    def items(self):
        for g in self.genes:
            yield g, self[g]

    def __iter__(self):
        return iter(self.genes)


def _reconstruct_compact_v2(data: dict) -> dict:
    """Reconstruct compact_v2 search result to dense format in-place.

    Kept for backward compatibility; prefer SearchResults for new code.
    """
    if not isinstance(data, dict) or data.get('_format') != 'compact_v2':
        return data

    global_cells = data.get('_global_cells', [])
    global_samples = data.get('_global_samples', [])
    results = data.get('results', {})

    for gene, gdata in results.items():
        if gene.startswith('_') or not isinstance(gdata, dict):
            continue
        ci = gdata.pop('_ci', None)
        if ci:
            gdata['cell'] = [global_cells[i] for i in ci]
            gdata['sample'] = [global_samples[i] for i in ci]
        else:
            gdata['cell'] = []
            gdata['sample'] = []
        ed = gdata.get('expression_data')
        if isinstance(ed, dict) and ed.get('_fmt') == 'col':
            gdata['expression_data'] = _expr_data_from_columnar(ed)

    data.pop('_format', None)
    data.pop('_global_cells', None)
    data.pop('_global_samples', None)
    return data

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class MalvaClient:
    """
    Main client for interacting with the Malva API
    
    Provides methods for searching, retrieving samples, and downloading data files.
    """

    def __init__(self, base_url: str = None, api_token: str = None, timeout: int = 300,
                 verify_ssl: bool = True):
        """
        Initialize the Malva client
        
        Args:
            base_url: Base URL of the Malva API  
            api_token: API token for authentication
            timeout: Request timeout in seconds
            verify_ssl: Whether to verify SSL certificates
        """
        # Load config if no explicit parameters provided
        if base_url is None or api_token is None:
            config = Config.load()
            base_url = base_url or config.server_url
            api_token = api_token or config.api_token
            if verify_ssl is None:
                verify_ssl = config.verify_ssl
        
        if not base_url:
            raise ValueError("base_url must be provided either directly or via config/environment")
        
        self.base_url = base_url.rstrip('/')
        self.api_token = api_token
        self.timeout = timeout
        self.verify_ssl = verify_ssl
        
        self.session = requests.Session()
        self.session.verify = verify_ssl
        
        if api_token:
            self.session.headers.update({'Authorization': f'Bearer {api_token}'})
        
        self._test_connection()

    
    def _test_connection(self):
        """Test connection and authentication with better error handling"""
        try:
            response = self._request('GET', '/health')
            logger.info(f"Connected to Malva API: {response.get('status', 'unknown')}")
            
            # Test authentication if token provided
            if self.api_token:
                try:
                    quota_response = self._request('GET', '/api/quota-status')
                    logger.info(f"Authenticated successfully: {quota_response.get('account_type', 'unknown')} account")
                except AuthenticationError:
                    logger.warning("API token authentication failed")
                    raise
                except MalvaAPIError as e:
                    if "404" in str(e) or "not found" in str(e).lower():
                        logger.info("Connected but quota endpoint not available")
                    else:
                        raise
                        
        except MalvaAPIError as e:
            if "not connect" in str(e) or "404" in str(e):
                raise MalvaAPIError(f"Could not connect to Malva API at {self.base_url}. Please check the URL and ensure the server is running.")
            raise

    def get_quota_status(self) -> Dict[str, Any]:
        """Get current quota status"""
        if not self.api_token:
            raise AuthenticationError("API token required to check quota status")
        return self._request('GET', '/api/quota-status')
    
    def is_authenticated(self) -> bool:
        """Check if the client is properly authenticated"""
        try:
            self.get_quota_status()
            return True
        except (AuthenticationError, MalvaAPIError):
            return False
    
    def _request(self, method: str, endpoint: str, **kwargs) -> Dict[str, Any]:
        """Make a request to the API with enhanced error handling"""
        url = urljoin(self.base_url, endpoint.lstrip('/'))
        
        try:
            response = self.session.request(method, url, timeout=self.timeout, **kwargs)
            
            logger.debug(f"{method} {url} -> {response.status_code}")
            logger.debug(f"Response content type: {response.headers.get('content-type', 'unknown')}")
            
            if response.status_code == 404:
                content_type = response.headers.get('content-type', '')
                if 'text/html' in content_type:
                    raise MalvaAPIError(f"Endpoint not found: {endpoint}. Check if the Malva API is running at {self.base_url}")

                try:
                    error_data = response.json()
                    error_msg = error_data.get('detail', error_data.get('error', 'Endpoint not found'))
                    raise MalvaAPIError(f"Not found: {error_msg}")
                except (json.JSONDecodeError, ValueError):
                    raise MalvaAPIError(f"Endpoint not found: {endpoint}")

            # Handle quota exceeded
            elif response.status_code == 429:
                try:
                    error_data = response.json()
                    quota_info = {
                        'account_type': error_data.get('account_type'),
                        'quota_exceeded': error_data.get('quota_exceeded', True)
                    }
                    raise QuotaExceededError(
                        error_data.get('error', 'Search quota exceeded'),
                        quota_info
                    )
                except (json.JSONDecodeError, ValueError):
                    raise QuotaExceededError("Search quota exceeded")

            elif response.status_code in [401, 403]:
                try:
                    error_data = response.json()
                    error_msg = error_data.get('error', error_data.get('detail', 'Authentication failed'))
                except (json.JSONDecodeError, ValueError):
                    error_msg = "Authentication failed"
                raise AuthenticationError(error_msg)

            elif not response.ok:
                try:
                    error_data = response.json()
                    error_msg = error_data.get('detail', error_data.get('error', f'HTTP {response.status_code}'))
                    raise MalvaAPIError(f"API error: {error_msg}")
                except (json.JSONDecodeError, ValueError):
                    body = response.text[:200] if response.text else '(empty body)'
                    raise MalvaAPIError(f"HTTP {response.status_code}: {body}")

            try:
                result = response.json()
                # Return a lazy SearchResults wrapper for compact_v2 to avoid
                # O(N_cells × N_genes) list expansion at parse time.
                if isinstance(result, dict) and result.get('results', {}).get('_format') == 'compact_v2':
                    return SearchResults(result)
                return result
            except (json.JSONDecodeError, ValueError) as e:
                content_type = response.headers.get('content-type', 'unknown')
                body = response.text[:200] if response.text else '(empty body)'
                raise MalvaAPIError(
                    f"Server returned non-JSON response (HTTP {response.status_code}, "
                    f"{content_type}): {body}"
                )
            
        except requests.exceptions.Timeout:
            raise MalvaAPIError(f"Request timed out after {self.timeout} seconds")
        except requests.exceptions.ConnectionError:
            raise MalvaAPIError(f"Could not connect to {self.base_url}. Check if the server is running.")
        except requests.exceptions.RequestException as e:
            raise MalvaAPIError(f"Request failed: {e}")

    def _build_search_payload(self, query: str, wait_for_completion: bool,
                              max_wait: int, poll_interval: int,
                              aggregate_expression: bool,
                              stranded: Optional[bool],
                              min_kmer_presence: int,
                              max_kmer_presence: int) -> dict:
        """Build the POST /search payload with correctly mapped parameters.

        Always submits the job asynchronously (wait_for_completion=False at
        the HTTP level) to avoid holding a long-lived connection that would
        hit nginx proxy_read_timeout on slow queries.  The client-level
        wait_for_completion flag is handled by polling after submission.
        """
        payload = {
            'query': query,
            'wait_for_completion': False,   # never block the HTTP POST
            'aggregate_expression': aggregate_expression,
            'min_kmer_presence': min_kmer_presence,
            'max_kmer_presence': max_kmer_presence,
        }
        # stranded=True  → forward strand only  → unstranded_mode=False
        # stranded=False → both strands         → unstranded_mode=True
        # stranded=None  → use server default   → omit (server defaults to forward only)
        if stranded is True:
            payload['unstranded_mode'] = False
        elif stranded is False:
            payload['unstranded_mode'] = True
        return payload

    def _get_job_status_lightweight(self, job_id: str) -> Dict[str, Any]:
        """Fetch job status without forcing the server to return result data."""
        try:
            return self._request('GET', f'/search/{job_id}/status')
        except MalvaAPIError as exc:
            if self._is_not_found_error(exc):
                return self._request('GET', f'/search/{job_id}')
            raise

    @staticmethod
    def _adaptive_poll_delay(elapsed: float, poll_interval: int) -> float:
        max_delay = max(float(poll_interval), 0.0)
        if max_delay <= 0:
            return 0.0
        if elapsed < 2:
            return min(max_delay, 0.2)
        if elapsed < 5:
            return min(max_delay, 0.5)
        if elapsed < 15:
            return min(max_delay, 1.0)
        return max_delay

    def _poll_until_complete(self, job_id: str, poll_interval: int,
                             max_wait: int, context: str = "search") -> SearchResult:
        """Poll GET /search/<job_id> until completed, return lazy SearchResult."""
        start_time = time.time()
        while time.time() - start_time < max_wait:
            try:
                resp = self._get_job_status_lightweight(job_id)
                status = resp.get('status')
                logger.info(f"{context} status: {status}")
                if status == 'completed':
                    logger.info(f"{context} completed: {job_id}")
                    return SearchResult({'job_id': job_id, 'status': 'completed'}, self)
                elif status == 'error':
                    raise MalvaAPIError(f"{context} failed: {resp.get('error', 'unknown')}")
            except MalvaAPIError as e:
                if "404" in str(e) or "not found" in str(e).lower():
                    pass  # job not yet visible, keep polling
                else:
                    raise
            delay = self._adaptive_poll_delay(time.time() - start_time, poll_interval)
            if delay > 0:
                time.sleep(delay)
        raise MalvaAPIError(f"{context} timed out after {max_wait} seconds")

    def search(self, query: str,
               wait_for_completion: bool = True,
               poll_interval: int = 2, max_wait: int = 300,
               aggregate_expression: bool = True,
               stranded: Optional[bool] = None,
               min_kmer_presence: int = 0,
               max_kmer_presence: int = 100000) -> SearchResult:
        """
        Search by gene name, natural-language description, or DNA sequence.

        Args:
            query: Gene symbol (``"BRCA1"``), comma-separated list of genes
                (``"BRCA1, TP53"``), natural-language description
                (``"immune cell marker"``), or a DNA sequence string.
                For multiple sequences or sequences >500 nt, prefer
                :meth:`search_sequences`.
            wait_for_completion: Block until the search finishes.  If ``False``,
                returns immediately with a pending :class:`SearchResult` whose
                data is fetched lazily on first ``.df`` access.
            poll_interval: Seconds between status polls (only used when
                ``wait_for_completion=True`` and the job does not complete
                synchronously).
            max_wait: Maximum seconds to wait before raising a timeout error.
            aggregate_expression: If ``True`` (default), returns per-(sample,
                cell-type) aggregate scores suitable for plotting.  If
                ``False``, returns raw per-cell hit arrays. For downstream
                cell-level analysis, use :meth:`retrieve_cells` on the search
                result.
            stranded: Strand mode for k-mer matching.
                ``True``  → forward strand only.
                ``False`` → both strands (reverse-complement included).
                ``None``  (default) → server default (currently forward only).
            min_kmer_presence: Exclude k-mers that appear in fewer than this
                many cells.  Useful to filter out ultra-rare k-mers (noise).
                Default ``0`` (no lower filter).
            max_kmer_presence: Exclude k-mers that appear in more than this
                many cells.  Useful to remove highly repetitive/ubiquitous
                k-mers.  Default ``100000``.

        Returns:
            :class:`SearchResult` with per-(sample, cell-type) expression
            aggregates.  Access ``.df`` to get a pandas DataFrame.

        Example::

            results = client.search("ACTA2")
            print(results.df.head())

            # Both strands, excluding very common k-mers
            results = client.search("GATTACA" * 10,
                                    stranded=False,
                                    max_kmer_presence=50000)
        """
        payload = self._build_search_payload(
            query, wait_for_completion, max_wait, poll_interval,
            aggregate_expression, stranded, min_kmer_presence, max_kmer_presence,
        )
        response = self._request('POST', '/search', json=payload)

        job_id = response.get('job_id')
        if not job_id:
            raise MalvaAPIError("No job_id received from server")

        status = response.get('status')
        logger.info(f"Search submitted: job_id={job_id}, status={status}")

        if status == 'error':
            raise MalvaAPIError(f"Search failed: {response.get('error', 'unknown')}")

        if status == 'completed':
            return SearchResult({
                'job_id': job_id,
                'status': 'completed',
                'aggregate_expression': aggregate_expression,
            }, self)

        if not wait_for_completion:
            return SearchResult({
                'job_id': job_id,
                'status': 'pending',
                'aggregate_expression': aggregate_expression,
            }, self)

        result = self._poll_until_complete(job_id, poll_interval, max_wait, "search")
        result.raw_data['aggregate_expression'] = aggregate_expression
        return result
    
    def print_dict_summary(self, d, max_items=3, indent=0):
        """
        Print a summary of a dictionary structure, truncating large lists.
        
        Args:
            d: The dictionary to summarize
            max_items: Maximum number of items to show from lists (default 3)
            indent: Current indentation level (used for recursion)
        """
        spacing = "  " * indent
        
        for key, value in d.items():
            if isinstance(value, dict):
                print(f"{spacing}{key}: {{")
                self.print_dict_summary(value, max_items, indent + 1)
                print(f"{spacing}}}")
            elif isinstance(value, list):
                if len(value) <= max_items:
                    print(f"{spacing}{key}: {value}")
                else:
                    preview = value[:max_items]
                    print(f"{spacing}{key}: {preview} ... ({len(value)} total items)")
            else:
                print(f"{spacing}{key}: {value}")

    @staticmethod
    def _is_not_found_error(exc: Exception) -> bool:
        msg = str(exc).lower()
        return "endpoint not found" in msg or "not found" in msg or "404" in msg

    def _request_first_available(self, method: str, endpoints: List[str], **kwargs) -> Tuple[Dict[str, Any], str]:
        """Try equivalent endpoint paths for explorer-backed and search-backed URLs."""
        last_error = None
        for endpoint in endpoints:
            try:
                return self._request(method, endpoint, **kwargs), endpoint
            except MalvaAPIError as exc:
                last_error = exc
                if self._is_not_found_error(exc):
                    continue
                raise
        if last_error:
            raise last_error
        raise MalvaAPIError("No endpoint candidates were provided")

    def _request_raw(self, method: str, endpoint: str, **kwargs) -> requests.Response:
        """Make a request and return the raw response for binary endpoints."""
        if endpoint.startswith(('http://', 'https://')):
            url = endpoint
        else:
            url = urljoin(self.base_url, endpoint.lstrip('/'))
        try:
            response = self.session.request(method, url, timeout=self.timeout, **kwargs)
        except requests.exceptions.Timeout:
            raise MalvaAPIError(f"Request timed out after {self.timeout} seconds")
        except requests.exceptions.ConnectionError:
            raise MalvaAPIError(f"Could not connect to {self.base_url}. Check if the server is running.")
        except requests.exceptions.RequestException as exc:
            raise MalvaAPIError(f"Request failed: {exc}")

        logger.debug(f"{method} {url} -> {response.status_code}")
        if response.status_code == 404:
            raise MalvaAPIError(f"Endpoint not found: {endpoint}")
        if response.status_code in [401, 403]:
            raise AuthenticationError("Authentication failed")
        if response.status_code == 429:
            raise QuotaExceededError("Search quota exceeded")
        if not response.ok:
            try:
                error_data = response.json()
                msg = error_data.get('detail', error_data.get('error', f'HTTP {response.status_code}'))
            except (json.JSONDecodeError, ValueError):
                msg = response.text[:200] if response.text else f'HTTP {response.status_code}'
            raise MalvaAPIError(f"API error: {msg}")
        return response

    @staticmethod
    def _normalise_id_list(values: Optional[Union[int, str, List[Union[int, str]]]]) -> Optional[List[int]]:
        if values is None:
            return None
        if isinstance(values, (int, np.integer, str)):
            values = [values]
        return [int(v) for v in values]

    @staticmethod
    def _normalise_feature_list(features: Optional[Union[str, List[str]]]) -> Optional[List[str]]:
        if features is None:
            return None
        if isinstance(features, str):
            return [features]
        return [str(f) for f in features]

    def _resolve_search_job(self, search_job: Union[str, SearchResult, Dict[str, Any]],
                            poll_interval: int,
                            max_wait: int) -> Tuple[str, Optional[SearchResult]]:
        if isinstance(search_job, SearchResult):
            aggregate_expression = (
                search_job.raw_data.get('aggregate_expression')
                if isinstance(search_job.raw_data, dict)
                else None
            )
            if search_job.status in {'pending', 'running'}:
                search_job = self.wait_for_job(
                    search_job.job_id,
                    poll_interval=poll_interval,
                    max_wait=max_wait,
                )
                if aggregate_expression is not None:
                    search_job.raw_data['aggregate_expression'] = aggregate_expression
            return search_job.job_id, search_job
        if isinstance(search_job, dict):
            job_id = search_job.get('job_id')
            if not job_id:
                raise ValueError("search_job dictionary must include job_id")
            return str(job_id), SearchResult(search_job, self)
        if isinstance(search_job, str):
            if not search_job.strip():
                raise ValueError("search_job must be a SearchResult or a non-empty job_id")
            return search_job.strip(), None
        raise TypeError("search_job must be a SearchResult, job_id string, or job response dictionary")

    def _sample_metadata_frame(self, sample_ids: List[int]) -> pd.DataFrame:
        if not sample_ids:
            return pd.DataFrame()
        try:
            metadata = self.get_sample_metadata(sample_ids)
        except Exception as exc:
            logger.warning(f"Sample metadata lookup failed during cell retrieval: {exc}")
            return pd.DataFrame()

        rows = []
        for sample_id, values in metadata.items():
            if isinstance(values, dict):
                row = {'sample_id': int(sample_id)}
                row.update(values)
                rows.append(row)
        return pd.DataFrame(rows) if rows else pd.DataFrame()

    def _cell_barcodes_frame(self, cells: pd.DataFrame, batch_size: int = 100000) -> pd.DataFrame:
        if cells.empty:
            return pd.DataFrame()
        rows = []
        for start in range(0, len(cells), batch_size):
            batch = cells.iloc[start:start + batch_size]
            payload = {
                'sample_ids': batch['sample_id'].astype(int).tolist(),
                'cell_ids': batch['cell_id'].astype(int).tolist(),
            }
            try:
                response = self._request('POST', '/api/cells/barcodes', json=payload)
            except Exception as exc:
                logger.warning(f"Barcode lookup failed during cell retrieval: {exc}")
                return pd.DataFrame()
            barcodes = response.get('barcodes', [])
            for row_index, sample_id, cell_id, barcode in zip(
                batch['row_index'].tolist(),
                payload['sample_ids'],
                payload['cell_ids'],
                barcodes,
            ):
                rows.append({
                    'row_index': int(row_index),
                    'sample_id': int(sample_id),
                    'cell_id': int(cell_id),
                    'barcode': barcode,
                })
        return pd.DataFrame(rows) if rows else pd.DataFrame()

    def _result_from_entry_rows(self, job_id: str,
                                entries: pd.DataFrame,
                                features: pd.DataFrame,
                                source: str,
                                include_barcodes: bool,
                                include_sample_metadata: bool) -> CellExpressionMatrixResult:
        if entries.empty:
            cells = pd.DataFrame(columns=['row_index', 'sample_id', 'cell_id'])
            matrix_entries = pd.DataFrame(columns=['row_index', 'feature_index', 'value'])
        else:
            cells = (
                entries[['sample_id', 'cell_id']]
                .drop_duplicates()
                .sort_values(['sample_id', 'cell_id'])
                .reset_index(drop=True)
            )
            cells.insert(0, 'row_index', np.arange(1, len(cells) + 1, dtype=np.int64))
            matrix_entries = entries.merge(cells, on=['sample_id', 'cell_id'], how='left')
            matrix_entries = matrix_entries[['row_index', 'feature_index', 'value']]
            matrix_entries = matrix_entries.sort_values(['feature_index', 'row_index']).reset_index(drop=True)

        sample_ids = cells['sample_id'].astype(int).drop_duplicates().tolist() if not cells.empty else []
        sample_metadata = self._sample_metadata_frame(sample_ids) if include_sample_metadata else pd.DataFrame()
        normalization = pd.DataFrame()
        barcodes = self._cell_barcodes_frame(cells) if include_barcodes else pd.DataFrame()

        return CellExpressionMatrixResult(
            cells=cells,
            features=features,
            matrix_entries=matrix_entries,
            normalization_factors=normalization,
            sample_metadata=sample_metadata,
            barcodes=barcodes,
            client=self,
            job_id=job_id,
            source=source,
        )

    def _retrieve_cells_from_full_search_endpoint(self, job_id: str,
                                                  features: Optional[List[str]],
                                                  sample_ids: Optional[List[int]],
                                                  include_barcodes: bool,
                                                  include_sample_metadata: bool) -> CellExpressionMatrixResult:
        response = self._request('GET', f'/search/{job_id}', params={'max_cells': 0})
        status = response.get('status') if hasattr(response, 'get') else None
        if status and status != 'completed':
            raise MalvaAPIError(f"Search job {job_id} is {status}")

        feature_filter = set(features) if features else None
        sample_filter = set(sample_ids) if sample_ids else None
        feature_rows = []
        entry_frames = []
        feature_index = 1

        if isinstance(response, dict):
            results = response.get('results', {}) if isinstance(response, dict) else {}
            iterator = (
                (feature, data) for feature, data in results.items()
                if isinstance(feature, str) and not feature.startswith('_') and isinstance(data, dict)
            )
        elif hasattr(response, 'items'):
            iterator = response.items()
        else:
            iterator = iter(())

        for feature, data in iterator:
            if feature_filter is not None and feature not in feature_filter:
                continue

            if hasattr(data, 'cells') and hasattr(data, 'samples'):
                cells = np.asarray(data.cells, dtype=np.int64)
                samples = np.asarray(data.samples, dtype=np.int64)
                values = np.asarray(getattr(data, 'expression', np.ones(len(cells))), dtype=np.float32)
            else:
                cells = np.asarray(data.get('cell', []), dtype=np.int64)
                samples = np.asarray(data.get('sample', []), dtype=np.int64)
                values = np.asarray(data.get('expression', np.ones(len(cells))), dtype=np.float32)

            n = min(len(cells), len(samples), len(values))
            if n == 0:
                continue
            cells = cells[:n]
            samples = samples[:n]
            values = values[:n]

            if sample_filter is not None:
                mask = np.isin(samples, np.fromiter(sample_filter, dtype=np.int64))
                cells = cells[mask]
                samples = samples[mask]
                values = values[mask]
            if len(cells) == 0:
                continue

            feature_rows.append({
                'feature_index': feature_index,
                'job_id': job_id,
                'feature': feature,
                'label': feature,
                'source': 'search_full',
            })
            entry_frames.append(pd.DataFrame({
                'sample_id': samples,
                'cell_id': cells,
                'feature_index': feature_index,
                'value': values,
            }))
            feature_index += 1

        if not entry_frames:
            raise MalvaAPIError("No cells matched the requested feature or sample filter")

        return self._result_from_entry_rows(
            job_id=job_id,
            entries=pd.concat(entry_frames, ignore_index=True),
            features=pd.DataFrame(feature_rows),
            source='search_full',
            include_barcodes=include_barcodes,
            include_sample_metadata=include_sample_metadata,
        )

    def _retrieve_cells_from_cell_ids_endpoint(self, job_id: str,
                                               features: Optional[List[str]],
                                               sample_ids: Optional[List[int]],
                                               include_barcodes: bool,
                                               include_sample_metadata: bool) -> CellExpressionMatrixResult:
        response = self._request('GET', f'/search/{job_id}/cell-ids')
        results = response.get('results', {}) if isinstance(response, dict) else {}
        if not results:
            raise MalvaAPIError("No cell IDs returned for this search job")

        feature_filter = set(features) if features else None
        sample_filter = set(sample_ids) if sample_ids else None
        feature_rows = []
        entry_frames = []
        feature_index = 1

        for feature, data in results.items():
            if feature_filter is not None and feature not in feature_filter:
                continue
            if not isinstance(data, dict):
                continue
            cells = np.asarray(data.get('cell', []), dtype=np.int64)
            samples = np.asarray(data.get('sample', []), dtype=np.int64)
            n = min(len(cells), len(samples))
            if n == 0:
                continue
            cells = cells[:n]
            samples = samples[:n]
            values_raw = data.get('expression')
            if values_raw is not None and len(values_raw) >= n:
                values = np.asarray(values_raw[:n], dtype=np.float32)
            else:
                values = np.ones(n, dtype=np.float32)

            if sample_filter is not None:
                mask = np.isin(samples, np.fromiter(sample_filter, dtype=np.int64))
                cells = cells[mask]
                samples = samples[mask]
                values = values[mask]
            if len(cells) == 0:
                continue

            feature_rows.append({
                'feature_index': feature_index,
                'job_id': job_id,
                'feature': feature,
                'label': feature,
                'source': 'search_cell_ids',
            })
            entry_frames.append(pd.DataFrame({
                'sample_id': samples,
                'cell_id': cells,
                'feature_index': feature_index,
                'value': values,
            }))
            feature_index += 1

        if not entry_frames:
            raise MalvaAPIError("No cells matched the requested feature or sample filter")

        return self._result_from_entry_rows(
            job_id=job_id,
            entries=pd.concat(entry_frames, ignore_index=True),
            features=pd.DataFrame(feature_rows),
            source='search_cell_ids',
            include_barcodes=include_barcodes,
            include_sample_metadata=include_sample_metadata,
        )

    def _retrieve_cells_from_global_ids_endpoint(self, job_id: str,
                                                features: Optional[List[str]],
                                                sample_ids: Optional[List[int]],
                                                include_barcodes: bool,
                                                include_sample_metadata: bool) -> CellExpressionMatrixResult:
        params = {}
        if sample_ids:
            params['sample_ids'] = ','.join(str(s) for s in sample_ids)
        response = self._request_raw('GET', f'/search/{job_id}/global-cell-ids', params=params)
        decoded = msgpack.unpackb(response.content, raw=False, strict_map_key=False)
        cells = np.frombuffer(decoded.get('c', b''), dtype=np.int64)
        samples = np.frombuffer(decoded.get('s', b''), dtype=np.int64)
        n = min(len(cells), len(samples))
        if n == 0:
            raise MalvaAPIError("No cell IDs returned for this search job")
        cells = cells[:n]
        samples = samples[:n]
        if sample_ids:
            mask = np.isin(samples, np.asarray(sample_ids, dtype=np.int64))
            cells = cells[mask]
            samples = samples[mask]
        if len(cells) == 0:
            raise MalvaAPIError("No cells matched the requested sample filter")

        label = features[0] if features and len(features) == 1 else 'positive_cells'
        entries = pd.DataFrame({
            'sample_id': samples,
            'cell_id': cells,
            'feature_index': 1,
            'value': np.ones(len(cells), dtype=np.float32),
        })
        feature_df = pd.DataFrame([{
            'feature_index': 1,
            'job_id': job_id,
            'feature': label,
            'label': label,
            'source': 'global_cell_ids',
        }])
        return self._result_from_entry_rows(
            job_id=job_id,
            entries=entries,
            features=feature_df,
            source='global_cell_ids',
            include_barcodes=include_barcodes,
            include_sample_metadata=include_sample_metadata,
        )

    def _retrieve_cells_from_explorer_endpoint(self, job_id: str,
                                               features: Optional[List[str]],
                                               sample_ids: Optional[List[int]],
                                               include_barcodes: bool,
                                               include_sample_metadata: bool) -> CellExpressionMatrixResult:
        params = {}
        if sample_ids:
            params['sample_ids'] = ','.join(str(s) for s in sample_ids)
        if include_barcodes:
            params['barcodes'] = '1'

        response = self._request('GET', f'/api/expression/results/{job_id}/cells', params=params)
        samples_payload = response.get('samples', {}) if isinstance(response, dict) else {}
        entries = []
        barcode_rows = []
        for sample_id, payload in samples_payload.items():
            if not isinstance(payload, dict):
                continue
            cell_ids = payload.get('cell_ids') or []
            barcodes = payload.get('barcodes') or []
            for idx, cell_id in enumerate(cell_ids):
                entries.append({
                    'sample_id': int(sample_id),
                    'cell_id': int(cell_id),
                    'feature_index': 1,
                    'value': 1.0,
                })
                if include_barcodes and idx < len(barcodes):
                    barcode_rows.append({
                        'sample_id': int(sample_id),
                        'cell_id': int(cell_id),
                        'barcode': barcodes[idx],
                    })

        if not entries:
            raise MalvaAPIError("No cells returned for this search job")

        label = features[0] if features and len(features) == 1 else 'positive_cells'
        result = self._result_from_entry_rows(
            job_id=job_id,
            entries=pd.DataFrame(entries),
            features=pd.DataFrame([{
                'feature_index': 1,
                'job_id': job_id,
                'feature': label,
                'label': label,
                'source': 'expression_cells',
            }]),
            source='expression_cells',
            include_barcodes=False,
            include_sample_metadata=include_sample_metadata,
        )
        if barcode_rows:
            barcode_df = pd.DataFrame(barcode_rows)
            result._barcodes_df = result.cells.merge(
                barcode_df,
                on=['sample_id', 'cell_id'],
                how='inner',
            )
        return result

    def retrieve_cells(self, search_job: Union[str, SearchResult, Dict[str, Any]],
                       features: Optional[Union[str, List[str]]] = None,
                       sample_ids: Optional[Union[int, str, List[Union[int, str]]]] = None,
                       include_barcodes: bool = False,
                       include_sample_metadata: bool = True,
                       poll_interval: int = 2,
                       max_wait: int = 300) -> CellExpressionMatrixResult:
        """
        Retrieve positive cells for an existing search job.

        Run :meth:`search` or :meth:`search_sequences` first, then pass the
        returned SearchResult or its job_id here. This method does not start a
        new search and does not use the ZIP export workflow. Aggregate jobs use
        the fast binary cell-ID endpoint when feature-level values are not
        present. Jobs created with ``aggregate_expression=False`` can provide
        per-feature values when the search service exposes them.

        """
        job_id, result = self._resolve_search_job(search_job, poll_interval, max_wait)
        feature_list = self._normalise_feature_list(features)
        if feature_list is None and result is not None:
            raw_results = result.raw_data.get('results', {}) if isinstance(result.raw_data, dict) else {}
            feature_list = [
                key for key, value in raw_results.items()
                if isinstance(key, str) and not key.startswith('_') and isinstance(value, dict)
            ] or None
        sample_id_list = self._normalise_id_list(sample_ids)

        errors = []
        aggregate_job = None
        if result is not None and isinstance(result.raw_data, dict):
            aggregate_job = result.raw_data.get('aggregate_expression')

        if aggregate_job is True and feature_list is None:
            fetchers = (
                self._retrieve_cells_from_global_ids_endpoint,
                self._retrieve_cells_from_explorer_endpoint,
                self._retrieve_cells_from_cell_ids_endpoint,
                self._retrieve_cells_from_full_search_endpoint,
            )
        else:
            fetchers = (
                self._retrieve_cells_from_full_search_endpoint,
                self._retrieve_cells_from_cell_ids_endpoint,
                self._retrieve_cells_from_global_ids_endpoint,
                self._retrieve_cells_from_explorer_endpoint,
            )

        for fetcher in fetchers:
            try:
                return fetcher(
                    job_id,
                    feature_list,
                    sample_id_list,
                    include_barcodes,
                    include_sample_metadata,
                )
            except MalvaAPIError as exc:
                message = str(exc)
                errors.append(message)
                recoverable = (
                    self._is_not_found_error(exc)
                    or "No cell IDs returned" in message
                    or "No cells returned" in message
                    or "No cells matched" in message
                )
                if not recoverable:
                    raise

        detail = '; '.join(errors[-3:]) if errors else 'No retrieval endpoint was available'
        raise MalvaAPIError(f"Could not retrieve cells for search job {job_id}: {detail}")
    
    def search_sequences(self, sequences: Union[str, List[str]],
                         wait_for_completion: bool = True,
                         poll_interval: int = 2, max_wait: int = 300,
                         aggregate_expression: bool = True,
                         stranded: Optional[bool] = None,
                         min_kmer_presence: int = 0,
                         max_kmer_presence: int = 100000) -> SearchResult:
        """
        Search for one or more DNA sequences.

        Accepts either a single sequence string or a list of sequences.
        Multiple sequences are submitted as a single FASTA batch request
        (only one search quota slot consumed).

        Args:
            sequences: A single DNA sequence string **or** a list of DNA
                sequence strings.

                - Single string: up to 500,000 nt.
                - List: each entry up to 500,000 nt; total across all
                  sequences must not exceed 100,000 nt.

                Valid nucleotide characters: ``A``, ``T``, ``G``, ``C``,
                ``U``, ``N`` (case-insensitive).
            wait_for_completion: Block until the search finishes.
            poll_interval: Seconds between status polls.
            max_wait: Maximum seconds to wait before raising a timeout error.
            aggregate_expression: If ``True`` (default), returns per-(sample,
                cell-type) aggregate scores. For downstream cell-level
                analysis, use :meth:`retrieve_cells` on the search result.
            stranded: Strand mode for k-mer matching.
                ``True``  → forward strand only.
                ``False`` → both strands (reverse-complement included).
                ``None``  (default) → server default (currently forward only).
            min_kmer_presence: Exclude k-mers that appear in fewer than this
                many cells.  Default ``0`` (no lower filter).
            max_kmer_presence: Exclude k-mers that appear in more than this
                many cells.  Default ``100000``.

        Returns:
            :class:`SearchResult` with per-(sample, cell-type) aggregates.

        Examples::

            # Single sequence
            results = client.search_sequences("ATGCATGCATGC")

            # Multiple sequences (FASTA batch)
            results = client.search_sequences([
                "ATGCATGCATGC",
                "GCTAGCTAGCTA",
            ])

            # Both strands, filter common k-mers
            results = client.search_sequences("GATTACA" * 20,
                                              stranded=False,
                                              max_kmer_presence=50000)
        """
        valid_nucleotides = set('ATGCUN')

        # ── Single sequence ────────────────────────────────────────────────
        if isinstance(sequences, str):
            sequence = sequences
            if not sequence:
                raise ValueError("Sequence must not be empty")
            if len(sequence) > 500_000:
                raise ValueError("Sequence length cannot exceed 500,000 nucleotides")
            if not all(c.upper() in valid_nucleotides for c in sequence):
                raise ValueError(
                    "Sequence contains invalid nucleotides. "
                    "Only A, T, G, C, U, N are allowed."
                )
            return self.search(
                sequence,
                wait_for_completion=wait_for_completion,
                poll_interval=poll_interval,
                max_wait=max_wait,
                aggregate_expression=aggregate_expression,
                stranded=stranded,
                min_kmer_presence=min_kmer_presence,
                max_kmer_presence=max_kmer_presence,
            )

        # ── Batch (list of sequences) ──────────────────────────────────────
        if not sequences:
            raise ValueError("At least one sequence is required")

        total_nt = 0
        for i, seq in enumerate(sequences):
            if not seq:
                raise ValueError(f"Sequence {i + 1} is empty")
            if not all(c.upper() in valid_nucleotides for c in seq):
                raise ValueError(
                    f"Sequence {i + 1} contains invalid nucleotides. "
                    "Only A, T, G, C, U, N are allowed."
                )
            total_nt += len(seq)

        if total_nt > 100_000:
            raise ValueError(
                f"Total nucleotides ({total_nt:,}) exceed the 100,000 nt batch limit. "
                "Submit fewer or shorter sequences, or call search_sequences() "
                "separately for each sequence."
            )

        fasta_lines = []
        for i, seq in enumerate(sequences, 1):
            fasta_lines.append(f">seq_{i}")
            fasta_lines.append(seq)
        fasta_content = "\n".join(fasta_lines) + "\n"

        payload = self._build_search_payload(
            'batch sequence search', wait_for_completion, max_wait, poll_interval,
            aggregate_expression, stranded, min_kmer_presence, max_kmer_presence,
        )
        payload['fasta_sequences'] = fasta_content

        logger.info(
            f"Submitting batch search: {len(sequences)} sequences, "
            f"{total_nt:,} total nucleotides"
        )

        response = self._request('POST', '/search', json=payload)

        job_id = response.get('job_id')
        if not job_id:
            raise MalvaAPIError("No job_id received from server")

        status = response.get('status')
        logger.info(f"Batch search submitted: job_id={job_id}, status={status}")

        if status == 'error':
            raise MalvaAPIError(f"Batch search failed: {response.get('error', 'unknown')}")

        if status == 'completed':
            return SearchResult({
                'job_id': job_id,
                'status': 'completed',
                'aggregate_expression': aggregate_expression,
            }, self)

        if not wait_for_completion:
            return SearchResult({
                'job_id': job_id,
                'status': 'pending',
                'aggregate_expression': aggregate_expression,
            }, self)

        result = self._poll_until_complete(job_id, poll_interval, max_wait, "batch search")
        result.raw_data['aggregate_expression'] = aggregate_expression
        return result

    def search_genes(self, genes: List[str],
                     wait_for_completion: bool = True,
                     poll_interval: int = 2, max_wait: int = 300,
                     aggregate_expression: bool = True,
                     stranded: Optional[bool] = None,
                     min_kmer_presence: int = 0,
                     max_kmer_presence: int = 100000) -> SearchResult:
        """
        Search for multiple gene symbols in a single request.

        Args:
            genes: List of gene symbols (e.g., ``["BRCA1", "TP53"]``).
                Maximum 10 genes per request.
            wait_for_completion: Block until finished.
            poll_interval: Seconds between status polls.
            max_wait: Maximum seconds to wait before raising a timeout error.
            aggregate_expression: If ``True`` (default), returns per-(sample,
                cell-type) aggregate scores.
            stranded: Strand mode for k-mer matching.
                ``True``  → forward strand only.
                ``False`` → both strands (reverse-complement included).
                ``None``  (default) → server default.
            min_kmer_presence: Exclude k-mers appearing in fewer than this
                many cells.  Default ``0``.
            max_kmer_presence: Exclude k-mers appearing in more than this
                many cells.  Default ``100000``.

        Returns:
            :class:`SearchResult`.
        """
        if not genes:
            raise ValueError("At least one gene is required")
        if len(genes) > 10:
            raise ValueError(f"Too many genes ({len(genes)}). Maximum is 10 per request.")

        gene_query = ", ".join(genes)
        logger.info(f"Submitting gene search: {gene_query}")

        payload = self._build_search_payload(
            gene_query, wait_for_completion, max_wait, poll_interval,
            aggregate_expression, stranded, min_kmer_presence, max_kmer_presence,
        )
        response = self._request('POST', '/search', json=payload)

        job_id = response.get('job_id')
        if not job_id:
            raise MalvaAPIError("No job_id received from server")

        status = response.get('status')
        logger.info(f"Gene search submitted: job_id={job_id}, status={status}")

        if status == 'error':
            raise MalvaAPIError(f"Gene search failed: {response.get('error', 'unknown')}")

        if status == 'completed':
            return SearchResult({
                'job_id': job_id,
                'status': 'completed',
                'aggregate_expression': aggregate_expression,
            }, self)

        if not wait_for_completion:
            return SearchResult({
                'job_id': job_id,
                'status': 'pending',
                'aggregate_expression': aggregate_expression,
            }, self)

        result = self._poll_until_complete(job_id, poll_interval, max_wait, "gene search")
        result.raw_data['aggregate_expression'] = aggregate_expression
        return result
    
    def get_samples(self, page: int = 1, page_size: int = 50, 
                   filters: Dict[str, Any] = None, search_query: str = "") -> Dict[str, Any]:
        """
        Retrieve samples with optional filtering
        
        Args:
            page: Page number (1-based)
            page_size: Number of samples per page
            filters: Dictionary of filters (e.g., {'organ': 'brain', 'species': 'human'})
            search_query: Text search query
            
        Returns:
            Dictionary containing samples and pagination info
        """
        params = {
            'page': page,
            'page_size': page_size,
            'q': search_query
        }
        
        if filters:
            for key, value in filters.items():
                if isinstance(value, list):
                    params[key] = value
                else:
                    params[key] = [value] if value else []
        
        return self._request('GET', '/samples/api/search', params=params)
    
    def search_samples(self, query: str = "", organ: List[str] = None, 
                      species: List[str] = None, study: List[str] = None,
                      disease: List[str] = None, **kwargs) -> pd.DataFrame:
        """
        Search samples and return as a DataFrame
        
        Args:
            query: Text search query
            organ: List of organs to filter by
            species: List of species to filter by
            study: List of studies to filter by
            disease: List of diseases to filter by
            **kwargs: Additional filters
            
        Returns:
            DataFrame containing sample metadata
        """
        filters = {}
        if organ:
            filters['organ'] = organ
        if species:
            filters['species'] = species
        if study:
            filters['study'] = study
        if disease:
            filters['disease'] = disease
        
        # Add any additional filters
        filters.update(kwargs)
        
        # Get all pages
        all_samples = []
        page = 1
        
        while True:
            response = self.get_samples(
                page=page, 
                page_size=100,
                filters=filters,
                search_query=query
            )
            
            samples = response.get('samples', [])
            if not samples:
                break
                
            all_samples.extend(samples)
            
            if not response.get('has_more', False):
                break
                
            page += 1
        
        return pd.DataFrame(all_samples)
    
    def get_studies(self) -> pd.DataFrame:
        """
        Get all studies with summary information
        
        Returns:
            DataFrame containing study information
        """
        response = self._request('GET', '/samples/api/studies')
        studies = response.get('studies', [])
        return pd.DataFrame(studies)

    def _normalise_sample_download_response(self, sample_uuid: str, response: Dict[str, Any]) -> Dict[str, Any]:
        """Normalize explorer and metadata-service sample download responses."""
        download_url = response.get('download_url')
        token = response.get('download_token')
        if not download_url and token:
            download_url = f'/api/download/{token}'
        if not download_url:
            raise MalvaAPIError("Download URL not provided by server")
        parsed_download = urlparse(download_url)
        if parsed_download.scheme and parsed_download.netloc:
            if parsed_download.path.startswith('/api/download/'):
                parsed_base = urlparse(self.base_url)
                download_url = urlunparse((
                    parsed_base.scheme,
                    parsed_base.netloc,
                    parsed_download.path,
                    '',
                    parsed_download.query,
                    parsed_download.fragment,
                ))
        else:
            download_url = urljoin(self.base_url, download_url.lstrip('/'))

        expires_in = response.get('expires_in_seconds')
        if expires_in is None and response.get('expires_at') is not None:
            expires_in = max(0, int(float(response['expires_at']) - time.time()))
        if expires_in is None:
            expires_in = 0

        return {
            'download_url': download_url,
            'filename': response.get('filename', f'{sample_uuid}.h5ad'),
            'file_size': response.get('file_size', 0),
            'expires_in_seconds': expires_in,
        }

    def _request_sample_download(self, sample_uuid: str) -> Dict[str, Any]:
        """Request sample download metadata from either compatibility API shape."""
        try:
            response = self._request('GET', f'/api/samples/{sample_uuid}/download')
        except MalvaAPIError as exc:
            if not self._is_not_found_error(exc):
                raise
            response = self._request('POST', f'/api/samples/{sample_uuid}/request-download')
        return self._normalise_sample_download_response(sample_uuid, response)
    
    def download_sample(self, sample_uuid: str, output_path: Optional[str] = None, 
                    chunk_size: int = 8192, show_progress: bool = True) -> Union[str, 'ad.AnnData']:
        """
        Download a sample's data file with enhanced error handling
        
        Args:
            sample_uuid: UUID of the sample to download
            output_path: Path to save the file. If None, loads into memory as AnnData
            chunk_size: Size of chunks to download (bytes)
            show_progress: Whether to show download progress
            
        Returns:
            File path if saved to disk, or AnnData object if loaded in memory
        """
        try:
            response = self._request_sample_download(sample_uuid)
            download_url = response['download_url']
            filename = response.get('filename', f'{sample_uuid}.h5ad')
            file_size = response.get('file_size', 0)
            expires_in = response.get('expires_in_seconds', 0)
            
            if expires_in <= 0:
                raise MalvaAPIError("Download token has expired")
            
            logger.info(f"Downloading {filename} ({file_size:,} bytes)")
            
            # Download the file with progress tracking
            file_response = self.session.get(download_url, stream=True, timeout=self.timeout)
            
            # Handle download errors
            if file_response.status_code == 404:
                raise MalvaAPIError("Download token expired or file not found")
            elif file_response.status_code == 403:
                raise MalvaAPIError("Access denied to file")
            file_response.raise_for_status()
            
            if output_path:
                # Save to disk
                output_path = Path(output_path)
                output_path.parent.mkdir(parents=True, exist_ok=True)
                
                downloaded_bytes = 0
                with open(output_path, 'wb') as f:
                    for chunk in file_response.iter_content(chunk_size=chunk_size):
                        if chunk:
                            f.write(chunk)
                            downloaded_bytes += len(chunk)
                            
                            if show_progress and file_size > 0:
                                progress = (downloaded_bytes / file_size) * 100
                                print(f"\rDownloading: {progress:.1f}% ({downloaded_bytes:,}/{file_size:,} bytes)", end='')
                
                if show_progress:
                    print()  # New line after progress
                
                logger.info(f"Downloaded sample {sample_uuid} to {output_path}")
                return str(output_path)
            else:
                # Load into memory
                try:
                    import tempfile
                    with tempfile.NamedTemporaryFile(suffix='.h5ad', delete=False) as tmp_file:
                        downloaded_bytes = 0
                        for chunk in file_response.iter_content(chunk_size=chunk_size):
                            if chunk:
                                tmp_file.write(chunk)
                                downloaded_bytes += len(chunk)
                                
                                if show_progress and file_size > 0:
                                    progress = (downloaded_bytes / file_size) * 100
                                    print(f"\rDownloading: {progress:.1f}% ({downloaded_bytes:,}/{file_size:,} bytes)", end='')
                        
                        tmp_path = tmp_file.name
                    
                    if show_progress:
                        print()  # New line after progress
                    
                    # Load with anndata
                    import anndata as ad
                    adata = ad.read_h5ad(tmp_path)
                    
                    # Clean up temp file
                    Path(tmp_path).unlink()
                    
                    logger.info(f"Loaded sample {sample_uuid} into memory")
                    return adata
                    
                except ImportError:
                    raise ImportError("anndata is required to load H5AD files into memory. "
                                    "Install with: pip install anndata")
                    
        except requests.exceptions.RequestException as e:
            if "404" in str(e):
                raise MalvaAPIError(f"Sample {sample_uuid} not found")
            else:
                raise MalvaAPIError(f"Download failed: {e}")
        except Exception as e:
            raise MalvaAPIError(f"Download failed: {e}")

    def check_sample_availability(self, sample_uuid: str) -> Dict[str, Any]:
        """
        Check if a sample is available for download without creating a download token
        
        Args:
            sample_uuid: UUID of the sample to check
            
        Returns:
            Dictionary with availability information
        """
        try:
            response = self._request('GET', f'/api/samples/{sample_uuid}/download/info')
            return response
        except MalvaAPIError as e:
            if "404" in str(e) or "not found" in str(e).lower():
                return {
                    'sample_uuid': sample_uuid,
                    'downloadable': False,
                    'error': 'Sample not found'
                }
            raise
    
    def get_sample_metadata(self, sample_ids: List[int]) -> Dict[int, Dict[str, Any]]:
        """
        Get metadata for specific samples
        
        Args:
            sample_ids: List of sample IDs
            
        Returns:
            Dictionary mapping sample IDs to their metadata
        """
        data = {'ids': sample_ids}
        response = self._request('POST', '/api/sample-metadata', json=data)
        
        # Convert string keys back to integers if needed
        metadata = response.get('metadata', {})
        converted_metadata = {}
        
        for key, value in metadata.items():
            try:
                # Try to convert key to int
                int_key = int(key)
                converted_metadata[int_key] = value
            except (ValueError, TypeError):
                # If conversion fails, keep original key
                converted_metadata[key] = value
        
        return converted_metadata
    
    def get_cells_by_metadata(self,
                              sample_ids: Optional[Union[int, List[int]]] = None,
                              cell_types: Optional[Union[str, List[str]]] = None,
                              include_cell_metadata: bool = True,
                              include_sample_metadata: bool = False,
                              include_all_database_cells: bool = False,
                              poll_interval: int = 2,
                              max_wait: int = 900) -> pd.DataFrame:
        """Fetch all cells matching metadata filters, independent of search jobs.

        This is the denominator table for downstream calculations such as
        fraction expressing or mean expression including zeros. Call it once for
        a sample/cell-type selection and reuse it with any number of search
        results. Database-wide retrieval requires the explicit
        ``include_all_database_cells=True`` opt-in because it can return tens
        of millions of rows.
        """
        sample_id_list = self._normalise_id_list(sample_ids)
        if isinstance(cell_types, str):
            cell_type_list = [cell_types]
        elif cell_types is None:
            cell_type_list = []
        else:
            cell_type_list = [str(ct) for ct in cell_types]
        if not sample_id_list and not cell_type_list and not include_all_database_cells:
            raise ValueError("Provide sample_ids/cell_types or set include_all_database_cells=True")

        request_payload = {
            'sample_ids': sample_id_list or [],
            'cell_types': cell_type_list,
            'include_cell_metadata': include_cell_metadata,
            'include_all_database_cells': include_all_database_cells,
        }
        if include_all_database_cells and not sample_id_list and not cell_type_list:
            queued = self._request('POST', '/exports/all-cells', json=request_payload)
            export_id = queued.get('export_id')
            if not export_id:
                raise MalvaAPIError("No export_id received from server")
            start_time = time.time()
            status_payload = queued
            while time.time() - start_time < max_wait:
                status_payload = self._request('GET', f'/exports/all-cells/{export_id}')
                status = status_payload.get('status')
                if status == 'completed':
                    break
                if status == 'error':
                    raise MalvaAPIError(f"All-cells export failed: {status_payload.get('error', 'unknown')}")
                delay = self._adaptive_poll_delay(time.time() - start_time, poll_interval)
                if delay > 0:
                    time.sleep(delay)
            else:
                raise MalvaAPIError(f"All-cells export timed out after {max_wait} seconds")
            response = self._request_raw('GET', f'/exports/all-cells/{export_id}/download')
        else:
            response = self._request_raw('POST', '/api/cells/by-filter', json=request_payload)
        decoded = msgpack.unpackb(response.content, raw=False, strict_map_key=False)
        cells = np.frombuffer(decoded.get('c', b''), dtype=np.int64)
        samples = np.frombuffer(decoded.get('s', b''), dtype=np.int64)
        n = min(len(cells), len(samples))
        df = pd.DataFrame({'sample_id': samples[:n], 'cell_id': cells[:n]})

        if include_cell_metadata and n:
            dense_valid = bool(decoded.get('_vi_dense', False))
            if decoded.get('_bin_ct_i16'):
                ct_ids = np.frombuffer(decoded.get('_bin_ct_i16', b''), dtype=np.int16)
            else:
                ct_ids = np.frombuffer(decoded.get('_bin_ct', b''), dtype=np.int32)
            if decoded.get('_bin_tc_u16'):
                total_counts_u16 = np.frombuffer(decoded.get('_bin_tc_u16', b''), dtype=np.uint16)
                total_counts = None
            else:
                total_counts = np.frombuffer(decoded.get('_bin_tc', b''), dtype=np.float32)
                total_counts_u16 = None
            valid_all = None if dense_valid else np.frombuffer(decoded.get('_bin_vi', b''), dtype=np.int64)
            mapping = {int(k): v for k, v in (decoded.get('cell_type_mapping') or {}).items()}

            m = min(n, len(ct_ids) if dense_valid else len(valid_all), len(ct_ids))
            if m:
                if dense_valid:
                    valid_ct = None
                    ids = ct_ids[:m]
                else:
                    valid_ct = valid_all[:m]
                    ids = ct_ids[:m]
                    in_bounds = (valid_ct >= 0) & (valid_ct < n)
                    valid_ct = valid_ct[in_bounds]
                    ids = ids[in_bounds]

                sorted_ids = sorted(mapping)
                categories = [mapping[i] for i in sorted_ids]
                codes = np.full(n, -1, dtype=np.int16 if len(sorted_ids) <= np.iinfo(np.int16).max else np.int32)
                if sorted_ids and len(ids):
                    max_id = max(max(sorted_ids), int(ids.max(initial=0)))
                    if max_id <= 1_000_000:
                        remap = np.full(max_id + 1, -1, dtype=codes.dtype)
                        remap[np.asarray(sorted_ids, dtype=np.int64)] = np.arange(len(sorted_ids), dtype=codes.dtype)
                        known = ids >= 0
                        known &= ids <= max_id
                        if dense_valid:
                            codes[np.nonzero(known)[0]] = remap[ids[known]]
                        else:
                            codes[valid_ct[known]] = remap[ids[known]]
                    else:
                        id_to_code = {ct_id: i for i, ct_id in enumerate(sorted_ids)}
                        mapped = np.fromiter(
                            (id_to_code.get(int(x), -1) for x in ids),
                            dtype=codes.dtype,
                            count=len(ids),
                        )
                        if dense_valid:
                            codes[:len(mapped)] = mapped
                        else:
                            codes[valid_ct] = mapped
                df['cell_type'] = pd.Categorical.from_codes(codes, categories=categories)

            counts_u16 = np.zeros(n, dtype=np.uint16)
            if total_counts_u16 is not None:
                mt = min(n, len(total_counts_u16) if dense_valid else len(valid_all), len(total_counts_u16))
                if mt:
                    if dense_valid:
                        counts_u16[:mt] = total_counts_u16[:mt]
                    else:
                        valid_tc = valid_all[:mt]
                        in_bounds = (valid_tc >= 0) & (valid_tc < n)
                        counts_u16[valid_tc[in_bounds]] = total_counts_u16[:mt][in_bounds]
            else:
                mt = min(n, len(total_counts) if dense_valid else len(valid_all), len(total_counts))
                if mt:
                    if dense_valid:
                        counts = np.nan_to_num(total_counts[:mt], nan=0.0, posinf=65535.0, neginf=0.0)
                        counts_u16[:mt] = np.rint(np.clip(counts, 0, 65535)).astype(np.uint16)
                    else:
                        valid_tc = valid_all[:mt]
                        in_bounds = (valid_tc >= 0) & (valid_tc < n)
                        valid_tc = valid_tc[in_bounds]
                        counts = np.nan_to_num(total_counts[:mt][in_bounds], nan=0.0, posinf=65535.0, neginf=0.0)
                        counts_u16[valid_tc] = np.rint(np.clip(counts, 0, 65535)).astype(np.uint16)
            df['total_counts'] = counts_u16

        if include_sample_metadata and n:
            meta = self._sample_metadata_frame(df['sample_id'].drop_duplicates().astype(int).tolist())
            if not meta.empty:
                df = df.merge(meta, on='sample_id', how='left')
        return df

    def get_database_stats(self) -> Dict[str, Any]:
        """Get database statistics"""
        return self._request('GET', '/api/stats')
    
    def get_available_filters(self) -> Dict[str, Any]:
        """Get available filter options for samples"""
        return self._request('GET', '/samples/api/filters')

    def submit_search(self, query: str, aggregate_expression: bool = True,
                      window_size: Optional[int] = None,
                      threshold: Optional[float] = None,
                      stranded: Optional[bool] = None) -> str:
        """
        Submit search without waiting for completion

        Args:
            query: Search query
            aggregate_expression: Whether to aggregate expression by cell type
            window_size: Sliding window size in k-mers
            threshold: Match threshold between 0.0 and 1.0
            stranded: If True, restrict search to the forward strand only

        Returns:
            Job ID for tracking
        """
        if window_size is not None and window_size < 1:
            raise ValueError("window_size must be a positive integer")
        if threshold is not None and not (0.0 <= threshold <= 1.0):
            raise ValueError("threshold must be between 0.0 and 1.0")

        from malva_client.storage import save_search, save_job_status

        data = {
            'query': query,
            'wait_for_completion': False,
            'aggregate_expression': aggregate_expression
        }
        if window_size is not None:
            data['window_size'] = window_size
        if threshold is not None:
            data['threshold'] = threshold
        if stranded is not None:
            data['stranded'] = stranded

        response = self._request('POST', '/search', json=data)
        job_id = response['job_id']

        # Save to local storage
        search_id = save_search(
            query=query,
            job_id=job_id,
            server_url=self.base_url,
            status='submitted'
        )

        save_job_status(
            job_id=job_id,
            status='pending',
            query=query,
            server_url=self.base_url
        )

        logger.info(f"Search submitted with job ID: {job_id}")
        return job_id

    def get_job_status(self, job_id: str) -> Dict[str, Any]:
        """
        Get status of an async job

        Args:
            job_id: Job ID to check

        Returns:
            Dictionary with job status information
        """
        from malva_client.storage import save_job_status

        try:
            response = self._get_job_status_lightweight(job_id)

            # Update local storage
            save_job_status(
                job_id=job_id,
                status=response.get('status', 'unknown'),
                server_url=self.base_url,
                results_data=None,
                error_message=response.get('error') if response.get('status') == 'error' else None
            )

            return response

        except MalvaAPIError as e:
            if "404" in str(e) or "not found" in str(e).lower():
                # Job not found on server, check local storage
                from malva_client.storage import get_job_status as get_local_status
                local_status = get_local_status(job_id)
                if local_status:
                    return {
                        'job_id': job_id,
                        'status': local_status['status'],
                        'error': 'Job not found on server (may have expired)'
                    }
                else:
                    return {
                        'job_id': job_id,
                        'status': 'not_found',
                        'error': 'Job not found'
                    }
            raise

    def get_job_results(self, job_id: str) -> 'SearchResult':
        """
        Get results of a completed job

        Args:
            job_id: Job ID

        Returns:
            SearchResult object
        """
        from malva_client.storage import get_job_status as get_local_status, save_job_status

        # Check local cache — always use lazy loading via job_id
        local_status = get_local_status(job_id)
        if local_status and local_status.get('status') == 'completed':
            logger.info(f"Job {job_id} found in local cache, using lazy load")
            return SearchResult({'job_id': job_id, 'status': 'completed'}, self)

        # Otherwise fetch from server
        status_response = self.get_job_status(job_id)

        if status_response['status'] == 'completed':
            # Use lazy loading — SearchResult will fetch /api/expression-data/{job_id}/results
            save_job_status(
                job_id=job_id,
                status='completed',
                server_url=self.base_url,
                results_data={}
            )
            return SearchResult({'job_id': job_id, 'status': 'completed'}, self)

        elif status_response['status'] == 'error':
            error_msg = status_response.get('error', 'Unknown error')
            raise SearchError(f"Job {job_id} failed: {error_msg}")

        elif status_response['status'] in ['pending', 'running']:
            raise SearchError(f"Job {job_id} is still {status_response['status']}. Use get_job_status() to check progress.")

        else:
            raise SearchError(f"Job {job_id} has unexpected status: {status_response['status']}")

    def list_jobs(self, limit: int = 20) -> List[Dict[str, Any]]:
        """
        List recent jobs

        Args:
            limit: Maximum number of jobs to return

        Returns:
            List of job information
        """
        from malva_client.storage import get_storage

        storage = get_storage()
        history = storage.get_search_history(limit=limit, server_url=self.base_url)

        jobs = []
        for entry in history:
            if entry.get('job_id'):
                jobs.append({
                    'job_id': entry['job_id'],
                    'query': entry['query'],
                    'status': entry['status'],
                    'timestamp': entry['timestamp'],
                    'results_count': entry.get('results_count', 0),
                    'error': entry.get('error_message')
                })

        return jobs

    def wait_for_job(self, job_id: str, poll_interval: int = 2,
                    max_wait: int = 300) -> 'SearchResult':
        """
        Wait for a job to complete and return results

        Args:
            job_id: Job ID to wait for
            poll_interval: How often to check status (seconds)
            max_wait: Maximum time to wait (seconds)

        Returns:
            SearchResult object
        """
        start_time = time.time()

        while time.time() - start_time < max_wait:
            elapsed = time.time() - start_time
            status_response = self.get_job_status(job_id)
            status = status_response.get('status')

            if status == 'completed':
                return self.get_job_results(job_id)
            elif status == 'error':
                error_msg = status_response.get('error', 'Unknown error')
                raise SearchError(f"Job {job_id} failed: {error_msg}")
            elif status in ['pending', 'running']:
                logger.info(f"Job {job_id} status: {status}")
                delay = self._adaptive_poll_delay(elapsed, poll_interval)
                if delay > 0:
                    time.sleep(delay)
            else:
                raise SearchError(f"Job {job_id} has unexpected status: {status}")

        raise MalvaAPIError(f"Job {job_id} timed out after {max_wait} seconds")

    def cancel_job(self, job_id: str) -> bool:
        """
        Attempt to cancel a running job

        Args:
            job_id: Job ID to cancel

        Returns:
            True if cancellation was successful
        """
        try:
            response = self._request('DELETE', f'/search/{job_id}')

            # Update local storage
            from malva_client.storage import save_job_status
            save_job_status(
                job_id=job_id,
                status='cancelled',
                server_url=self.base_url
            )

            return response.get('success', True)

        except MalvaAPIError as e:
            if "404" in str(e):
                logger.warning(f"Job {job_id} not found for cancellation")
                return False
            raise

    # ── Coverage Analysis Methods ──

    def get_coverage(self, chr: str, start: int, end: int,
                     strand: str = '+', zoom: int = 1,
                     metadata_filters: Optional[Dict[str, Any]] = None,
                     poll_interval: int = 2, max_wait: int = 300) -> 'CoverageResult':
        """
        Get genomic region coverage across cell types.

        Args:
            chr: Chromosome name (e.g., 'chr1')
            start: Start position
            end: End position
            strand: Strand ('+', '-', or 'both')
            zoom: Zoom level for binning
            metadata_filters: Optional metadata filters (e.g., {'organ': ['brain']})
            poll_interval: How often to check for completion (seconds)
            max_wait: Maximum time to wait for completion (seconds)

        Returns:
            CoverageResult object
        """
        # Map +/- notation to server's expected format
        strand_map = {'+': 'forward', '-': 'reverse', 'both': 'both',
                      'forward': 'forward', 'reverse': 'reverse'}
        data = {
            'chromosome': chr,
            'start': start,
            'end': end,
            'strand': strand_map.get(strand, 'both'),
            'zoom_level': zoom,
        }
        if metadata_filters:
            data['metadata_filters'] = metadata_filters

        response = self._request('POST', '/api/genome-browser/search', json=data)
        job_id = response.get('job_id')

        if not job_id:
            raise MalvaAPIError("No job_id received from coverage search")

        region = f"{chr}:{start}-{end}({strand})"
        return self._poll_coverage(
            job_id,
            region,
            positions=response.get('positions'),
            poll_interval=poll_interval,
            max_wait=max_wait,
        )

    def get_sequence_coverage(self, sequence: str, sequence_name: str = '',
                              metadata_filters: Optional[Dict[str, Any]] = None,
                              poll_interval: int = 2, max_wait: int = 300) -> 'CoverageResult':
        """
        Get coverage for an arbitrary DNA sequence.

        Args:
            sequence: DNA sequence
            sequence_name: Optional name for the sequence
            metadata_filters: Optional metadata filters
            poll_interval: How often to check for completion (seconds)
            max_wait: Maximum time to wait for completion (seconds)

        Returns:
            CoverageResult object
        """
        data = {
            'sequence': sequence,
        }
        if sequence_name:
            data['sequence_name'] = sequence_name
        if metadata_filters:
            data['metadata_filters'] = metadata_filters

        response = self._request('POST', '/api/genome-browser/search-sequence', json=data)
        job_id = response.get('job_id')

        if not job_id:
            raise MalvaAPIError("No job_id received from sequence coverage search")

        region = sequence_name or f"seq:{sequence[:20]}..."
        return self._poll_coverage(
            job_id,
            region,
            positions=response.get('positions'),
            poll_interval=poll_interval,
            max_wait=max_wait,
        )

    def _poll_coverage(self, job_id: str, region: str = '',
                       positions: Optional[List[int]] = None,
                       poll_interval: int = 2, max_wait: int = 300) -> 'CoverageResult':
        """
        Internal helper to poll for coverage results.

        Newer explorer-backed APIs expose the same /api/coverage/status and
        /api/coverage/results endpoints used by the web Coverage Explorer. Use
        those when submit returned probe positions, and fall back to the older
        /api/genome-browser/coverage endpoint for compatibility.
        """
        start_time = time.time()

        if positions is not None:
            try:
                while time.time() - start_time < max_wait:
                    response = self._request('GET', f'/api/coverage/status/{job_id}')
                    status = response.get('status')
                    if status == 'completed':
                        result = self._request(
                            'POST',
                            f'/api/coverage/results/{job_id}',
                            json={'positions': positions, 'filters': {}},
                        )
                        coverage_data = result.get('coverage', result.get('coverage_data', result))
                        coverage_data['job_id'] = job_id
                        coverage_data['region'] = region
                        if 'probe_data' in result:
                            coverage_data['probe_data'] = result['probe_data']
                        return CoverageResult(coverage_data, self)
                    if status == 'error':
                        error_msg = response.get('error', 'Unknown error')
                        raise MalvaAPIError(f"Coverage job failed: {error_msg}")
                    logger.info(f"Coverage job {job_id} status: {status}")
                    time.sleep(poll_interval)
                raise MalvaAPIError(f"Coverage job {job_id} timed out after {max_wait} seconds")
            except MalvaAPIError as exc:
                if not self._is_not_found_error(exc):
                    raise
                logger.info("Explorer coverage result endpoint unavailable; falling back to genome-browser coverage endpoint")

        while time.time() - start_time < max_wait:
            response = self._request('GET', f'/api/genome-browser/coverage/{job_id}')
            status = response.get('status')

            if status == 'completed':
                coverage_data = response.get('coverage_data', response)
                coverage_data['job_id'] = job_id
                coverage_data['region'] = region
                if 'probe_data' in response:
                    coverage_data['probe_data'] = response['probe_data']
                return CoverageResult(coverage_data, self)
            elif status == 'error':
                error_msg = response.get('error', 'Unknown error')
                raise MalvaAPIError(f"Coverage job failed: {error_msg}")

            logger.info(f"Coverage job {job_id} status: {status}")
            time.sleep(poll_interval)

        raise MalvaAPIError(f"Coverage job {job_id} timed out after {max_wait} seconds")

    def get_coverage_data(self, job_id: str, **filters) -> 'CoverageResult':
        """
        Retrieve coverage data for an existing job.

        Args:
            job_id: Coverage job ID
            **filters: Optional metadata filters to apply

        Returns:
            CoverageResult object
        """
        params = {}
        if filters:
            for key, value in filters.items():
                if isinstance(value, list):
                    params[key] = value
                else:
                    params[key] = [value] if value else []

        response = self._request('GET', f'/api/genome-browser/coverage/{job_id}', params=params)
        coverage_data = response.get('coverage_data', response)
        coverage_data['job_id'] = job_id
        return CoverageResult(coverage_data, self)

    def download_coverage_wig(self, job_id: str, output_path: str,
                              poll_interval: int = 2, max_wait: int = 120,
                              **filters) -> str:
        """
        Download coverage data as a WIG file.

        The endpoint returns 202 while the file is being generated,
        so we poll until the file is ready.

        Args:
            job_id: Coverage job ID
            output_path: Path to save the WIG file
            poll_interval: How often to check for completion (seconds)
            max_wait: Maximum time to wait (seconds)
            **filters: Optional metadata filters

        Returns:
            Path to the saved file
        """
        url = urljoin(self.base_url, f'/api/genome-browser/coverage-file/{job_id}')

        params = {}
        if filters:
            for key, value in filters.items():
                if isinstance(value, list):
                    params[key] = value
                else:
                    params[key] = [value] if value else []

        start_time = time.time()
        while time.time() - start_time < max_wait:
            response = self.session.get(url, params=params, timeout=self.timeout)

            if response.status_code == 200:
                output = Path(output_path)
                output.parent.mkdir(parents=True, exist_ok=True)
                output.write_text(response.text)
                logger.info(f"Coverage WIG saved to {output_path}")
                return str(output)
            elif response.status_code == 202:
                logger.info(f"WIG file still generating for job {job_id}...")
                time.sleep(poll_interval)
            else:
                raise MalvaAPIError(f"Failed to download WIG file: HTTP {response.status_code}")

        raise MalvaAPIError(f"WIG download timed out after {max_wait} seconds")

    def get_coverage_filter_options(self, job_id: str) -> Dict[str, Any]:
        """
        Get available filter options for a coverage result.

        Args:
            job_id: Coverage job ID

        Returns:
            Dictionary with available filter options
        """
        return self._request('GET', f'/api/genome-browser/filter-options/{job_id}')

    # ── Dataset Discovery Methods ──

    def get_datasets_hierarchy(self) -> Dict[str, Any]:
        """
        Get the full dataset → study → sample hierarchy.

        Returns:
            Dictionary with the full hierarchy tree
        """
        return self._request('GET', '/samples/api/datasets/hierarchy')

    def get_dataset_studies(self, dataset_id: str,
                            page: int = 1, page_size: int = 50) -> Dict[str, Any]:
        """
        Get studies belonging to a dataset.

        Args:
            dataset_id: Dataset identifier
            page: Page number (1-based)
            page_size: Number of studies per page

        Returns:
            Dictionary with studies and pagination info
        """
        params = {'page': page, 'page_size': page_size}
        return self._request('GET', f'/samples/api/datasets/{dataset_id}/studies', params=params)

    def get_study_samples(self, dataset_id: str, study_name: str,
                          page: int = 1, page_size: int = 50) -> Dict[str, Any]:
        """
        Get samples belonging to a study within a dataset.

        Args:
            dataset_id: Dataset identifier
            study_name: Study name
            page: Page number (1-based)
            page_size: Number of samples per page

        Returns:
            Dictionary with samples and pagination info
        """
        params = {'page': page, 'page_size': page_size}
        return self._request('GET', f'/samples/api/datasets/{dataset_id}/studies/{study_name}/samples',
                             params=params)

    def get_sample_details(self, sample_uuid: str) -> Dict[str, Any]:
        """
        Get full metadata for a single sample.

        Args:
            sample_uuid: UUID of the sample

        Returns:
            Dictionary with complete sample metadata
        """
        return self._request('GET', f'/samples/api/sample/{sample_uuid}')

    def get_filter_values(self, column: str, search: str = '',
                          limit: int = 50) -> List[str]:
        """
        Get unique values for a filter column, with optional search.

        Args:
            column: Filter column name (e.g., 'organ', 'species')
            search: Optional text to filter values by
            limit: Maximum number of values to return

        Returns:
            List of unique values for the column
        """
        params = {'limit': limit}
        if search:
            params['search'] = search
        response = self._request('GET', f'/samples/api/filters/{column}/values', params=params)
        return response.get('values', [])

    def get_overview_stats(self) -> Dict[str, Any]:
        """
        Get platform-wide statistics.

        Returns:
            Dictionary with overview statistics (total samples, studies, etc.)
        """
        return self._request('GET', '/samples/api/overview')

    # ── Coexpression Analysis Methods ──

    def get_umap_coordinates(self, dataset_id: str) -> 'UMAPCoordinates':
        """
        Get UMAP coordinates for a dataset's metacells.

        Args:
            dataset_id: Dataset identifier

        Returns:
            UMAPCoordinates object with x, y, cluster, and sample data
        """
        response = self._request('GET', f'/api/coexpression/umap/{dataset_id}')
        return UMAPCoordinates(response, self)

    def _coexpression_from_response(self, response: Dict[str, Any],
                                    task_endpoint: str,
                                    poll_interval: int = 2,
                                    max_wait: int = 300) -> 'CoexpressionResult':
        """Poll async coexpression task responses and return a result object."""
        status = response.get('status')
        task_id = response.get('task_id')
        if status not in {'pending', 'running'} or not task_id:
            return CoexpressionResult(response, self)

        start_time = time.time()
        while time.time() - start_time < max_wait:
            task = self._request('GET', f'{task_endpoint.rstrip("/")}/{task_id}')
            task_status = task.get('status')
            if task_status == 'completed':
                return CoexpressionResult(task, self)
            if task_status in {'failed', 'error'}:
                raise MalvaAPIError(f"Coexpression task failed: {task.get('error', 'unknown')}")
            logger.info(f"Coexpression task {task_id} status: {task_status}")
            time.sleep(poll_interval)
        raise MalvaAPIError(f"Coexpression task {task_id} timed out after {max_wait} seconds")

    @staticmethod
    def _coexpression_task_endpoint(endpoint: str) -> str:
        if endpoint.startswith('/api/coexpression/'):
            return '/api/coexpression/tasks'
        return '/tasks'

    def project_cells(self, cell_ids: List[int], sample_ids: List[int],
                      dataset_id: str,
                      poll_interval: int = 2,
                      max_wait: int = 300,
                      **kwargs) -> 'CoexpressionResult':
        """
        Project explicit cells onto a coexpression index.

        Args:
            cell_ids: Cell IDs to project.
            sample_ids: Encoded sample IDs corresponding to cell_ids.
            dataset_id: Coexpression index or dataset identifier.
            poll_interval: Seconds between task polls when the service runs
                asynchronously.
            max_wait: Maximum seconds to wait for completion.
            **kwargs: Coexpression parameters such as top_n_genes,
                include_go_enrichment, include_umap_scores, and
                min_abs_correlation.
        """
        if len(cell_ids) != len(sample_ids):
            raise ValueError("cell_ids and sample_ids must have the same length")
        if not cell_ids:
            raise ValueError("At least one cell_id is required")

        data = {
            'dataset_id': dataset_id,
            'cell_ids': [int(c) for c in cell_ids],
            'sample_ids': [int(s) for s in sample_ids],
        }
        data.update(kwargs)

        response, endpoint = self._request_first_available(
            'POST',
            ['/api/coexpression/query', '/query'],
            json=data,
        )
        return self._coexpression_from_response(
            response,
            task_endpoint=self._coexpression_task_endpoint(endpoint),
            poll_interval=poll_interval,
            max_wait=max_wait,
        )

    def get_coexpression(self, job_id: str, dataset_id: str,
                         filter_sample_ids: Optional[List[int]] = None,
                         poll_interval: int = 2,
                         max_wait: int = 300,
                         **kwargs) -> 'CoexpressionResult':
        """
        Run a full coexpression analysis for a search job on a dataset.

        Calls ``POST /api/coexpression/query-by-job`` with all enrichment
        options enabled (GO enrichment, UMAP scores, cell-type enrichment,
        tissue breakdown).

        Args:
            job_id: Job ID from a previous search
            dataset_id: Dataset identifier to run coexpression against
            filter_sample_ids: Optional encoded sample IDs. When provided,
                only cells from those samples are projected.
            poll_interval: Seconds between task polls when the service runs
                asynchronously.
            max_wait: Maximum seconds to wait for completion.
            **kwargs: Additional coexpression parameters.

        Returns:
            CoexpressionResult with correlated genes, UMAP scores, GO
            enrichment, cell-type enrichment, and tissue breakdown
        """
        data = {
            'job_id': job_id,
            'dataset_id': dataset_id,
        }
        if filter_sample_ids:
            data['filter_sample_ids'] = [int(s) for s in filter_sample_ids]
        data.update(kwargs)
        response, endpoint = self._request_first_available(
            'POST',
            ['/api/coexpression/query-by-job', '/query-by-job'],
            json=data,
        )
        return self._coexpression_from_response(
            response,
            task_endpoint=self._coexpression_task_endpoint(endpoint),
            poll_interval=poll_interval,
            max_wait=max_wait,
        )

    def get_coexpression_genes(self, job_id: str, dataset_id: str,
                               filter_sample_ids: Optional[List[int]] = None,
                               poll_interval: int = 2,
                               max_wait: int = 300,
                               **kwargs) -> 'CoexpressionResult':
        """
        Lightweight coexpression query returning only correlated genes.

        Same endpoint as :meth:`get_coexpression` but disables GO
        enrichment, UMAP scores, cell-type enrichment, and tissue
        breakdown for faster response.

        Args:
            job_id: Job ID from a previous search
            dataset_id: Dataset identifier
            filter_sample_ids: Optional encoded sample IDs. When provided,
                only cells from those samples are projected.
            poll_interval: Seconds between task polls when the service runs
                asynchronously.
            max_wait: Maximum seconds to wait for completion.
            **kwargs: Additional coexpression parameters.

        Returns:
            CoexpressionResult with only the correlated_genes field populated
        """
        data = {
            'job_id': job_id,
            'dataset_id': dataset_id,
            'include_go_enrichment': False,
            'include_umap_scores': False,
            'include_cell_type_enrichment': False,
            'include_tissue_breakdown': False,
        }
        if filter_sample_ids:
            data['filter_sample_ids'] = [int(s) for s in filter_sample_ids]
        data.update(kwargs)
        response, endpoint = self._request_first_available(
            'POST',
            ['/api/coexpression/query-by-job', '/query-by-job'],
            json=data,
        )
        return self._coexpression_from_response(
            response,
            task_endpoint=self._coexpression_task_endpoint(endpoint),
            poll_interval=poll_interval,
            max_wait=max_wait,
        )


# Convenience functions for common operations
def search_gene(gene: str, base_url: str = "https://malva.mdc-berlin.de") -> SearchResult:
    """
    Quick gene search

    Args:
        gene: Gene symbol to search for
        base_url: Malva API base URL

    Returns:
        SearchResult object
    """
    client = MalvaClient(base_url)
    return client.search(gene)


def search_sequence(sequence: str, base_url: str = "https://malva.mdc-berlin.de") -> SearchResult:
    """
    Quick sequence search (convenience wrapper around :meth:`MalvaClient.search_sequences`).

    Args:
        sequence: DNA sequence to search for
        base_url: Malva API base URL

    Returns:
        SearchResult object
    """
    client = MalvaClient(base_url)
    return client.search_sequences(sequence)
