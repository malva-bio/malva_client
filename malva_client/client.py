import json
import time
import requests
import numpy as np
import pandas as pd
from typing import Dict, List, Any, Optional, Union, Tuple
from pathlib import Path
import logging
from urllib.parse import urljoin
import h5py
import anndata as ad

from malva_client.exceptions import MalvaAPIError, AuthenticationError, SearchError, QuotaExceededError
from malva_client.models import SearchResult, SingleCellResult, CoverageResult, CoexpressionResult, UMAPCoordinates
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

    def __init__(self, base_url: str = None, api_token: str = None, timeout: int = 300, verify_ssl: bool = True):
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
                except json.JSONDecodeError:
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
                except json.JSONDecodeError:
                    raise QuotaExceededError("Search quota exceeded")
            
            elif response.status_code in [401, 403]:
                try:
                    error_data = response.json()
                    error_msg = error_data.get('error', error_data.get('detail', 'Authentication failed'))
                except json.JSONDecodeError:
                    error_msg = "Authentication failed"
                raise AuthenticationError(error_msg)
            
            elif not response.ok:
                try:
                    error_data = response.json()
                    error_msg = error_data.get('detail', error_data.get('error', f'HTTP {response.status_code}'))
                    raise MalvaAPIError(f"API error: {error_msg}")
                except json.JSONDecodeError:
                    raise MalvaAPIError(f"HTTP {response.status_code}: {response.text[:200]}")
            
            try:
                result = response.json()
                # Return a lazy SearchResults wrapper for compact_v2 to avoid
                # O(N_cells × N_genes) list expansion at parse time.
                if isinstance(result, dict) and result.get('results', {}).get('_format') == 'compact_v2':
                    return SearchResults(result)
                return result
            except json.JSONDecodeError as e:
                content_type = response.headers.get('content-type', 'unknown')
                raise MalvaAPIError(f"Expected JSON response but got {content_type}. Content: {response.text[:200]}")
            
        except requests.exceptions.Timeout:
            raise MalvaAPIError(f"Request timed out after {self.timeout} seconds")
        except requests.exceptions.ConnectionError:
            raise MalvaAPIError(f"Could not connect to {self.base_url}. Check if the server is running.")
        except requests.exceptions.RequestException as e:
            raise MalvaAPIError(f"Request failed: {e}")

    def search(self, query: str, wait_for_completion: bool = True,
               poll_interval: int = 2, max_wait: int = 300, aggregate_expression: bool = True,
               window_size: Optional[int] = None,
               threshold: Optional[float] = None,
               stranded: Optional[bool] = None) -> SearchResult:
        """
        Perform a natural language or gene search

        Args:
            query: Natural language query or gene symbol (e.g., "BRCA1 expression in cancer")
            wait_for_completion: Whether to wait for results or return immediately
            poll_interval: How often to check for completion (seconds)
            max_wait: Maximum time to wait for completion (seconds)
            aggregate_expression: Whether to aggregate expression by cell type
            window_size: Sliding window size in k-mers. Controls how many consecutive
                k-mers are evaluated per window. Larger values suit transcript-level
                queries; smaller values (e.g., 24) suit junction or SNV detection.
            threshold: Match threshold between 0.0 and 1.0. Fraction of k-mers in
                a window that must match. Lower values allow more mismatches.
            stranded: If True, restrict search to the forward strand only.

        Returns:
            SearchResult object containing the results
        """
        if window_size is not None and window_size < 1:
            raise ValueError("window_size must be a positive integer")
        if threshold is not None and not (0.0 <= threshold <= 1.0):
            raise ValueError("threshold must be between 0.0 and 1.0")

        data = {
            'query': query,
            'wait_for_completion': wait_for_completion,
            'max_wait': max_wait,
            'poll_interval': poll_interval,
            'aggregate_expression': aggregate_expression
        }
        if window_size is not None:
            data['window_size'] = window_size
        if threshold is not None:
            data['threshold'] = threshold
        if stranded is not None:
            data['stranded'] = stranded
        
        response = self._request('POST', '/search', json=data)
        
        if response.get('status') == 'completed':
            job_id = response.get('job_id')
            logger.info(f"Search completed with job ID: {job_id}")
            return SearchResult({'job_id': job_id, 'status': 'completed'}, self)
        
        job_id = response['job_id']
        logger.info(f"Search submitted with job ID: {job_id}")
        
        if not wait_for_completion or response.get('status') == 'pending':
            return SearchResult({'job_id': job_id, 'status': 'pending'}, self)
        
        if not response.get('success', True):
            error_msg = response.get('error', 'Unknown error')
            raise MalvaAPIError(f"Search failed: {error_msg}")
        
        start_time = time.time()
        while time.time() - start_time < max_wait:
            try:
                results_response = self._request('GET', f'/search/{job_id}')
                status = results_response.get('status')
                
                logger.info(f"Search status: {status}")
                
                if status == 'completed':
                    logger.info(f"Search completed with job ID: {job_id}")
                    # Use the results_response (has job_id + status); expression data
                    # is fetched lazily via /api/expression-data/<job_id>/results on
                    # first access to SearchResult.df
                    return SearchResult(
                        {'job_id': job_id, 'status': 'completed'}, self
                    )
                elif status == 'error':
                    error_msg = results_response.get('error', 'Unknown error')
                    raise MalvaAPIError(f"Search failed: {error_msg}")
                
            except MalvaAPIError as e:
                if "404" in str(e) or "not found" in str(e).lower():
                    pass
                else:
                    raise
            
            time.sleep(poll_interval)
        
        raise MalvaAPIError(f"Search timed out after {max_wait} seconds")
    
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

    def search_cells(self, query: str, wait_for_completion: bool = True,
                    poll_interval: int = 2, max_wait: int = 300,
                    window_size: Optional[int] = None,
                    threshold: Optional[float] = None,
                    stranded: Optional[bool] = None) -> SingleCellResult:
        """
        Perform a natural language or gene search and returns individual cells

        Args:
            query: Natural language query, gene symbol, or DNA sequence
            wait_for_completion: Whether to wait for results or return immediately
            poll_interval: How often to check for completion (seconds)
            max_wait: Maximum time to wait for completion (seconds)
            window_size: Sliding window size in k-mers (see :doc:`query_parameters`)
            threshold: Match threshold between 0.0 and 1.0
            stranded: If True, restrict search to the forward strand only

        Returns:
            SingleCellResult object with cell-level resolution
        """
        if window_size is not None and window_size < 1:
            raise ValueError("window_size must be a positive integer")
        if threshold is not None and not (0.0 <= threshold <= 1.0):
            raise ValueError("threshold must be between 0.0 and 1.0")

        data = {
            'query': query,
            'wait_for_completion': wait_for_completion,
            'max_wait': max_wait,
            'poll_interval': poll_interval,
            'aggregate_expression': False
        }
        if window_size is not None:
            data['window_size'] = window_size
        if threshold is not None:
            data['threshold'] = threshold
        if stranded is not None:
            data['stranded'] = stranded
        
        response = self._request('POST', '/search', json=data)
        
        if response.get('status') == 'completed' and 'results' in response:
            logger.info(f"Search completed immediately with job ID: {response.get('job_id', 'unknown')}")
            return SingleCellResult(response, self)
        
        # If not completed immediately, get the job_id for polling
        job_id = response.get('job_id')
        if not job_id:
            raise MalvaAPIError("No job_id received from server")
        
        logger.info(f"Search submitted with job ID: {job_id}")
        
        if not wait_for_completion:
            if response.get('status') == 'completed' and 'results' in response:
                return SingleCellResult(response, self)
            return SingleCellResult({'job_id': job_id, 'status': 'pending'}, self)
        
        if response.get('status') == 'completed' and 'results' in response:
            return SingleCellResult(response, self)

        start_time = time.time()
        while time.time() - start_time < max_wait:
            try:
                results_response = self._request('GET', f'/search/{job_id}')
                
                if not results_response:
                    if response.get('status') == 'completed' and 'results' in response:
                        logger.info(f"Using cached results for job {job_id}")
                        return SingleCellResult(response, self)
                    else:
                        logger.warning(f"No results found for job {job_id}")
                        time.sleep(poll_interval)
                        continue
                
                status = results_response.get('status')
                logger.info(f"Search status: {status}")
                
                if status == 'completed' and 'results' in results_response:
                    logger.info(f"Search completed with job ID: {job_id}")
                    return SingleCellResult(results_response, self)
                elif status == 'error':
                    error_msg = results_response.get('error', 'Unknown error')
                    raise MalvaAPIError(f"Search failed: {error_msg}")
                    
            except MalvaAPIError as e:
                if "404" in str(e) or "not found" in str(e).lower():
                    if response.get('status') == 'completed' and 'results' in response:
                        logger.info(f"Job {job_id} not found on server, but we have cached results")
                        return SingleCellResult(response, self)
                    pass
                else:
                    raise
            
            time.sleep(poll_interval)
        
        if response.get('status') == 'completed' and 'results' in response:
            logger.warning(f"Polling timed out but returning initial results for job {job_id}")
            return SingleCellResult(response, self)
        
        raise MalvaAPIError(f"Search timed out after {max_wait} seconds")
    
    def search_sequence(self, sequence: str, **kwargs) -> SearchResult:
        """
        Search for a specific DNA sequence

        Args:
            sequence: DNA sequence to search for
            **kwargs: Additional arguments passed to search(), including
                ``window_size``, ``threshold``, and ``stranded``

        Returns:
            SearchResult object
        """
        if len(sequence) > 500000:
            raise ValueError("Sequence length cannot exceed 500kb (500,000 nucleotides)")
        
        valid_nucleotides = set('ATGCUN')
        if not all(c.upper() in valid_nucleotides for c in sequence):
            raise ValueError("Sequence contains invalid nucleotides. Only A, T, G, C, U, N are allowed.")
        
        return self.search(sequence, **kwargs)

    def search_sequences(self, sequences: List[str],
                         wait_for_completion: bool = True,
                         poll_interval: int = 2, max_wait: int = 300,
                         aggregate_expression: bool = True) -> SearchResult:
        """
        Search for multiple DNA sequences in a single batch request

        Sends all sequences as FASTA content in one API call using the
        ``fasta_sequences`` field. Only one search quota slot is consumed.

        Args:
            sequences: List of DNA sequences to search for
            wait_for_completion: Whether to wait for results or return immediately
            poll_interval: How often to check for completion (seconds)
            max_wait: Maximum time to wait for completion (seconds)
            aggregate_expression: Whether to aggregate expression by cell type

        Returns:
            SearchResult object containing results for all sequences
        """
        if not sequences:
            raise ValueError("At least one sequence is required")

        valid_nucleotides = set('ATGCUN')
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
                f"Total nucleotides ({total_nt:,}) exceed the 100,000 limit"
            )

        # Build FASTA content
        fasta_lines = []
        for i, seq in enumerate(sequences, 1):
            fasta_lines.append(f">seq_{i}")
            fasta_lines.append(seq)
        fasta_content = "\n".join(fasta_lines) + "\n"

        data = {
            'query': 'batch sequence search',
            'fasta_sequences': fasta_content,
            'wait_for_completion': wait_for_completion,
            'max_wait': max_wait,
            'poll_interval': poll_interval,
            'aggregate_expression': aggregate_expression,
        }

        logger.info(
            f"Submitting batch search: {len(sequences)} sequences, "
            f"{total_nt:,} total nucleotides"
        )

        response = self._request('POST', '/search', json=data)

        if response.get('status') == 'completed':
            logger.info(f"Batch search completed with job ID: {response['job_id']}")
            return SearchResult(response['results'], self)

        job_id = response['job_id']
        logger.info(f"Batch search submitted with job ID: {job_id}")

        if not wait_for_completion or response.get('status') == 'pending':
            return SearchResult({'job_id': job_id, 'status': 'pending'}, self)

        if not response.get('success', True):
            error_msg = response.get('error', 'Unknown error')
            raise MalvaAPIError(f"Batch search failed: {error_msg}")

        start_time = time.time()
        while time.time() - start_time < max_wait:
            try:
                results_response = self._request('GET', f'/search/{job_id}')
                status = results_response.get('status')

                logger.info(f"Batch search status: {status}")

                if status == 'completed':
                    logger.info(f"Batch search completed with job ID: {job_id}")
                    return SearchResult(results_response.get('results', results_response), self)
                elif status == 'error':
                    error_msg = results_response.get('error', 'Unknown error')
                    raise MalvaAPIError(f"Batch search failed: {error_msg}")

            except MalvaAPIError as e:
                if "404" in str(e) or "not found" in str(e).lower():
                    pass
                else:
                    raise

            time.sleep(poll_interval)

        raise MalvaAPIError(f"Batch search timed out after {max_wait} seconds")

    def search_genes(self, genes: List[str],
                     wait_for_completion: bool = True,
                     poll_interval: int = 2, max_wait: int = 300,
                     aggregate_expression: bool = True) -> SearchResult:
        """
        Search for multiple genes in a single request

        Sends a comma-separated gene list as the query string.

        Args:
            genes: List of gene symbols (e.g., ``["BRCA1", "TP53"]``)
            wait_for_completion: Whether to wait for results or return immediately
            poll_interval: How often to check for completion (seconds)
            max_wait: Maximum time to wait for completion (seconds)
            aggregate_expression: Whether to aggregate expression by cell type

        Returns:
            SearchResult object
        """
        if not genes:
            raise ValueError("At least one gene is required")
        if len(genes) > 10:
            raise ValueError(
                f"Too many genes ({len(genes)}). Maximum is 10 per request."
            )

        gene_query = ", ".join(genes)

        data = {
            'query': gene_query,
            'wait_for_completion': wait_for_completion,
            'max_wait': max_wait,
            'poll_interval': poll_interval,
            'aggregate_expression': aggregate_expression,
        }

        logger.info(f"Submitting gene search: {gene_query}")

        response = self._request('POST', '/search', json=data)

        if response.get('status') == 'completed':
            logger.info(f"Gene search completed with job ID: {response['job_id']}")
            return SearchResult(response['results'], self)

        job_id = response['job_id']
        logger.info(f"Gene search submitted with job ID: {job_id}")

        if not wait_for_completion or response.get('status') == 'pending':
            return SearchResult({'job_id': job_id, 'status': 'pending'}, self)

        if not response.get('success', True):
            error_msg = response.get('error', 'Unknown error')
            raise MalvaAPIError(f"Gene search failed: {error_msg}")

        start_time = time.time()
        while time.time() - start_time < max_wait:
            try:
                results_response = self._request('GET', f'/search/{job_id}')
                status = results_response.get('status')

                logger.info(f"Gene search status: {status}")

                if status == 'completed':
                    logger.info(f"Gene search completed with job ID: {job_id}")
                    return SearchResult(results_response.get('results', results_response), self)
                elif status == 'error':
                    error_msg = results_response.get('error', 'Unknown error')
                    raise MalvaAPIError(f"Gene search failed: {error_msg}")

            except MalvaAPIError as e:
                if "404" in str(e) or "not found" in str(e).lower():
                    pass
                else:
                    raise

            time.sleep(poll_interval)

        raise MalvaAPIError(f"Gene search timed out after {max_wait} seconds")
    
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
            # Request the download URL from frontend
            response = self._request('GET', f'/api/samples/{sample_uuid}/download')
            
            if 'download_url' not in response:
                raise MalvaAPIError("Download URL not provided by server")
            
            download_url = response['download_url']
            filename = response.get('filename', f'{sample_uuid}.h5ad')
            file_size = response.get('file_size', 0)
            expires_in = response.get('expires_in_seconds', 0)
            
            if expires_in <= 0:
                raise MalvaAPIError("Download token has expired")
            
            logger.info(f"Downloading {filename} ({file_size:,} bytes)")
            
            # Download the file with progress tracking
            file_response = self.session.get(download_url, stream=True)
            file_response.raise_for_status()
            
            # Handle download errors
            if file_response.status_code == 404:
                raise MalvaAPIError("Download token expired or file not found")
            elif file_response.status_code == 403:
                raise MalvaAPIError("Access denied to file")
            
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
            response = self._request('GET', f'/search/{job_id}')

            # Update local storage
            save_job_status(
                job_id=job_id,
                status=response.get('status', 'unknown'),
                server_url=self.base_url,
                results_data=response.get('results') if response.get('status') == 'completed' else None,
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

        # First check if we have results cached locally
        local_status = get_local_status(job_id)
        if local_status and local_status.get('results_data'):
            logger.info(f"Using cached results for job {job_id}")
            return SearchResult(local_status['results_data'], self)

        # Otherwise fetch from server
        status_response = self.get_job_status(job_id)

        if status_response['status'] == 'completed':
            results_data = status_response.get('results', {})

            # Cache results locally
            save_job_status(
                job_id=job_id,
                status='completed',
                server_url=self.base_url,
                results_data=results_data
            )

            return SearchResult(results_data, self)

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
            status_response = self.get_job_status(job_id)
            status = status_response.get('status')

            if status == 'completed':
                return self.get_job_results(job_id)
            elif status == 'error':
                error_msg = status_response.get('error', 'Unknown error')
                raise SearchError(f"Job {job_id} failed: {error_msg}")
            elif status in ['pending', 'running']:
                logger.info(f"Job {job_id} status: {status}")
                time.sleep(poll_interval)
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
        return self._poll_coverage(job_id, region, poll_interval=poll_interval, max_wait=max_wait)

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
        return self._poll_coverage(job_id, region, poll_interval=poll_interval, max_wait=max_wait)

    def _poll_coverage(self, job_id: str, region: str = '',
                       poll_interval: int = 2, max_wait: int = 300) -> 'CoverageResult':
        """
        Internal helper to poll for coverage results.

        Args:
            job_id: Coverage job ID
            region: Region description for the result object
            poll_interval: How often to check for completion (seconds)
            max_wait: Maximum time to wait for completion (seconds)

        Returns:
            CoverageResult object
        """
        start_time = time.time()

        while time.time() - start_time < max_wait:
            response = self._request('GET', f'/api/genome-browser/coverage/{job_id}')
            status = response.get('status')

            if status == 'completed':
                coverage_data = response.get('coverage_data', response)
                coverage_data['job_id'] = job_id
                coverage_data['region'] = region
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

    def get_coexpression(self, job_id: str, dataset_id: str) -> 'CoexpressionResult':
        """
        Run a full coexpression analysis for a search job on a dataset.

        Calls ``POST /api/coexpression/query-by-job`` with all enrichment
        options enabled (GO enrichment, UMAP scores, cell-type enrichment,
        tissue breakdown).

        Args:
            job_id: Job ID from a previous search
            dataset_id: Dataset identifier to run coexpression against

        Returns:
            CoexpressionResult with correlated genes, UMAP scores, GO
            enrichment, cell-type enrichment, and tissue breakdown
        """
        data = {
            'job_id': job_id,
            'dataset_id': dataset_id,
        }
        response = self._request('POST', '/api/coexpression/query-by-job', json=data)
        return CoexpressionResult(response, self)

    def get_coexpression_genes(self, job_id: str, dataset_id: str) -> 'CoexpressionResult':
        """
        Lightweight coexpression query returning only correlated genes.

        Same endpoint as :meth:`get_coexpression` but disables GO
        enrichment, UMAP scores, cell-type enrichment, and tissue
        breakdown for faster response.

        Args:
            job_id: Job ID from a previous search
            dataset_id: Dataset identifier

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
        response = self._request('POST', '/api/coexpression/query-by-job', json=data)
        return CoexpressionResult(response, self)


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
    Quick sequence search

    Args:
        sequence: DNA sequence to search for
        base_url: Malva API base URL

    Returns:
        SearchResult object
    """
    client = MalvaClient(base_url)
    return client.search_sequence(sequence)