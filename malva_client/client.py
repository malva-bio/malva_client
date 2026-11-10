import json
import logging
import os
import tarfile
import tempfile
import time
import warnings
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple, Union
from urllib.parse import urljoin, urlparse, urlunparse

import anndata as ad
import msgpack
import numpy as np
import pandas as pd
import requests

from malva_client import __version__ as _CLIENT_VERSION
from malva_client.config import Config
from malva_client.exceptions import (
    AuthenticationError,
    ClientVersionError,
    MalvaAPIError,
    QuotaExceededError,
    SearchError,
)
from malva_client.models import (
    CellExpressionMatrixResult,
    CellExpressionStream,
    CoexpressionResult,
    CoverageResult,
    DifferentialExpressionResult,
    MetacellMatchResult,
    SearchResult,
    UMAPCoordinates,
)

# Stable empty-result schema for federated searches. Kept next to the client so
# an empty (but successful) federation still returns the same columns as a
# populated result instead of a bare provider-only frame.
_FEDERATED_RESULT_COLUMNS = [
    "provider",
    "provider_id",
    "sample_idx",
    "cell_type_idx",
    "sample_id",
    "cell_type",
    "gene_sequence",
    "rel",
    "exp",
    "pct",
    "raw_kmers",
    "cell_count",
    "sample_celltype_id",
]

# One logical operation can have a canonical and a legacy HTTP path. Keeping
# both variants side by side lets every call site use one readable operation
# name instead of branching on the endpoint style itself.
_ENDPOINT_PATHS = {
    "health": ("/api/health", "/health"),
    "quota": ("/api/account/quota", "/api/quota-status"),
    "datasets": ("/api/datasets", "/datasets"),
    "samples": ("/api/samples", "/samples/api/search"),
    "sample_details": ("/api/samples/{sample_uuid}", "/samples/api/sample/{sample_uuid}"),
    "sample_metadata": ("/api/samples/metadata", "/api/sample-metadata"),
    "studies": ("/api/metadata/studies", "/samples/api/studies"),
    "overview": ("/api/metadata/overview", "/samples/api/overview"),
    "available_filters": ("/api/metadata/available-filters", "/samples/api/filters"),
    "filter_values": (
        "/api/metadata/filters/{column}/values",
        "/samples/api/filters/{column}/values",
    ),
    "datasets_hierarchy": (
        "/api/metadata/datasets/hierarchy",
        "/samples/api/datasets/hierarchy",
    ),
    "dataset_studies": (
        "/api/metadata/datasets/{dataset_id}/studies",
        "/samples/api/datasets/{dataset_id}/studies",
    ),
    "study_samples": (
        "/api/metadata/datasets/{dataset_id}/studies/{study_name}/samples",
        "/samples/api/datasets/{dataset_id}/studies/{study_name}/samples",
    ),
    "search": ("/api/search", "/search"),
    "search_job": ("/api/search/{job_id}", "/search/{job_id}"),
    "search_job_status": ("/api/search/{job_id}/status", "/search/{job_id}/status"),
    "search_job_results": (
        "/api/search/{job_id}/results",
        "/api/expression-data/{job_id}/results",
    ),
    "search_job_cells": (
        "/api/search/{job_id}/cells",
        "/api/expression/results/{job_id}/cells",
    ),
    "search_job_cell_ids": (
        "/api/search/{job_id}/cell-ids",
        "/search/{job_id}/cell-ids",
    ),
    "search_job_global_cell_ids": (
        "/api/search/{job_id}/global-cell-ids",
        "/search/{job_id}/global-cell-ids",
    ),
    "search_federated": ("/api/search", "/api/expression/submit"),
    "search_federated_job": (
        "/api/search/federated/{job_id}",
        "/api/expression/federated/{job_id}",
    ),
    "search_shard_results": (
        "/api/search/{job_id}/results",
        "/api/expression/results/{job_id}",
    ),
    "cells_by_metadata": ("/api/cells/by-metadata", "/api/cells/by-filter"),
    "cell_expression_export": (
        "/api/exports/cell-expression",
        "/exports/cell-expression",
    ),
    "cell_expression_export_job": (
        "/api/exports/cell-expression/{export_id}",
        "/exports/cell-expression/{export_id}",
    ),
    "cell_expression_export_download": (
        "/api/exports/cell-expression/{export_id}/download",
        "/exports/cell-expression/{export_id}/download",
    ),
    "all_cells_export": ("/api/exports/all-cells", "/exports/all-cells"),
    "all_cells_export_job": (
        "/api/exports/all-cells/{export_id}",
        "/exports/all-cells/{export_id}",
    ),
    "all_cells_export_download": (
        "/api/exports/all-cells/{export_id}/download",
        "/exports/all-cells/{export_id}/download",
    ),
    "download_token": ("/api/downloads/{token}", "/api/download/{token}"),
}

_ENDPOINT_STYLES = {"auto", "canonical", "legacy"}


def _expr_data_from_columnar(expr_col: dict) -> dict:
    """Reconstruct row-major expression_data from compact_v2 columnar format."""
    si = expr_col.get("si", [])
    n = len(si)
    if n == 0:
        return {
            "data": [],
            "samples": expr_col.get("samples", []),
            "cell_types": expr_col.get("cell_types", []),
            "columns": expr_col.get("columns", []),
            "celltype_sample_counts": expr_col.get("celltype_sample_counts", {}),
        }
    ci = expr_col.get("ci", [])
    v0 = expr_col.get("v0", [])
    v1 = expr_col.get("v1", [])
    v2 = expr_col.get("v2", [])
    v3 = expr_col.get("v3", [])
    v4 = expr_col.get("v4", [])
    return {
        "data": [[si[j], ci[j], v0[j], v1[j], v2[j], v3[j], v4[j]] for j in range(n)],
        "samples": expr_col.get("samples", []),
        "cell_types": expr_col.get("cell_types", []),
        "columns": expr_col.get("columns", []),
        "celltype_sample_counts": expr_col.get("celltype_sample_counts", {}),
    }


def _looks_like_raw_counts(X, sample_size: int = 10000) -> bool:
    """Heuristic: non-negative, finite, integer-valued expression values.

    Works for dense arrays and scipy sparse matrices without densifying.
    Used only to warn the user; the data is never modified.
    """
    from scipy import sparse

    if sparse.issparse(X):
        data = X.data
        if data.size == 0:
            return False
        sample = data[:sample_size]
    else:
        arr = np.asarray(X)
        if arr.size == 0:
            return False
        sample = arr.ravel()[:sample_size]
    sample = np.asarray(sample)
    if sample.size == 0 or not np.all(np.isfinite(sample)):
        return False
    if float(sample.min()) < 0 or float(sample.max()) <= 0:
        return False
    return bool(np.allclose(sample, np.round(sample)))


def _cluster_pseudobulk(X, codes, n_groups: int):
    """Mean expression per cluster using sparse group-indicator multiplication.

    Returns a CSR matrix (n_groups × n_vars).  Only non-zero entries are kept.
    """
    from scipy import sparse

    n_obs = X.shape[0]
    keep = codes >= 0
    kept_codes = codes[keep]
    kept_rows = np.flatnonzero(keep)
    indicator = sparse.csr_matrix(
        (np.ones(len(kept_rows), dtype=np.float64), (kept_codes, kept_rows)),
        shape=(n_groups, n_obs),
    )
    Xcsr = X.tocsr() if sparse.issparse(X) else sparse.csr_matrix(np.asarray(X))
    sums = indicator @ Xcsr
    counts = np.bincount(kept_codes, minlength=n_groups).astype(np.float64)
    inv_counts = np.zeros(n_groups, dtype=np.float64)
    nonzero = counts > 0
    inv_counts[nonzero] = 1.0 / counts[nonzero]
    means = sparse.diags(inv_counts) @ sums
    return means.tocsr()


class GeneResult:
    """Lazy per-gene view over compact_v2 data; per-cell arrays expanded on access."""

    def __init__(self, gdata: dict, global_cells, global_samples):
        self._d = gdata
        self._gc = global_cells
        self._gs = global_samples

    @property
    def ncells(self) -> int:
        return len(self._d.get("_ci") or self._d.get("cell") or [])

    @property
    def expression_data(self) -> dict:
        ed = self._d.get("expression_data", {})
        if isinstance(ed, dict) and ed.get("_fmt") == "col":
            return _expr_data_from_columnar(ed)
        return ed

    @property
    def cells(self) -> np.ndarray:
        ci = self._d.get("_ci")
        if ci is not None and self._gc is not None:
            return np.asarray(self._gc, dtype=np.uint32)[np.asarray(ci, dtype=np.int32)]
        return np.asarray(self._d.get("cell", []), dtype=np.uint64)

    @property
    def samples(self) -> np.ndarray:
        ci = self._d.get("_ci")
        if ci is not None and self._gs is not None:
            return np.asarray(self._gs, dtype=np.uint64)[np.asarray(ci, dtype=np.int32)]
        return np.asarray(self._d.get("sample", []), dtype=np.uint64)

    @property
    def expression(self) -> np.ndarray:
        return np.asarray(self._d.get("expression", []), dtype=np.float32)

    # Allow dict-style access for backward-compatible code that does result['cell']
    def __getitem__(self, key: str):
        if key == "cell":
            return self.cells.tolist()
        if key == "sample":
            return self.samples.tolist()
        if key == "expression":
            return self.expression.tolist()
        if key == "expression_data":
            return self.expression_data
        if key == "ncells":
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
        self._data = data  # full top-level response dict (job_id, status, …)
        r = data.get("results", data)
        self._compact = isinstance(r, dict) and r.get("_format") == "compact_v2"
        self._gc = r.get("_global_cells") if self._compact else None
        self._gs = r.get("_global_samples") if self._compact else None
        self._r = r
        self.sample_metadata = r.get("_sample_metadata", {})

    # Allow dict-style field access for top-level metadata (status, job_id, …)
    def get(self, key, default=None):
        return self._data.get(key, default)

    @property
    def genes(self) -> list:
        return [k for k, v in self._r.items() if not k.startswith("_") and isinstance(v, dict)]

    def __getitem__(self, key):
        # Gene access (primary use case)
        if key in self._r and not key.startswith("_") and isinstance(self._r.get(key), dict):
            return GeneResult(self._r[key], self._gc, self._gs)
        # Top-level metadata access (job_id, status, …)
        return self._data[key]

    def __contains__(self, gene: str) -> bool:
        return gene in self._r and not gene.startswith("_")

    def items(self):
        for g in self.genes:
            yield g, self[g]

    def __iter__(self):
        return iter(self.genes)


def _reconstruct_compact_v2(data: dict) -> dict:
    """Reconstruct compact_v2 search result to dense format in-place.

    Kept for backward compatibility; prefer SearchResults for new code.
    """
    if not isinstance(data, dict) or data.get("_format") != "compact_v2":
        return data

    global_cells = data.get("_global_cells", [])
    global_samples = data.get("_global_samples", [])
    results = data.get("results", {})

    for gene, gdata in results.items():
        if gene.startswith("_") or not isinstance(gdata, dict):
            continue
        ci = gdata.pop("_ci", None)
        if ci:
            gdata["cell"] = [global_cells[i] for i in ci]
            gdata["sample"] = [global_samples[i] for i in ci]
        else:
            gdata["cell"] = []
            gdata["sample"] = []
        ed = gdata.get("expression_data")
        if isinstance(ed, dict) and ed.get("_fmt") == "col":
            gdata["expression_data"] = _expr_data_from_columnar(ed)

    data.pop("_format", None)
    data.pop("_global_cells", None)
    data.pop("_global_samples", None)
    return data


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class MalvaClient:
    """
    Main client for interacting with the Malva API

    Provides methods for searching, retrieving samples, and downloading data files.

    The client prefers canonical API paths and automatically falls back to the
    historical paths when the server does not expose them. Detection is a
    one-time authenticated probe of ``GET /api/account/quota`` whose result is
    cached per client; a 200 selects canonical mode and a 404 selects legacy
    mode. A 401 raises :class:`AuthenticationError`.

    Set ``MALVA_ENDPOINT_STYLE`` to ``canonical``, ``legacy`` or ``auto`` (the
    default) to override detection. The resolved style is exposed through the
    read-only :attr:`endpoint_style` property, and
    :meth:`reset_endpoint_style` clears a cached ``auto`` result.
    """

    def __init__(
        self,
        base_url: str = None,
        api_token: str = None,
        timeout: int = 300,
        verify_ssl: bool = True,
        workspace: str = None,
        dataset_id: str = None,
        provider_id: str = None,
    ):
        """
        Initialize the Malva client

        Args:
            base_url: Base URL of the Malva API
            api_token: API token for authentication
            timeout: Request timeout in seconds
            verify_ssl: Whether to verify SSL certificates
            workspace: Active ``project:index`` workspace (optional)
            dataset_id: Dataset/index used when a method omits one (optional)
            provider_id: Default federation provider for provider-scoped calls
                (optional). Per-call ``provider=`` arguments take precedence.
        """
        # Load config when any connection or workspace value is missing.
        if (
            base_url is None
            or api_token is None
            or workspace is None
            or dataset_id is None
            or provider_id is None
        ):
            config = Config.load()

            def _config_str(name: str) -> Optional[str]:
                value = getattr(config, name, None)
                return value if isinstance(value, str) and value else None

            base_url = base_url or config.server_url
            api_token = api_token or config.api_token
            workspace = workspace or _config_str("workspace")
            dataset_id = dataset_id or _config_str("dataset_id")
            provider_id = provider_id or _config_str("provider_id")
            if verify_ssl is None:
                verify_ssl = config.verify_ssl

        if not base_url:
            raise ValueError("base_url must be provided either directly or via config/environment")

        self.base_url = base_url.rstrip("/")
        self.api_token = api_token
        self.timeout = timeout
        self.verify_ssl = verify_ssl
        self.workspace = workspace or None
        self.dataset_id = dataset_id or None
        self.provider_id = provider_id or None

        # ``None`` means auto-detection is still pending. An explicit
        # ``MALVA_ENDPOINT_STYLE`` value skips the network probe entirely.
        self._endpoint_style = self._configured_endpoint_style()

        self.session = requests.Session()
        self.session.verify = verify_ssl
        self.session.headers.update(
            {
                "X-Malva-Client-Version": _CLIENT_VERSION,
                "User-Agent": f"malva-client/{_CLIENT_VERSION}",
            }
        )

        if api_token:
            self.session.headers.update({"Authorization": f"Bearer {api_token}"})
        self._apply_workspace_headers()
        self._test_connection()

    @staticmethod
    def _configured_endpoint_style() -> Optional[str]:
        """Return the explicit endpoint style from the environment, if any.

        ``auto`` (the default) returns ``None`` so the first request performs
        the capability probe and caches the result.
        """
        override = os.environ.get("MALVA_ENDPOINT_STYLE", "auto").strip().lower()
        override = override or "auto"
        if override not in _ENDPOINT_STYLES:
            raise ValueError(
                "MALVA_ENDPOINT_STYLE must be one of 'auto', 'canonical' or 'legacy' "
                f"(got {override!r})"
            )
        return None if override == "auto" else override

    def _detect_endpoint_style(self) -> str:
        """Resolve and cache whether the server exposes canonical endpoints."""
        if self._endpoint_style in {"canonical", "legacy"}:
            return self._endpoint_style

        explicit = self._configured_endpoint_style()
        if explicit is not None:
            self._endpoint_style = explicit
            return explicit

        try:
            # Authenticated capability probe. A canonical server answers 200;
            # a legacy-only server answers 404. Authentication failures must
            # surface exactly like any other request.
            self._request("GET", _ENDPOINT_PATHS["quota"][0])
        except AuthenticationError:
            raise
        except MalvaAPIError as exc:
            if self._is_not_found_error(exc):
                self._endpoint_style = "legacy"
                return "legacy"
            raise
        self._endpoint_style = "canonical"
        return "canonical"

    def reset_endpoint_style(self) -> None:
        """Forget the cached style so the next request probes again.

        An explicit ``MALVA_ENDPOINT_STYLE=canonical|legacy`` override is kept;
        only an ``auto`` result is discarded.
        """
        self._endpoint_style = self._configured_endpoint_style()

    @property
    def endpoint_style(self) -> str:
        """Resolved endpoint style: ``"canonical"`` or ``"legacy"``."""
        return self._detect_endpoint_style()

    def _endpoint(self, operation: str, **values: Any) -> str:
        """Build the canonical or legacy path for one logical operation."""
        try:
            canonical, legacy = _ENDPOINT_PATHS[operation]
        except KeyError as exc:
            raise ValueError(f"Unknown endpoint operation: {operation}") from exc
        template = canonical if self._detect_endpoint_style() == "canonical" else legacy
        return template.format(**values)

    def _apply_workspace_headers(self) -> None:
        """Send the active workspace context as Explorer-compatible headers."""
        for header, value in (
            ("X-Malva-Workspace", self.workspace),
            ("X-Malva-Dataset", self.dataset_id),
            ("X-Malva-Provider", self.provider_id),
        ):
            if value:
                self.session.headers[header] = str(value)
            else:
                self.session.headers.pop(header, None)

    def set_workspace(
        self,
        workspace: str = None,
        dataset_id: str = None,
        provider_id: str = None,
    ) -> None:
        """Select the workspace/dataset/provider used by subsequent requests.

        Passing an empty string clears the corresponding value. Per-call
        ``provider=`` arguments still win over ``provider_id``.
        """
        if workspace is not None:
            self.workspace = workspace or None
        if dataset_id is not None:
            self.dataset_id = dataset_id or None
        if provider_id is not None:
            self.provider_id = provider_id or None
        self._apply_workspace_headers()

    def get_workspace(self) -> Dict[str, Any]:
        """Return the active Explorer workspace, provider and dataset list.

        Used to validate requested federation provider ids before submitting a
        job, and to discover the dataset ids configured for this account.
        """
        return self._request("GET", "/api/workspace")

    def _default_dataset_id(self, dataset_id: str = None) -> Optional[str]:
        return str(dataset_id) if dataset_id else self.dataset_id

    def _provider_params(self, provider: Optional[str] = None) -> Dict[str, str]:
        """Return the ``?provider=`` query parameter for a metadata call.

        An explicit ``provider`` argument wins over the client-level
        ``provider_id``. Passing ``provider=""`` explicitly disables provider
        scoping for that call.
        """
        selected = provider if provider is not None else self.provider_id
        if selected is None:
            return {}
        selected = str(selected).strip()
        return {"provider": selected} if selected else {}

    def _test_connection(self):
        """Test connection and authentication with better error handling"""
        try:
            endpoint_style = self._detect_endpoint_style()
            logger.debug(f"Using {endpoint_style} endpoint style")

            response = self._request("GET", self._endpoint("health"))
            logger.info(f"Connected to Malva API: {response.get('status', 'unknown')}")

            # Test authentication if token provided
            if self.api_token:
                try:
                    quota_response = self._request("GET", self._endpoint("quota"))
                    logger.info(
                        f"Authenticated successfully: {quota_response.get('account_type', 'unknown')} account"
                    )
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
                raise MalvaAPIError(
                    f"Could not connect to Malva API at {self.base_url}. Please check the URL and ensure the server is running."
                )
            raise

    def get_quota_status(self) -> Dict[str, Any]:
        """Get current quota status"""
        if not self.api_token:
            raise AuthenticationError("API token required to check quota status")
        return self._request("GET", self._endpoint("quota"))

    def is_authenticated(self) -> bool:
        """Check if the client is properly authenticated"""
        try:
            self.get_quota_status()
            return True
        except (AuthenticationError, MalvaAPIError):
            return False

    def _request(self, method: str, endpoint: str, **kwargs) -> Dict[str, Any]:
        """Make a request to the API with enhanced error handling"""
        url = urljoin(self.base_url, endpoint.lstrip("/"))

        try:
            response = self.session.request(method, url, timeout=self.timeout, **kwargs)

            logger.debug(f"{method} {url} -> {response.status_code}")
            logger.debug(
                f"Response content type: {response.headers.get('content-type', 'unknown')}"
            )

            if response.status_code == 404:
                content_type = response.headers.get("content-type", "")
                if "text/html" in content_type:
                    raise MalvaAPIError(
                        f"Endpoint not found: {endpoint}. Check if the Malva API is running at {self.base_url}"
                    )

                try:
                    error_data = response.json()
                    error_msg = error_data.get(
                        "detail", error_data.get("error", "Endpoint not found")
                    )
                    raise MalvaAPIError(f"Not found: {error_msg}")
                except (json.JSONDecodeError, ValueError):
                    raise MalvaAPIError(f"Endpoint not found: {endpoint}")

            # Handle a server-side minimum client version
            elif response.status_code == 426:
                try:
                    error_data = response.json()
                except (json.JSONDecodeError, ValueError):
                    error_data = {}
                raise ClientVersionError(
                    error_data.get(
                        "error",
                        f"Client version {_CLIENT_VERSION} is no longer supported; "
                        "upgrade malva_client",
                    ),
                    min_version=error_data.get("min_client_version"),
                    client_version=error_data.get("client_version", _CLIENT_VERSION),
                )

            # Handle quota exceeded
            elif response.status_code == 429:
                try:
                    error_data = response.json()
                    quota_info = {
                        "account_type": error_data.get("account_type"),
                        "quota_exceeded": error_data.get("quota_exceeded", True),
                    }
                    raise QuotaExceededError(
                        error_data.get("error", "Search quota exceeded"), quota_info
                    )
                except (json.JSONDecodeError, ValueError):
                    raise QuotaExceededError("Search quota exceeded")

            elif response.status_code in [401, 403]:
                try:
                    error_data = response.json()
                    error_msg = error_data.get(
                        "error", error_data.get("detail", "Authentication failed")
                    )
                except (json.JSONDecodeError, ValueError):
                    error_msg = "Authentication failed"
                raise AuthenticationError(error_msg)

            elif not response.ok:
                try:
                    error_data = response.json()
                    error_msg = error_data.get(
                        "detail", error_data.get("error", f"HTTP {response.status_code}")
                    )
                    raise MalvaAPIError(f"API error: {error_msg}")
                except (json.JSONDecodeError, ValueError):
                    body = response.text[:200] if response.text else "(empty body)"
                    raise MalvaAPIError(f"HTTP {response.status_code}: {body}")

            try:
                result = response.json()
                # Return a lazy SearchResults wrapper for compact_v2 to avoid
                # O(N_cells × N_genes) list expansion at parse time.
                if (
                    isinstance(result, dict)
                    and result.get("results", {}).get("_format") == "compact_v2"
                ):
                    return SearchResults(result)
                return result
            except (json.JSONDecodeError, ValueError) as e:
                content_type = response.headers.get("content-type", "unknown")
                body = response.text[:200] if response.text else "(empty body)"
                raise MalvaAPIError(
                    f"Server returned non-JSON response (HTTP {response.status_code}, "
                    f"{content_type}): {body}"
                )

        except requests.exceptions.Timeout:
            raise MalvaAPIError(f"Request timed out after {self.timeout} seconds")
        except requests.exceptions.ConnectionError:
            raise MalvaAPIError(
                f"Could not connect to {self.base_url}. Check if the server is running."
            )
        except requests.exceptions.RequestException as e:
            raise MalvaAPIError(f"Request failed: {e}")

    def _build_search_payload(
        self,
        query: str,
        wait_for_completion: bool,
        max_wait: int,
        poll_interval: int,
        aggregate_expression: bool,
        stranded: Optional[bool],
        min_kmer_presence: int,
        max_kmer_presence: int,
    ) -> dict:
        """Build the POST /search payload with correctly mapped parameters.

        Always submits the job asynchronously (wait_for_completion=False at
        the HTTP level) to avoid holding a long-lived connection that would
        hit nginx proxy_read_timeout on slow queries.  The client-level
        wait_for_completion flag is handled by polling after submission.
        """
        payload = {
            "query": query,
            "wait_for_completion": False,  # never block the HTTP POST
            "aggregate_expression": aggregate_expression,
            "min_kmer_presence": min_kmer_presence,
            "max_kmer_presence": max_kmer_presence,
        }
        # stranded=True  → forward strand only  → unstranded_mode=False
        # stranded=False → both strands         → unstranded_mode=True
        # stranded=None  → use server default   → omit (server defaults to forward only)
        if stranded is True:
            payload["unstranded_mode"] = False
        elif stranded is False:
            payload["unstranded_mode"] = True
        if self.dataset_id:
            payload["dataset_id"] = self.dataset_id
        return payload

    def _get_job_status_lightweight(self, job_id: str) -> Dict[str, Any]:
        """Fetch job status without forcing the server to return result data."""
        try:
            return self._request("GET", self._endpoint("search_job_status", job_id=job_id))
        except MalvaAPIError as exc:
            if self._is_not_found_error(exc):
                return self._request("GET", self._endpoint("search_job", job_id=job_id))
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

    def _poll_until_complete(
        self, job_id: str, poll_interval: int, max_wait: int, context: str = "search"
    ) -> SearchResult:
        """Poll GET /search/<job_id> until completed, return lazy SearchResult."""
        start_time = time.time()
        while time.time() - start_time < max_wait:
            try:
                resp = self._get_job_status_lightweight(job_id)
                status = resp.get("status")
                logger.info(f"{context} status: {status}")
                if status == "completed":
                    logger.info(f"{context} completed: {job_id}")
                    return SearchResult({"job_id": job_id, "status": "completed"}, self)
                elif status == "error":
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

    def search(
        self,
        query: str,
        wait_for_completion: bool = True,
        poll_interval: int = 2,
        max_wait: int = 300,
        aggregate_expression: bool = True,
        stranded: Optional[bool] = None,
        min_kmer_presence: int = 0,
        max_kmer_presence: int = 100000,
    ) -> SearchResult:
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
            query,
            wait_for_completion,
            max_wait,
            poll_interval,
            aggregate_expression,
            stranded,
            min_kmer_presence,
            max_kmer_presence,
        )
        response = self._request("POST", self._endpoint("search"), json=payload)

        job_id = response.get("job_id")
        if not job_id:
            raise MalvaAPIError("No job_id received from server")

        status = response.get("status")
        logger.info(f"Search submitted: job_id={job_id}, status={status}")

        if status == "error":
            raise MalvaAPIError(f"Search failed: {response.get('error', 'unknown')}")

        if status == "completed":
            return SearchResult(
                {
                    "job_id": job_id,
                    "status": "completed",
                    "aggregate_expression": aggregate_expression,
                },
                self,
            )

        if not wait_for_completion:
            return SearchResult(
                {
                    "job_id": job_id,
                    "status": "pending",
                    "aggregate_expression": aggregate_expression,
                },
                self,
            )

        result = self._poll_until_complete(job_id, poll_interval, max_wait, "search")
        result.raw_data["aggregate_expression"] = aggregate_expression
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

    def _request_first_available(
        self, method: str, endpoints: List[str], **kwargs
    ) -> Tuple[Dict[str, Any], str]:
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

    @staticmethod
    def _emit_stage(
        callback: Optional[Callable[[str, Dict[str, Any]], None]], stage: str, **details
    ) -> None:
        payload = {"stage": stage}
        payload.update(details)
        msg = f"Coverage stage: {stage}"
        if details:
            compact = ", ".join(f"{k}={v}" for k, v in details.items() if v is not None)
            msg = f"{msg} ({compact})"
        logger.info(msg)
        if callback:
            callback(stage, payload)

    def _request_json_with_transfer_timing(
        self,
        method: str,
        endpoint: str,
        stage_callback: Optional[Callable[[str, Dict[str, Any]], None]] = None,
        stage_prefix: str = "coverage_result",
        **kwargs,
    ) -> Dict[str, Any]:
        """Request JSON while reporting server-generation, download, and parse stages.

        With requests, time-to-first-byte approximates server-side result
        generation. The streamed body read then measures the actual transfer
        time, which is what users notice for large coverage responses.
        """
        if endpoint.startswith(("http://", "https://")):
            url = endpoint
        else:
            url = urljoin(self.base_url, endpoint.lstrip("/"))

        request_started = time.time()
        try:
            response = self.session.request(
                method,
                url,
                timeout=self.timeout,
                stream=True,
                **kwargs,
            )
        except requests.exceptions.Timeout:
            raise MalvaAPIError(f"Request timed out after {self.timeout} seconds")
        except requests.exceptions.ConnectionError:
            raise MalvaAPIError(
                f"Could not connect to {self.base_url}. Check if the server is running."
            )
        except requests.exceptions.RequestException as exc:
            raise MalvaAPIError(f"Request failed: {exc}")

        first_byte_elapsed = time.time() - request_started
        self._emit_stage(
            stage_callback,
            f"{stage_prefix}_response_started",
            elapsed_seconds=round(first_byte_elapsed, 3),
            status_code=response.status_code,
        )

        if response.status_code == 404:
            response.close()
            raise MalvaAPIError(f"Endpoint not found: {endpoint}")
        if response.status_code in [401, 403]:
            response.close()
            raise AuthenticationError("Authentication failed")
        if response.status_code == 426:
            try:
                error_data = response.json()
            except (json.JSONDecodeError, ValueError):
                error_data = {}
            response.close()
            raise ClientVersionError(
                error_data.get(
                    "error",
                    f"Client version {_CLIENT_VERSION} is no longer supported; "
                    "upgrade malva_client",
                ),
                min_version=error_data.get("min_client_version"),
                client_version=error_data.get("client_version", _CLIENT_VERSION),
            )
        if response.status_code == 429:
            response.close()
            raise QuotaExceededError("Search quota exceeded")
        if not response.ok:
            try:
                error_data = response.json()
                msg = error_data.get(
                    "detail", error_data.get("error", f"HTTP {response.status_code}")
                )
            except (json.JSONDecodeError, ValueError):
                msg = response.text[:200] if response.text else f"HTTP {response.status_code}"
            finally:
                response.close()
            raise MalvaAPIError(f"API error: {msg}")

        chunks = []
        bytes_read = 0
        download_started = time.time()
        last_report = download_started
        try:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if not chunk:
                    continue
                chunks.append(chunk)
                bytes_read += len(chunk)
                now = time.time()
                if now - last_report >= 5:
                    self._emit_stage(
                        stage_callback,
                        f"{stage_prefix}_downloading",
                        bytes=bytes_read,
                        elapsed_seconds=round(now - download_started, 3),
                    )
                    last_report = now
        finally:
            encoding = response.encoding or "utf-8"
            response.close()

        download_elapsed = time.time() - download_started
        self._emit_stage(
            stage_callback,
            f"{stage_prefix}_downloaded",
            bytes=bytes_read,
            elapsed_seconds=round(download_elapsed, 3),
        )

        parse_started = time.time()
        body = b"".join(chunks)
        try:
            result = json.loads(body.decode(encoding))
        except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
            preview = body[:200].decode(encoding, errors="replace") if body else "(empty body)"
            raise MalvaAPIError(f"Server returned invalid JSON response: {preview}") from exc
        self._emit_stage(
            stage_callback,
            f"{stage_prefix}_parsed",
            elapsed_seconds=round(time.time() - parse_started, 3),
        )
        return result

    def _request_msgpack_with_transfer_timing(
        self,
        method: str,
        endpoint: str,
        stage_callback: Optional[Callable[[str, Dict[str, Any]], None]] = None,
        stage_prefix: str = "coverage_result",
        **kwargs,
    ) -> Dict[str, Any]:
        """Request msgpack while reporting response, download, and parse stages."""
        if endpoint.startswith(("http://", "https://")):
            url = endpoint
        else:
            url = urljoin(self.base_url, endpoint.lstrip("/"))

        request_started = time.time()
        try:
            response = self.session.request(
                method,
                url,
                timeout=self.timeout,
                stream=True,
                **kwargs,
            )
        except requests.exceptions.Timeout:
            raise MalvaAPIError(f"Request timed out after {self.timeout} seconds")
        except requests.exceptions.ConnectionError:
            raise MalvaAPIError(
                f"Could not connect to {self.base_url}. Check if the server is running."
            )
        except requests.exceptions.RequestException as exc:
            raise MalvaAPIError(f"Request failed: {exc}")

        first_byte_elapsed = time.time() - request_started
        self._emit_stage(
            stage_callback,
            f"{stage_prefix}_response_started",
            elapsed_seconds=round(first_byte_elapsed, 3),
            status_code=response.status_code,
            encoding="msgpack",
        )

        if response.status_code == 404:
            response.close()
            raise MalvaAPIError(f"Endpoint not found: {endpoint}")
        if response.status_code in [401, 403]:
            response.close()
            raise AuthenticationError("Authentication failed")
        if response.status_code == 429:
            response.close()
            raise QuotaExceededError("Search quota exceeded")
        if not response.ok:
            try:
                error_data = response.json()
                msg = error_data.get(
                    "detail", error_data.get("error", f"HTTP {response.status_code}")
                )
            except (json.JSONDecodeError, ValueError):
                msg = response.text[:200] if response.text else f"HTTP {response.status_code}"
            finally:
                response.close()
            raise MalvaAPIError(f"API error: {msg}")

        chunks = []
        bytes_read = 0
        download_started = time.time()
        last_report = download_started
        try:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if not chunk:
                    continue
                chunks.append(chunk)
                bytes_read += len(chunk)
                now = time.time()
                if now - last_report >= 5:
                    self._emit_stage(
                        stage_callback,
                        f"{stage_prefix}_downloading",
                        bytes=bytes_read,
                        elapsed_seconds=round(now - download_started, 3),
                    )
                    last_report = now
        finally:
            response.close()

        self._emit_stage(
            stage_callback,
            f"{stage_prefix}_downloaded",
            bytes=bytes_read,
            elapsed_seconds=round(time.time() - download_started, 3),
            encoding="msgpack",
        )

        parse_started = time.time()
        body = b"".join(chunks)
        try:
            result = msgpack.unpackb(body, raw=False, strict_map_key=False)
        except Exception as exc:
            raise MalvaAPIError("Server returned invalid msgpack response") from exc
        self._emit_stage(
            stage_callback,
            f"{stage_prefix}_parsed",
            elapsed_seconds=round(time.time() - parse_started, 3),
            encoding="msgpack",
        )
        return result

    def _request_raw(self, method: str, endpoint: str, **kwargs) -> requests.Response:
        """Make a request and return the raw response for binary endpoints."""
        if endpoint.startswith(("http://", "https://")):
            url = endpoint
        else:
            url = urljoin(self.base_url, endpoint.lstrip("/"))
        try:
            response = self.session.request(method, url, timeout=self.timeout, **kwargs)
        except requests.exceptions.Timeout:
            raise MalvaAPIError(f"Request timed out after {self.timeout} seconds")
        except requests.exceptions.ConnectionError:
            raise MalvaAPIError(
                f"Could not connect to {self.base_url}. Check if the server is running."
            )
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
                msg = error_data.get(
                    "detail", error_data.get("error", f"HTTP {response.status_code}")
                )
            except (json.JSONDecodeError, ValueError):
                msg = response.text[:200] if response.text else f"HTTP {response.status_code}"
            raise MalvaAPIError(f"API error: {msg}")
        return response

    @staticmethod
    def _normalise_id_list(
        values: Optional[Union[int, str, List[Union[int, str]]]],
    ) -> Optional[List[int]]:
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

    def _federated_shard_selection(
        self, frame: pd.DataFrame, provider: Optional[str] = None
    ) -> Tuple[str, Optional[str]]:
        """Select a completed shard's search job and owning provider id.

        The federated frame records each provider's shard under
        ``frame.attrs["shards"]``. ``provider`` may match either the shard's
        ``provider_id`` or its human-readable ``provider_name``.
        """
        shards = frame.attrs.get("shards")
        if not shards:
            raise MalvaAPIError(
                "Federated search result has no shard metadata "
                "(df.attrs['shards'] is missing or empty)"
            )

        requested = str(provider) if provider is not None else None
        matched_provider = False
        for shard in shards:
            if not isinstance(shard, dict):
                continue
            if requested is not None:
                identifiers = {
                    str(shard.get("provider_id")),
                    str(shard.get("provider_name")),
                }
                if requested not in identifiers:
                    continue
                matched_provider = True
            if shard.get("status") == "completed" and shard.get("search_job_id"):
                shard_provider_id = shard.get("provider_id")
                return (
                    str(shard["search_job_id"]),
                    str(shard_provider_id) if shard_provider_id is not None else None,
                )
            if requested is not None:
                raise MalvaAPIError(
                    f"Federated shard {requested!r} has no completed search job "
                    f"(status={shard.get('status')!r}, "
                    f"search_job_id={shard.get('search_job_id')!r})"
                )

        if requested is not None and not matched_provider:
            raise MalvaAPIError(f"No federated shard matches provider {requested!r}")
        raise MalvaAPIError(
            "Federated search result has no completed shard with a search_job_id; "
            "cannot retrieve cells"
        )

    def _search_job_from_federated_frame(
        self, frame: pd.DataFrame, provider: Optional[str] = None
    ) -> str:
        """Return only the search job for a federated frame shard."""
        return self._federated_shard_selection(frame, provider=provider)[0]

    def _resolve_search_job(
        self,
        search_job: Union[str, SearchResult, Dict[str, Any]],
        poll_interval: int,
        max_wait: int,
    ) -> Tuple[str, Optional[SearchResult]]:
        if isinstance(search_job, SearchResult):
            aggregate_expression = (
                search_job.raw_data.get("aggregate_expression")
                if isinstance(search_job.raw_data, dict)
                else None
            )
            if search_job.status in {"pending", "running"}:
                search_job = self.wait_for_job(
                    search_job.job_id,
                    poll_interval=poll_interval,
                    max_wait=max_wait,
                )
                if aggregate_expression is not None:
                    search_job.raw_data["aggregate_expression"] = aggregate_expression
            return search_job.job_id, search_job
        if isinstance(search_job, dict):
            job_id = search_job.get("job_id")
            if not job_id:
                raise ValueError("search_job dictionary must include job_id")
            return str(job_id), SearchResult(search_job, self)
        if isinstance(search_job, str):
            if not search_job.strip():
                raise ValueError("search_job must be a SearchResult or a non-empty job_id")
            return search_job.strip(), None
        raise TypeError(
            "search_job must be a SearchResult, job_id string, or job response dictionary"
        )

    def _sample_metadata_frame(
        self, sample_ids: List[int], provider: Optional[str] = None
    ) -> pd.DataFrame:
        if not sample_ids:
            return pd.DataFrame()
        try:
            metadata = self.get_sample_metadata(sample_ids, provider=provider)
        except Exception as exc:
            logger.warning(f"Sample metadata lookup failed during cell retrieval: {exc}")
            return pd.DataFrame()

        rows = []
        for sample_id, values in metadata.items():
            if isinstance(values, dict):
                row = {"sample_id": int(sample_id)}
                row.update(values)
                rows.append(row)
        return pd.DataFrame(rows) if rows else pd.DataFrame()

    def _cell_barcodes_frame(self, cells: pd.DataFrame, batch_size: int = 100000) -> pd.DataFrame:
        if cells.empty:
            return pd.DataFrame()
        rows = []
        for start in range(0, len(cells), batch_size):
            batch = cells.iloc[start : start + batch_size]
            payload = {
                "sample_ids": batch["sample_id"].astype(int).tolist(),
                "cell_ids": batch["cell_id"].astype(int).tolist(),
            }
            try:
                response = self._request("POST", "/api/cells/barcodes", json=payload)
            except Exception as exc:
                logger.warning(f"Barcode lookup failed during cell retrieval: {exc}")
                return pd.DataFrame()
            barcodes = response.get("barcodes", [])
            for row_index, sample_id, cell_id, barcode in zip(
                batch["row_index"].tolist(),
                payload["sample_ids"],
                payload["cell_ids"],
                barcodes,
            ):
                rows.append(
                    {
                        "row_index": int(row_index),
                        "sample_id": int(sample_id),
                        "cell_id": int(cell_id),
                        "barcode": barcode,
                    }
                )
        return pd.DataFrame(rows) if rows else pd.DataFrame()

    def _result_from_entry_rows(
        self,
        job_id: str,
        entries: pd.DataFrame,
        features: pd.DataFrame,
        source: str,
        include_barcodes: bool,
        include_sample_metadata: bool,
        provider_id: Optional[str] = None,
    ) -> CellExpressionMatrixResult:
        if entries.empty:
            cells = pd.DataFrame(columns=["row_index", "sample_id", "cell_id"])
            matrix_entries = pd.DataFrame(columns=["row_index", "feature_index", "value"])
        else:
            cells = (
                entries[["sample_id", "cell_id"]]
                .drop_duplicates()
                .sort_values(["sample_id", "cell_id"])
                .reset_index(drop=True)
            )
            cells.insert(0, "row_index", np.arange(1, len(cells) + 1, dtype=np.int64))
            matrix_entries = entries.merge(cells, on=["sample_id", "cell_id"], how="left")
            matrix_entries = matrix_entries[["row_index", "feature_index", "value"]]
            matrix_entries = matrix_entries.sort_values(["feature_index", "row_index"]).reset_index(
                drop=True
            )

        sample_ids = (
            cells["sample_id"].astype(int).drop_duplicates().tolist() if not cells.empty else []
        )
        sample_metadata = (
            self._sample_metadata_frame(sample_ids, provider=provider_id)
            if include_sample_metadata
            else pd.DataFrame()
        )
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

    def _retrieve_cells_from_values_export_endpoint(
        self,
        job_id: str,
        features: Optional[List[str]],
        sample_ids: Optional[List[int]],
        include_barcodes: bool,
        include_sample_metadata: bool,
        poll_interval: int = 2,
        max_wait: int = 300,
    ) -> CellExpressionMatrixResult:
        """Retrieve raw per-cell values through the direct msgpack values export.

        This uses the same per-cell value resolver as the Explorer matrix export
        but asks the server for compact arrays instead of a ZIP/MatrixMarket
        archive. The Celery worker does any required aggregate=False rerun, so
        the API request path is not blocked.
        """
        if not features:
            raise MalvaAPIError(
                "No feature names are available for direct per-cell value retrieval"
            )

        items = [
            {"job_id": job_id, "feature": str(feature), "label": str(feature)}
            for feature in features
        ]
        queued = self._request(
            "POST",
            self._endpoint("cell_expression_export"),
            json={
                "items": items,
                "sample_ids": sample_ids or [],
                "include_barcodes": False,
                "output_format": "msgpack",
            },
        )
        export_id = queued.get("export_id")
        if not export_id:
            raise MalvaAPIError("No export_id received from cell-expression values endpoint")

        started = time.time()
        while time.time() - started < max_wait:
            status_payload = self._request(
                "GET", self._endpoint("cell_expression_export_job", export_id=export_id)
            )
            status = status_payload.get("status")
            if status == "completed":
                break
            if status == "error":
                raise MalvaAPIError(
                    f"Cell-expression values export failed: {status_payload.get('error', 'unknown')}"
                )
            delay = self._adaptive_poll_delay(time.time() - started, poll_interval)
            if delay > 0:
                time.sleep(delay)
        else:
            raise MalvaAPIError(f"Cell-expression values export timed out after {max_wait} seconds")

        response = self._request_raw(
            "GET", self._endpoint("cell_expression_export_download", export_id=export_id)
        )
        decoded = msgpack.unpackb(response.content, raw=False, strict_map_key=False)
        if decoded.get("format") != "cell_expression_values_v1":
            raise MalvaAPIError("Unexpected cell-expression values payload format")

        samples = np.frombuffer(decoded.get("s", b""), dtype=np.int64)
        cells = np.frombuffer(decoded.get("c", b""), dtype=np.int64)
        feature_indices = np.frombuffer(decoded.get("fi", b""), dtype=np.int32)
        values = np.frombuffer(decoded.get("v", b""), dtype=np.float32)
        n = min(
            len(samples),
            len(cells),
            len(feature_indices),
            len(values),
            int(decoded.get("n", 0) or 0),
        )
        if n == 0:
            raise MalvaAPIError("No per-cell expression values returned for this search job")

        entries = pd.DataFrame(
            {
                "sample_id": samples[:n],
                "cell_id": cells[:n],
                "feature_index": feature_indices[:n].astype(np.int64, copy=False),
                "value": values[:n],
            }
        )
        feature_rows = decoded.get("features") or []
        features_df = pd.DataFrame(feature_rows)
        if features_df.empty:
            features_df = pd.DataFrame(
                [
                    {
                        "feature_index": idx + 1,
                        "job_id": job_id,
                        "feature": feature,
                        "label": feature,
                        "source": "cell_expression_values",
                    }
                    for idx, feature in enumerate(features)
                ]
            )

        return self._result_from_entry_rows(
            job_id=job_id,
            entries=entries,
            features=features_df,
            source="cell_expression_values",
            include_barcodes=include_barcodes,
            include_sample_metadata=include_sample_metadata,
        )

    def _retrieve_cells_from_full_search_endpoint(
        self,
        job_id: str,
        features: Optional[List[str]],
        sample_ids: Optional[List[int]],
        include_barcodes: bool,
        include_sample_metadata: bool,
        provider_id: Optional[str] = None,
    ) -> CellExpressionMatrixResult:
        params = {"max_cells": 0}
        if provider_id is not None:
            params["provider"] = provider_id
        response = self._request("GET", self._endpoint("search_job", job_id=job_id), params=params)
        status = response.get("status") if hasattr(response, "get") else None
        if status and status != "completed":
            raise MalvaAPIError(f"Search job {job_id} is {status}")

        feature_filter = set(features) if features else None
        sample_filter = set(sample_ids) if sample_ids else None
        feature_rows = []
        entry_frames = []
        feature_index = 1

        if isinstance(response, dict):
            results = response.get("results", {}) if isinstance(response, dict) else {}
            iterator = (
                (feature, data)
                for feature, data in results.items()
                if isinstance(feature, str)
                and not feature.startswith("_")
                and isinstance(data, dict)
            )
        elif hasattr(response, "items"):
            iterator = response.items()
        else:
            iterator = iter(())

        for feature, data in iterator:
            if feature_filter is not None and feature not in feature_filter:
                continue

            if hasattr(data, "cells") and hasattr(data, "samples"):
                cells = np.asarray(data.cells, dtype=np.int64)
                samples = np.asarray(data.samples, dtype=np.int64)
                values = np.asarray(
                    getattr(data, "expression", np.ones(len(cells))), dtype=np.float32
                )
            else:
                cells = np.asarray(data.get("cell", []), dtype=np.int64)
                samples = np.asarray(data.get("sample", []), dtype=np.int64)
                values = np.asarray(data.get("expression", np.ones(len(cells))), dtype=np.float32)

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

            feature_rows.append(
                {
                    "feature_index": feature_index,
                    "job_id": job_id,
                    "feature": feature,
                    "label": feature,
                    "source": "search_full",
                }
            )
            entry_frames.append(
                pd.DataFrame(
                    {
                        "sample_id": samples,
                        "cell_id": cells,
                        "feature_index": feature_index,
                        "value": values,
                    }
                )
            )
            feature_index += 1

        if not entry_frames:
            raise MalvaAPIError("No cells matched the requested feature or sample filter")

        return self._result_from_entry_rows(
            job_id=job_id,
            entries=pd.concat(entry_frames, ignore_index=True),
            features=pd.DataFrame(feature_rows),
            source="search_full",
            include_barcodes=include_barcodes,
            include_sample_metadata=include_sample_metadata,
            provider_id=provider_id,
        )

    def _retrieve_cells_from_cell_ids_endpoint(
        self,
        job_id: str,
        features: Optional[List[str]],
        sample_ids: Optional[List[int]],
        include_barcodes: bool,
        include_sample_metadata: bool,
        provider_id: Optional[str] = None,
    ) -> CellExpressionMatrixResult:
        params = {"provider": provider_id} if provider_id is not None else {}
        response = self._request(
            "GET",
            self._endpoint("search_job_cell_ids", job_id=job_id),
            params=params,
        )
        results = response.get("results", {}) if isinstance(response, dict) else {}
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
            cells = np.asarray(data.get("cell", []), dtype=np.int64)
            samples = np.asarray(data.get("sample", []), dtype=np.int64)
            n = min(len(cells), len(samples))
            if n == 0:
                continue
            cells = cells[:n]
            samples = samples[:n]
            values_raw = data.get("expression")
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

            feature_rows.append(
                {
                    "feature_index": feature_index,
                    "job_id": job_id,
                    "feature": feature,
                    "label": feature,
                    "source": "search_cell_ids",
                }
            )
            entry_frames.append(
                pd.DataFrame(
                    {
                        "sample_id": samples,
                        "cell_id": cells,
                        "feature_index": feature_index,
                        "value": values,
                    }
                )
            )
            feature_index += 1

        if not entry_frames:
            raise MalvaAPIError("No cells matched the requested feature or sample filter")

        return self._result_from_entry_rows(
            job_id=job_id,
            entries=pd.concat(entry_frames, ignore_index=True),
            features=pd.DataFrame(feature_rows),
            source="search_cell_ids",
            include_barcodes=include_barcodes,
            include_sample_metadata=include_sample_metadata,
            provider_id=provider_id,
        )

    def _retrieve_cells_from_global_ids_endpoint(
        self,
        job_id: str,
        features: Optional[List[str]],
        sample_ids: Optional[List[int]],
        include_barcodes: bool,
        include_sample_metadata: bool,
        provider_id: Optional[str] = None,
    ) -> CellExpressionMatrixResult:
        params = {}
        if sample_ids:
            params["sample_ids"] = ",".join(str(s) for s in sample_ids)
        if provider_id is not None:
            params["provider"] = provider_id
        response = self._request_raw(
            "GET",
            self._endpoint("search_job_global_cell_ids", job_id=job_id),
            params=params,
        )
        decoded = msgpack.unpackb(response.content, raw=False, strict_map_key=False)
        cells = np.frombuffer(decoded.get("c", b""), dtype=np.int64)
        samples = np.frombuffer(decoded.get("s", b""), dtype=np.int64)
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

        label = features[0] if features and len(features) == 1 else "positive_cells"
        entries = pd.DataFrame(
            {
                "sample_id": samples,
                "cell_id": cells,
                "feature_index": 1,
                "value": np.ones(len(cells), dtype=np.float32),
            }
        )
        feature_df = pd.DataFrame(
            [
                {
                    "feature_index": 1,
                    "job_id": job_id,
                    "feature": label,
                    "label": label,
                    "source": "global_cell_ids",
                }
            ]
        )
        return self._result_from_entry_rows(
            job_id=job_id,
            entries=entries,
            features=feature_df,
            source="global_cell_ids",
            include_barcodes=include_barcodes,
            include_sample_metadata=include_sample_metadata,
            provider_id=provider_id,
        )

    def _retrieve_cells_from_explorer_endpoint(
        self,
        job_id: str,
        features: Optional[List[str]],
        sample_ids: Optional[List[int]],
        include_barcodes: bool,
        include_sample_metadata: bool,
        provider_id: Optional[str] = None,
    ) -> CellExpressionMatrixResult:
        params = {}
        if sample_ids:
            params["sample_ids"] = ",".join(str(s) for s in sample_ids)
        if include_barcodes:
            params["barcodes"] = "1"
        if provider_id is not None:
            params["provider"] = provider_id

        response = self._request(
            "GET", self._endpoint("search_job_cells", job_id=job_id), params=params
        )
        samples_payload = response.get("samples", {}) if isinstance(response, dict) else {}
        entries = []
        barcode_rows = []
        for sample_id, payload in samples_payload.items():
            if not isinstance(payload, dict):
                continue
            cell_ids = payload.get("cell_ids") or []
            barcodes = payload.get("barcodes") or []
            for idx, cell_id in enumerate(cell_ids):
                entries.append(
                    {
                        "sample_id": int(sample_id),
                        "cell_id": int(cell_id),
                        "feature_index": 1,
                        "value": 1.0,
                    }
                )
                if include_barcodes and idx < len(barcodes):
                    barcode_rows.append(
                        {
                            "sample_id": int(sample_id),
                            "cell_id": int(cell_id),
                            "barcode": barcodes[idx],
                        }
                    )

        if not entries:
            raise MalvaAPIError("No cells returned for this search job")

        label = features[0] if features and len(features) == 1 else "positive_cells"
        result = self._result_from_entry_rows(
            job_id=job_id,
            entries=pd.DataFrame(entries),
            features=pd.DataFrame(
                [
                    {
                        "feature_index": 1,
                        "job_id": job_id,
                        "feature": label,
                        "label": label,
                        "source": "expression_cells",
                    }
                ]
            ),
            source="expression_cells",
            include_barcodes=False,
            include_sample_metadata=include_sample_metadata,
            provider_id=provider_id,
        )
        if barcode_rows:
            barcode_df = pd.DataFrame(barcode_rows)
            result._barcodes_df = result.cells.merge(
                barcode_df,
                on=["sample_id", "cell_id"],
                how="inner",
            )
        return result

    def stream_cell_expression(
        self,
        sample_ids: Optional[Union[int, str, List[Union[int, str]]]] = None,
        gene_ids: Optional[Union[str, List[str]]] = None,
        dataset_id: Optional[str] = None,
        max_batch_bytes: int = 4 * 1024 * 1024,
        chunk_size: int = 64 * 1024,
        cells: Optional[Iterable[Any]] = None,
    ) -> CellExpressionStream:
        """Stream sparse per-cell expression batches for selected samples.

        By default this returns a :class:`~malva_client.models.CellExpressionStream`
        that yields one pandas ``DataFrame`` per server batch, exactly like the
        previous generator API::

            for batch in client.stream_cell_expression(
                sample_ids=[30064771072],
                gene_ids=["CD3D", "CD3E", "CD247"],
                dataset_id="Human_Cell_Atlas",
            ):
                print(batch.head())

        DataFrames contain ``sample_id``, ``cell_id``, ``store_index``, ``gene``
        and ``value``. Only non-zero values are returned; absent
        ``(sample_id, cell_id, gene)`` combinations are zero. Passing
        ``gene_ids`` is strongly recommended: without it the server must scan
        every non-zero entry in the selected petal stores.

        The wrapper also provides ``to_anndata()``, ``write_h5ad(path)`` and
        ``iter_anndata()`` for standard AnnData workflows.
        """
        dataset_id = self._default_dataset_id(dataset_id)
        gene_list = None
        if gene_ids is not None:
            if isinstance(gene_ids, str):
                gene_list = [gene_ids]
            else:
                gene_list = [str(gene).strip() for gene in gene_ids if str(gene).strip()]
        if cells is not None:
            if sample_ids is not None:
                raise ValueError("Pass either sample_ids or cells, not both")
            known_samples = None
            if isinstance(cells, (list, tuple)):
                known_samples = sorted(
                    {
                        int(cell.get("sample_id"))
                        if isinstance(cell, dict)
                        else int(cell[0])
                        if isinstance(cell, (tuple, list)) and len(cell) == 2
                        else int(getattr(cell, "sample_id"))
                        for cell in cells
                    }
                )
        else:
            if sample_ids is None:
                raise ValueError("Pass sample_ids or cells")
            if isinstance(sample_ids, (int, str)):
                sample_ids = [sample_ids]
            known_samples = [int(sample_id) for sample_id in sample_ids]
            if not known_samples:
                raise ValueError("sample_ids must contain at least one sample id")

        return CellExpressionStream(
            lambda: self._iter_cell_expression_batches(
                sample_ids=None if cells is not None else sample_ids,
                gene_ids=gene_list,
                dataset_id=dataset_id,
                max_batch_bytes=max_batch_bytes,
                chunk_size=chunk_size,
                cells=cells,
            ),
            client=self,
            dataset_id=dataset_id,
            sample_ids=known_samples,
            gene_ids=gene_list,
        )

    def _iter_cell_expression_batches(
        self,
        sample_ids: Optional[Union[int, str, List[Union[int, str]]]] = None,
        gene_ids: Optional[Union[str, List[str]]] = None,
        dataset_id: Optional[str] = None,
        max_batch_bytes: int = 4 * 1024 * 1024,
        chunk_size: int = 64 * 1024,
        cells: Optional[Iterable[Any]] = None,
    ):
        """Yield the raw pandas batches from the server frame stream."""
        if cells is not None:
            if sample_ids is not None:
                raise ValueError("Pass either sample_ids or cells, not both")
            cell_payload: List[Dict[str, int]] = []
            for cell in cells:
                if isinstance(cell, dict):
                    sample_id = cell.get("sample_id")
                    cell_id = cell.get("cell_id")
                elif isinstance(cell, (tuple, list)) and len(cell) == 2:
                    sample_id, cell_id = cell
                else:
                    sample_id = getattr(cell, "sample_id", None)
                    cell_id = getattr(cell, "cell_id", None)
                if sample_id is None or cell_id is None:
                    raise ValueError(
                        "Each cell must be (sample_id, cell_id), a dict, or an object with those fields"
                    )
                cell_payload.append({"sample_id": int(sample_id), "cell_id": int(cell_id)})
            if not cell_payload:
                raise ValueError("cells must contain at least one cell")
            payload: Dict[str, Any] = {
                "cells": cell_payload,
                "max_batch_bytes": int(max_batch_bytes),
            }
        else:
            if sample_ids is None:
                raise ValueError("Pass sample_ids or cells")
            if isinstance(sample_ids, (int, str)):
                sample_ids = [sample_ids]
            samples = [int(sample_id) for sample_id in sample_ids]
            if not samples:
                raise ValueError("sample_ids must contain at least one sample id")
            payload = {
                "sample_ids": samples,
                "max_batch_bytes": int(max_batch_bytes),
            }
        if dataset_id:
            payload["dataset_id"] = str(dataset_id)
        if gene_ids is not None:
            if isinstance(gene_ids, str):
                gene_ids = [gene_ids]
            genes = [str(gene).strip() for gene in gene_ids if str(gene).strip()]
            if genes:
                payload["gene_ids"] = genes

        url = urljoin(self.base_url + "/", "api/expression/cell-expression-stream")
        response = self.session.post(
            url,
            json=payload,
            stream=True,
            timeout=(10, self.timeout),
        )
        try:
            if response.status_code == 401:
                raise AuthenticationError("API token was rejected by the cell expression stream")
            if response.status_code == 403:
                raise AuthenticationError("API token is not allowed to stream cell expression")
            if not response.ok:
                detail = response.text[:500] if response.text else ""
                raise MalvaAPIError(
                    f"Cell expression stream failed ({response.status_code}): {detail}"
                )

            features_by_store: Dict[int, List[str]] = {}
            buffer = bytearray()
            for chunk in response.iter_content(chunk_size=chunk_size):
                if not chunk:
                    continue
                buffer.extend(chunk)
                while len(buffer) >= 8:
                    frame_size = struct.unpack(">Q", buffer[:8])[0]
                    if len(buffer) < 8 + frame_size:
                        break
                    frame = msgpack.unpackb(bytes(buffer[8 : 8 + frame_size]), raw=False)
                    del buffer[: 8 + frame_size]
                    frame_type = frame.get("type")
                    if frame_type == "features":
                        store_index = int(frame.get("store_index", 0))
                        features_by_store[store_index] = list(frame.get("features") or [])
                        continue
                    if frame_type != "data":
                        continue

                    count = int(frame.get("n") or 0)
                    if count == 0:
                        continue
                    store_index = int(frame.get("store_index", 0))
                    feature_index = np.frombuffer(frame["feature_index"], dtype="<i4")
                    features = features_by_store.get(store_index, [])
                    gene_names = [
                        features[int(index)]
                        if 0 <= int(index) < len(features)
                        else f"feature_{int(index)}"
                        for index in feature_index
                    ]
                    yield pd.DataFrame(
                        {
                            "sample_id": np.frombuffer(frame["sample_ids"], dtype="<i8"),
                            "cell_id": np.frombuffer(frame["cell_ids"], dtype="<i8"),
                            "store_index": np.full(count, store_index, dtype=np.int32),
                            "gene": gene_names,
                            "value": np.frombuffer(frame["values"], dtype="<f4"),
                        }
                    )
        finally:
            response.close()

    def retrieve_cells(
        self,
        search_job: Union[str, SearchResult, Dict[str, Any], pd.DataFrame],
        features: Optional[Union[str, List[str]]] = None,
        sample_ids: Optional[Union[int, str, List[Union[int, str]]]] = None,
        include_barcodes: bool = False,
        include_sample_metadata: bool = True,
        poll_interval: int = 2,
        max_wait: int = 300,
        provider: Optional[str] = None,
    ) -> CellExpressionMatrixResult:
        """
        Retrieve positive cells for an existing search job.

        Run :meth:`search` or :meth:`search_sequences` first, then pass the
        returned SearchResult or its job_id here. A DataFrame returned by
        :meth:`search_federated` is also accepted: the client reads
        ``frame.attrs["shards"]`` and retrieves cells from one completed shard.
        Pass ``provider`` (provider id or name) to choose a specific shard;
        otherwise the first completed shard with a search job is used. When a
        provider is selected, its id is forwarded as ``?provider=`` to every
        search-service call used for retrieval, and the local values export is
        skipped. This method does not start a new search and does not use the
        ZIP export workflow. The client first asks the search service for
        per-feature cell/sample/expression arrays so the ``value`` column
        contains raw per-cell expression/k-mer counts when available. The fast
        global cell-ID endpoint is used only as a fallback; it cannot carry
        per-cell expression values and therefore reports positive-cell
        indicators.

        """
        if isinstance(search_job, pd.DataFrame):
            job_id, shard_provider_id = self._federated_shard_selection(
                search_job, provider=provider
            )
            # Only an explicit provider selection scopes the search-service
            # calls; implicit selection keeps the historical local behavior.
            provider_id = shard_provider_id if provider is not None else None
            result = None
        else:
            job_id, result = self._resolve_search_job(search_job, poll_interval, max_wait)
            provider_id = str(provider) if provider is not None else None
        feature_list = self._normalise_feature_list(features)
        if feature_list is None and result is not None:
            try:
                raw_results = result.results
            except Exception:
                raw_results = (
                    result.raw_data.get("results", {}) if isinstance(result.raw_data, dict) else {}
                )
            feature_list = [
                key
                for key, value in raw_results.items()
                if isinstance(key, str) and not key.startswith("_") and isinstance(value, dict)
            ] or None
        sample_id_list = self._normalise_id_list(sample_ids)

        errors = []
        fetchers = (
            self._retrieve_cells_from_values_export_endpoint,
            self._retrieve_cells_from_full_search_endpoint,
            self._retrieve_cells_from_cell_ids_endpoint,
            self._retrieve_cells_from_global_ids_endpoint,
            self._retrieve_cells_from_explorer_endpoint,
        )

        for fetcher in fetchers:
            is_values_export = (
                getattr(fetcher, "__name__", "") == "_retrieve_cells_from_values_export_endpoint"
            )
            if is_values_export and provider_id is not None:
                # Exports only run against the local service; provider-scoped
                # retrieval must go through the search-service fallbacks.
                logger.debug("Skipping local values export for provider %s", provider_id)
                continue
            try:
                if is_values_export:
                    return fetcher(
                        job_id,
                        feature_list,
                        sample_id_list,
                        include_barcodes,
                        include_sample_metadata,
                        poll_interval=poll_interval,
                        max_wait=max_wait,
                    )
                return fetcher(
                    job_id,
                    feature_list,
                    sample_id_list,
                    include_barcodes,
                    include_sample_metadata,
                    provider_id=provider_id,
                )
            except MalvaAPIError as exc:
                message = str(exc)
                errors.append(message)
                timeout_message = "timed out" in message.lower() or "timeout" in message.lower()
                recoverable = (
                    self._is_not_found_error(exc)
                    or "No cell IDs returned" in message
                    or "No cells returned" in message
                    or "No cells matched" in message
                    or "No feature names are available" in message
                    or (is_values_export and timeout_message)
                )
                if not recoverable:
                    raise
                if is_values_export and timeout_message:
                    logger.warning(
                        "Cell-expression values export timed out; falling back to "
                        "the next cell retrieval endpoint: %s",
                        message,
                    )

        detail = "; ".join(errors[-3:]) if errors else "No retrieval endpoint was available"
        scope = f"search job {job_id}"
        if provider_id is not None:
            scope += f" for provider {provider_id!r}"
        raise MalvaAPIError(f"Could not retrieve cells for {scope}: {detail}")

    def search_sequences(
        self,
        sequences: Union[str, List[str]],
        wait_for_completion: bool = True,
        poll_interval: int = 2,
        max_wait: int = 300,
        aggregate_expression: bool = True,
        stranded: Optional[bool] = None,
        min_kmer_presence: int = 0,
        max_kmer_presence: int = 100000,
    ) -> SearchResult:
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
        valid_nucleotides = set("ATGCUN")

        # ── Single sequence ────────────────────────────────────────────────
        if isinstance(sequences, str):
            sequence = sequences
            if not sequence:
                raise ValueError("Sequence must not be empty")
            if len(sequence) > 500_000:
                raise ValueError("Sequence length cannot exceed 500,000 nucleotides")
            if not all(c.upper() in valid_nucleotides for c in sequence):
                raise ValueError(
                    "Sequence contains invalid nucleotides. Only A, T, G, C, U, N are allowed."
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
            "batch sequence search",
            wait_for_completion,
            max_wait,
            poll_interval,
            aggregate_expression,
            stranded,
            min_kmer_presence,
            max_kmer_presence,
        )
        payload["fasta_sequences"] = fasta_content

        logger.info(
            f"Submitting batch search: {len(sequences)} sequences, {total_nt:,} total nucleotides"
        )

        response = self._request("POST", self._endpoint("search"), json=payload)

        job_id = response.get("job_id")
        if not job_id:
            raise MalvaAPIError("No job_id received from server")

        status = response.get("status")
        logger.info(f"Batch search submitted: job_id={job_id}, status={status}")

        if status == "error":
            raise MalvaAPIError(f"Batch search failed: {response.get('error', 'unknown')}")

        if status == "completed":
            return SearchResult(
                {
                    "job_id": job_id,
                    "status": "completed",
                    "aggregate_expression": aggregate_expression,
                },
                self,
            )

        if not wait_for_completion:
            return SearchResult(
                {
                    "job_id": job_id,
                    "status": "pending",
                    "aggregate_expression": aggregate_expression,
                },
                self,
            )

        result = self._poll_until_complete(job_id, poll_interval, max_wait, "batch search")
        result.raw_data["aggregate_expression"] = aggregate_expression
        return result

    def search_genes(
        self,
        genes: List[str],
        wait_for_completion: bool = True,
        poll_interval: int = 2,
        max_wait: int = 300,
        aggregate_expression: bool = True,
        stranded: Optional[bool] = None,
        min_kmer_presence: int = 0,
        max_kmer_presence: int = 100000,
    ) -> SearchResult:
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
            gene_query,
            wait_for_completion,
            max_wait,
            poll_interval,
            aggregate_expression,
            stranded,
            min_kmer_presence,
            max_kmer_presence,
        )
        response = self._request("POST", self._endpoint("search"), json=payload)

        job_id = response.get("job_id")
        if not job_id:
            raise MalvaAPIError("No job_id received from server")

        status = response.get("status")
        logger.info(f"Gene search submitted: job_id={job_id}, status={status}")

        if status == "error":
            raise MalvaAPIError(f"Gene search failed: {response.get('error', 'unknown')}")

        if status == "completed":
            return SearchResult(
                {
                    "job_id": job_id,
                    "status": "completed",
                    "aggregate_expression": aggregate_expression,
                },
                self,
            )

        if not wait_for_completion:
            return SearchResult(
                {
                    "job_id": job_id,
                    "status": "pending",
                    "aggregate_expression": aggregate_expression,
                },
                self,
            )

        result = self._poll_until_complete(job_id, poll_interval, max_wait, "gene search")
        result.raw_data["aggregate_expression"] = aggregate_expression
        return result

    def get_samples(
        self,
        page: int = 1,
        page_size: int = 50,
        filters: Dict[str, Any] = None,
        search_query: str = "",
        provider: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Retrieve samples with optional filtering

        Args:
            page: Page number (1-based)
            page_size: Number of samples per page
            filters: Dictionary of filters (e.g., {'organ': 'brain', 'species': 'human'})
            search_query: Text search query
            provider: Federation provider id. Defaults to the client-level
                provider; an explicit value wins.

        Returns:
            Dictionary containing samples and pagination info
        """
        params = {"page": page, "page_size": page_size, "q": search_query}
        params.update(self._provider_params(provider))

        if filters:
            for key, value in filters.items():
                if isinstance(value, list):
                    params[key] = value
                else:
                    params[key] = [value] if value else []

        return self._request("GET", self._endpoint("samples"), params=params)

    def search_samples(
        self,
        query: str = "",
        organ: List[str] = None,
        species: List[str] = None,
        study: List[str] = None,
        disease: List[str] = None,
        provider: Optional[str] = None,
        **kwargs,
    ) -> pd.DataFrame:
        """
        Search samples and return as a DataFrame

        Args:
            query: Text search query
            organ: List of organs to filter by
            species: List of species to filter by
            study: List of studies to filter by
            disease: List of diseases to filter by
            provider: Federation provider id. Defaults to the client-level
                provider; an explicit value wins.
            **kwargs: Additional filters

        Returns:
            DataFrame containing sample metadata
        """
        filters = {}
        if organ:
            filters["organ"] = organ
        if species:
            filters["species"] = species
        if study:
            filters["study"] = study
        if disease:
            filters["disease"] = disease

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
                search_query=query,
                provider=provider,
            )

            samples = response.get("samples", [])
            if not samples:
                break

            all_samples.extend(samples)

            if not response.get("has_more", False):
                break

            page += 1

        return pd.DataFrame(all_samples)

    def get_studies(self, provider: Optional[str] = None) -> pd.DataFrame:
        """
        Get all studies with summary information.

        Some public providers return an empty list from
        ``/samples/api/studies`` even though their dataset hierarchy contains
        studies. In that case the list is derived from
        ``/samples/api/datasets/hierarchy`` instead; hierarchy-derived rows also
        carry a ``dataset_id`` column so they can be passed to
        :meth:`get_study_samples`.

        Args:
            provider: Federation provider id. Defaults to the client-level
                provider; an explicit value wins.

        Returns:
            DataFrame containing study information
        """
        params = self._provider_params(provider)
        response = self._request("GET", self._endpoint("studies"), params=params or None)
        studies = response.get("studies", [])
        if studies:
            return pd.DataFrame(studies)

        hierarchy_rows = self._studies_from_hierarchy(provider=provider)
        if hierarchy_rows:
            logger.warning(
                "Study summary endpoint returned no studies; using the dataset "
                "hierarchy fallback (%d studies).",
                len(hierarchy_rows),
            )
            return pd.DataFrame(hierarchy_rows)
        return pd.DataFrame(studies)

    def _studies_from_hierarchy(self, provider: Optional[str] = None) -> List[Dict[str, Any]]:
        """Flatten dataset → study entries from the hierarchy endpoint."""
        try:
            hierarchy = self.get_datasets_hierarchy(provider=provider)
        except MalvaAPIError as exc:
            logger.warning("Could not derive studies from the dataset hierarchy: %s", exc)
            return []

        rows: List[Dict[str, Any]] = []
        seen = set()
        for dataset in hierarchy.get("datasets") or []:
            if not isinstance(dataset, dict):
                continue
            dataset_id = dataset.get("dataset_id") or dataset.get("id")
            for study in dataset.get("studies") or []:
                if not isinstance(study, dict):
                    continue
                study_name = study.get("study_name") or study.get("name")
                key = (str(dataset_id) if dataset_id else None, str(study_name))
                if key in seen:
                    continue
                seen.add(key)
                row = dict(study)
                if dataset_id:
                    row["dataset_id"] = dataset_id
                rows.append(row)
        return rows

    def _normalise_sample_download_response(
        self, sample_uuid: str, response: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Normalize explorer and metadata-service sample download responses."""
        download_url = response.get("download_url")
        token = response.get("download_token")
        if not download_url and token:
            download_url = self._endpoint("download_token", token=token)
        if not download_url:
            raise MalvaAPIError("Download URL not provided by server")
        parsed_download = urlparse(download_url)
        if parsed_download.scheme and parsed_download.netloc:
            if parsed_download.path.startswith(("/api/download/", "/api/downloads/")):
                parsed_base = urlparse(self.base_url)
                download_url = urlunparse(
                    (
                        parsed_base.scheme,
                        parsed_base.netloc,
                        parsed_download.path,
                        "",
                        parsed_download.query,
                        parsed_download.fragment,
                    )
                )
        else:
            download_url = urljoin(self.base_url, download_url.lstrip("/"))

        expires_in = response.get("expires_in_seconds")
        if expires_in is None and response.get("expires_at") is not None:
            expires_in = max(0, int(float(response["expires_at"]) - time.time()))
        if expires_in is None:
            expires_in = 0

        return {
            "download_url": download_url,
            "filename": response.get("filename", f"{sample_uuid}.h5ad"),
            "file_size": response.get("file_size", 0),
            "expires_in_seconds": expires_in,
        }

    def _request_sample_download(
        self, sample_uuid: str, provider: Optional[str] = None
    ) -> Dict[str, Any]:
        """Request sample download metadata from either compatibility API shape."""
        params = self._provider_params(provider) or None
        try:
            response = self._request("GET", f"/api/samples/{sample_uuid}/download", params=params)
        except MalvaAPIError as exc:
            if not self._is_not_found_error(exc):
                raise
            response = self._request(
                "POST", f"/api/samples/{sample_uuid}/request-download", params=params
            )
        return self._normalise_sample_download_response(sample_uuid, response)

    def download_sample(
        self,
        sample_uuid: str,
        output_path: Optional[str] = None,
        chunk_size: int = 8192,
        show_progress: bool = True,
        provider: Optional[str] = None,
    ) -> Union[str, "ad.AnnData"]:
        """
        Download a sample's data file with enhanced error handling

        Args:
            sample_uuid: UUID of the sample to download
            output_path: Path to save the file. If None, loads into memory as AnnData
            chunk_size: Size of chunks to download (bytes)
            show_progress: Whether to show download progress
            provider: Federation provider id. Defaults to the client-level
                provider; an explicit value wins.

        Returns:
            File path if saved to disk, or AnnData object if loaded in memory
        """
        try:
            response = self._request_sample_download(sample_uuid, provider=provider)
            download_url = response["download_url"]
            filename = response.get("filename", f"{sample_uuid}.h5ad")
            file_size = response.get("file_size", 0)
            expires_in = response.get("expires_in_seconds", 0)

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
                with open(output_path, "wb") as f:
                    for chunk in file_response.iter_content(chunk_size=chunk_size):
                        if chunk:
                            f.write(chunk)
                            downloaded_bytes += len(chunk)

                            if show_progress and file_size > 0:
                                progress = (downloaded_bytes / file_size) * 100
                                print(
                                    f"\rDownloading: {progress:.1f}% ({downloaded_bytes:,}/{file_size:,} bytes)",
                                    end="",
                                )

                if show_progress:
                    print()  # New line after progress

                logger.info(f"Downloaded sample {sample_uuid} to {output_path}")
                return str(output_path)
            else:
                # Load into memory
                try:
                    import tempfile

                    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp_file:
                        downloaded_bytes = 0
                        for chunk in file_response.iter_content(chunk_size=chunk_size):
                            if chunk:
                                tmp_file.write(chunk)
                                downloaded_bytes += len(chunk)

                                if show_progress and file_size > 0:
                                    progress = (downloaded_bytes / file_size) * 100
                                    print(
                                        f"\rDownloading: {progress:.1f}% ({downloaded_bytes:,}/{file_size:,} bytes)",
                                        end="",
                                    )

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
                    raise ImportError(
                        "anndata is required to load H5AD files into memory. "
                        "Install with: pip install anndata"
                    )

        except requests.exceptions.RequestException as e:
            if "404" in str(e):
                raise MalvaAPIError(f"Sample {sample_uuid} not found")
            else:
                raise MalvaAPIError(f"Download failed: {e}")
        except Exception as e:
            raise MalvaAPIError(f"Download failed: {e}")

    def check_sample_availability(
        self, sample_uuid: str, provider: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Check if a sample is available for download without creating a download token

        Args:
            sample_uuid: UUID of the sample to check
            provider: Federation provider id. Defaults to the client-level
                provider; an explicit value wins.

        Returns:
            Dictionary with availability information
        """
        try:
            response = self._request(
                "GET",
                f"/api/samples/{sample_uuid}/download/info",
                params=self._provider_params(provider) or None,
            )
            return response
        except MalvaAPIError as e:
            if "404" in str(e) or "not found" in str(e).lower():
                return {
                    "sample_uuid": sample_uuid,
                    "downloadable": False,
                    "error": "Sample not found",
                }
            raise

    def get_sample_metadata(
        self, sample_ids: List[int], provider: Optional[str] = None
    ) -> Dict[int, Dict[str, Any]]:
        """
        Get metadata for specific samples

        Args:
            sample_ids: List of sample IDs
            provider: Federation provider id. Defaults to the client-level
                provider; an explicit value wins.

        Returns:
            Dictionary mapping sample IDs to their metadata
        """
        data = {"ids": sample_ids}
        response = self._request(
            "POST",
            self._endpoint("sample_metadata"),
            json=data,
            params=self._provider_params(provider) or None,
        )

        # Convert string keys back to integers if needed
        metadata = response.get("metadata", {})
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

    def get_cells_by_metadata(
        self,
        sample_ids: Optional[Union[int, List[int]]] = None,
        cell_types: Optional[Union[str, List[str]]] = None,
        include_cell_metadata: bool = True,
        include_sample_metadata: bool = False,
        include_all_database_cells: bool = False,
        poll_interval: int = 2,
        max_wait: int = 900,
        provider: Optional[str] = None,
    ) -> pd.DataFrame:
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

        provider_params = self._provider_params(provider)
        request_payload = {
            "sample_ids": sample_id_list or [],
            "cell_types": cell_type_list,
            "include_cell_metadata": include_cell_metadata,
            "include_all_database_cells": include_all_database_cells,
        }
        if include_all_database_cells and not sample_id_list and not cell_type_list:
            queued = self._request(
                "POST",
                self._endpoint("all_cells_export"),
                json=request_payload,
                params=provider_params or None,
            )
            export_id = queued.get("export_id")
            if not export_id:
                raise MalvaAPIError("No export_id received from server")
            start_time = time.time()
            status_payload = queued
            while time.time() - start_time < max_wait:
                status_payload = self._request(
                    "GET",
                    self._endpoint("all_cells_export_job", export_id=export_id),
                    params=provider_params or None,
                )
                status = status_payload.get("status")
                if status == "completed":
                    break
                if status == "error":
                    raise MalvaAPIError(
                        f"All-cells export failed: {status_payload.get('error', 'unknown')}"
                    )
                delay = self._adaptive_poll_delay(time.time() - start_time, poll_interval)
                if delay > 0:
                    time.sleep(delay)
            else:
                raise MalvaAPIError(f"All-cells export timed out after {max_wait} seconds")
            response = self._request_raw(
                "GET",
                self._endpoint("all_cells_export_download", export_id=export_id),
                params=provider_params or None,
            )
        else:
            response = self._request_raw(
                "POST",
                self._endpoint("cells_by_metadata"),
                json=request_payload,
                params=provider_params or None,
            )
        decoded = msgpack.unpackb(response.content, raw=False, strict_map_key=False)
        cells = np.frombuffer(decoded.get("c", b""), dtype=np.int64)
        samples = np.frombuffer(decoded.get("s", b""), dtype=np.int64)
        n = min(len(cells), len(samples))
        df = pd.DataFrame({"sample_id": samples[:n], "cell_id": cells[:n]})

        if include_cell_metadata and n:
            dense_valid = bool(decoded.get("_vi_dense", False))
            if decoded.get("_bin_ct_i16"):
                ct_ids = np.frombuffer(decoded.get("_bin_ct_i16", b""), dtype=np.int16)
            else:
                ct_ids = np.frombuffer(decoded.get("_bin_ct", b""), dtype=np.int32)
            if decoded.get("_bin_tc_u16"):
                total_counts_u16 = np.frombuffer(decoded.get("_bin_tc_u16", b""), dtype=np.uint16)
                total_counts = None
            else:
                total_counts = np.frombuffer(decoded.get("_bin_tc", b""), dtype=np.float32)
                total_counts_u16 = None
            valid_all = (
                None if dense_valid else np.frombuffer(decoded.get("_bin_vi", b""), dtype=np.int64)
            )
            mapping = {int(k): v for k, v in (decoded.get("cell_type_mapping") or {}).items()}

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
                codes = np.full(
                    n, -1, dtype=np.int16 if len(sorted_ids) <= np.iinfo(np.int16).max else np.int32
                )
                if sorted_ids and len(ids):
                    max_id = max(max(sorted_ids), int(ids.max(initial=0)))
                    if max_id <= 1_000_000:
                        remap = np.full(max_id + 1, -1, dtype=codes.dtype)
                        remap[np.asarray(sorted_ids, dtype=np.int64)] = np.arange(
                            len(sorted_ids), dtype=codes.dtype
                        )
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
                            codes[: len(mapped)] = mapped
                        else:
                            codes[valid_ct] = mapped
                df["cell_type"] = pd.Categorical.from_codes(codes, categories=categories)

            counts_u16 = np.zeros(n, dtype=np.uint16)
            if total_counts_u16 is not None:
                mt = min(
                    n,
                    len(total_counts_u16) if dense_valid else len(valid_all),
                    len(total_counts_u16),
                )
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
                        counts = np.nan_to_num(
                            total_counts[:mt], nan=0.0, posinf=65535.0, neginf=0.0
                        )
                        counts_u16[:mt] = np.rint(np.clip(counts, 0, 65535)).astype(np.uint16)
                    else:
                        valid_tc = valid_all[:mt]
                        in_bounds = (valid_tc >= 0) & (valid_tc < n)
                        valid_tc = valid_tc[in_bounds]
                        counts = np.nan_to_num(
                            total_counts[:mt][in_bounds], nan=0.0, posinf=65535.0, neginf=0.0
                        )
                        counts_u16[valid_tc] = np.rint(np.clip(counts, 0, 65535)).astype(np.uint16)
            df["total_counts"] = counts_u16

        if include_sample_metadata and n:
            meta = self._sample_metadata_frame(
                df["sample_id"].drop_duplicates().astype(int).tolist(),
                provider=provider,
            )
            if not meta.empty:
                df = df.merge(meta, on="sample_id", how="left")
        return df

    def get_database_stats(self) -> Dict[str, Any]:
        """Get database statistics"""
        return self._request("GET", "/api/stats")

    def get_available_filters(self, provider: Optional[str] = None) -> Dict[str, Any]:
        """Get available filter options for samples."""
        return self._request(
            "GET",
            self._endpoint("available_filters"),
            params=self._provider_params(provider) or None,
        )

    @staticmethod
    def _workspace_provider_ids(workspace: Dict[str, Any]) -> List[str]:
        """Collect provider ids from the shapes returned by ``/api/workspace``."""
        ids: List[str] = []

        def _add(value: Any) -> None:
            if value is None:
                return
            text = str(value).strip()
            if text and text not in ids:
                ids.append(text)

        if not isinstance(workspace, dict):
            return ids
        for value in workspace.get("provider_ids") or []:
            _add(value)
        for group in workspace.get("providers") or []:
            if isinstance(group, dict):
                _add(group.get("id") or group.get("provider_id"))
        for row in workspace.get("datasets") or []:
            if isinstance(row, dict):
                _add(row.get("provider_id"))
        return ids

    def _validate_provider_ids(self, provider_ids: List[str]) -> None:
        """Reject provider ids that the active workspace does not advertise.

        The Explorer otherwise creates a zero-shard federated job that stays
        ``running`` for an unknown provider, which is indistinguishable from a
        slow search. Validation is skipped when the workspace endpoint is not
        available or advertises no providers.
        """
        requested = [str(provider_id) for provider_id in provider_ids if str(provider_id)]
        if not requested:
            return
        try:
            workspace = self.get_workspace()
        except MalvaAPIError as exc:
            if not self._is_not_found_error(exc):
                logger.warning("Could not validate federation providers: %s", exc)
            return

        known = self._workspace_provider_ids(workspace)
        if not known:
            return
        unknown = [provider_id for provider_id in requested if provider_id not in known]
        if unknown:
            raise MalvaAPIError(
                "Unknown federation provider(s): "
                + ", ".join(unknown)
                + ". Available providers: "
                + ", ".join(known)
            )

    def search_federated(
        self,
        query: str,
        provider_ids: Optional[List[str]] = None,
        dataset_id: Optional[str] = None,
        poll_interval: int = 2,
        max_wait: int = 300,
        strict: bool = False,
    ) -> pd.DataFrame:
        """Search every configured Explorer provider and keep the origin.

        This calls the Explorer's federated endpoints, not the monorepo
        ``/search`` endpoint. The returned DataFrame has ``provider`` and
        ``provider_id`` as its first two columns, followed by the usual
        expression columns.

        Provenance is always attached to ``df.attrs``:

        - ``job_id``: federated job id
        - ``status``: final aggregate status reported by the Explorer
        - ``shards``: final per-provider shard dictionaries
        - ``providers``: concise per-provider metadata (id, name, status, error)
        - ``errors``: human-readable per-provider errors
        - ``shard_errors``: the same errors as structured dictionaries

        Failed shards are never silently dropped: they emit a warning and are
        listed in ``attrs["errors"]``. Pass ``strict=True`` to raise instead.

        Args:
            query: Gene symbol, description, or DNA sequence.
            provider_ids: Optional provider ids to search. When omitted, the
                default list is resolved from the provider ids advertised by
                ``get_workspace()``. If that lookup fails or the workspace
                advertises no providers, the field is omitted. Explicit ids
                always take precedence and are validated against
                ``get_workspace()`` before submitting.
            dataset_id: Search index for this call. Defaults to the client's
                configured ``dataset_id``.
            poll_interval: Seconds between federated status polls.
            max_wait: Maximum seconds to wait for all shards.
            strict: Raise :class:`MalvaAPIError` if any shard failed, even when
                other shards returned rows.

        Returns:
            A pandas DataFrame with one row per provider, sample and cell type.
            Empty successful searches still carry the declared schema and
            provenance.

        Raises:
            MalvaAPIError: If every shard failed, a requested provider is
                unknown, any shard failed while ``strict=True``, or the job
                timed out (including a zero-shard job that never finishes).
        """
        requested_ids = [str(provider_id) for provider_id in (provider_ids or [])]
        if requested_ids:
            self._validate_provider_ids(requested_ids)
        else:
            # The Explorer decides whether a POST /api/search is federated by
            # the presence of provider_ids. Resolve the workspace defaults so
            # an implicit federated search does not become a normal search.
            try:
                workspace = self.get_workspace()
            except MalvaAPIError as exc:
                logger.warning(
                    "Could not resolve default federation providers; omitting provider_ids: %s",
                    exc,
                )
            else:
                requested_ids = self._workspace_provider_ids(workspace)
                if not requested_ids:
                    logger.warning(
                        "Workspace advertised no federation providers; "
                        "omitting provider_ids from the federated search"
                    )

        payload: Dict[str, Any] = {"query": query}
        dataset = self._default_dataset_id(dataset_id)
        if dataset:
            payload["dataset_id"] = dataset
        if requested_ids:
            payload["provider_ids"] = requested_ids

        submitted = self._request("POST", self._endpoint("search_federated"), json=payload)
        job_id = submitted.get("job_id")
        if not job_id:
            raise MalvaAPIError("No federated job_id received from Explorer")

        shards = list(submitted.get("shards") or [])
        job_status = submitted.get("status")
        job_errors: List[str] = []
        if submitted.get("error"):
            job_errors.append(str(submitted["error"]))

        terminal = {"completed", "failed", "error", "canceled", "interrupted"}
        deadline = time.time() + max_wait
        while True:
            if shards and all(shard.get("status") in terminal for shard in shards):
                break
            if not shards and job_status in terminal:
                break
            if time.time() >= deadline:
                if not shards:
                    raise MalvaAPIError(
                        f"Federated search job {job_id} returned no shards and "
                        f"never left status {job_status!r} within {max_wait} "
                        "seconds. The requested provider may be unknown or its "
                        "search service may be unreachable."
                    )
                raise MalvaAPIError(f"Federated search timed out after {max_wait} seconds")

            elapsed = max(0.0, max_wait - (deadline - time.time()))
            delay = self._adaptive_poll_delay(elapsed, poll_interval)
            if delay > 0:
                time.sleep(delay)
            state = self._request("GET", self._endpoint("search_federated_job", job_id=job_id))
            if isinstance(state, dict):
                if state.get("status"):
                    job_status = state.get("status")
                if state.get("shards"):
                    shards = list(state["shards"])
                if state.get("error"):
                    job_errors.append(str(state["error"]))
                elif job_status in {"failed", "error"} and not shards:
                    job_errors.append(f"Federated job {job_id} failed with no shards")

        if not shards and job_status in {"failed", "error"} and not job_errors:
            job_errors.append(f"Federated job {job_id} finished with status {job_status}")

        frames: List[pd.DataFrame] = []
        errors: List[str] = []
        shard_errors: List[Dict[str, Any]] = []
        provider_metadata: List[Dict[str, Any]] = []

        def _record_error(shard: Dict[str, Any], provider_name: str, message: str) -> None:
            errors.append(f"{provider_name}: {message}")
            shard_errors.append(
                {
                    "provider_id": shard.get("provider_id"),
                    "provider_name": provider_name,
                    "status": shard.get("status"),
                    "error": message,
                }
            )

        for shard in shards:
            if not isinstance(shard, dict):
                continue
            provider_id = shard.get("provider_id")
            provider_name = shard.get("provider_name") or provider_id or "unknown"
            provider_metadata.append(
                {
                    "provider_id": provider_id,
                    "provider_name": provider_name,
                    "status": shard.get("status"),
                    "error": shard.get("error"),
                    "search_job_id": shard.get("search_job_id"),
                }
            )
            if shard.get("status") != "completed" or not shard.get("search_job_id"):
                message = shard.get("error") or (
                    f"shard {shard.get('status') or 'unknown'}"
                    if shard.get("status") != "completed"
                    else "completed without a search_job_id"
                )
                _record_error(shard, provider_name, str(message))
                continue

            params = {"provider": provider_id} if provider_id else None
            try:
                data = self._request(
                    "GET",
                    self._endpoint("search_shard_results", job_id=shard["search_job_id"]),
                    params=params,
                )
            except MalvaAPIError as exc:
                _record_error(shard, provider_name, str(exc))
                continue

            genes = {}
            if isinstance(data, dict):
                genes = data.get("genes") or data.get("results") or {}
            if not genes:
                continue
            result = SearchResult(
                {
                    "results": genes,
                    "job_id": shard["search_job_id"],
                    "aggregate_expression": True,
                },
                self,
            )
            frame = result.df.copy()
            if frame.empty:
                continue
            frame.insert(0, "provider_id", provider_id)
            frame.insert(0, "provider", provider_name)
            frames.append(frame)

        if job_errors:
            for message in job_errors:
                if message not in errors:
                    errors.append(message)
                    shard_errors.append(
                        {
                            "provider_id": None,
                            "provider_name": "federated-job",
                            "status": job_status,
                            "error": message,
                        }
                    )

        if errors:
            message = "Federated search had failed shard(s): " + "; ".join(errors)
            logger.warning(message)
            warnings.warn(message, RuntimeWarning, stacklevel=2)

        if strict and errors:
            raise MalvaAPIError(
                message,
                response_data={"job_id": job_id, "shards": shards, "errors": errors},
            )

        attrs = {
            "job_id": job_id,
            "status": job_status,
            "shards": shards,
            "providers": provider_metadata,
            "provider_metadata": provider_metadata,
            "errors": errors,
            "shard_errors": shard_errors,
        }

        if not frames:
            if errors:
                raise MalvaAPIError(
                    "Federated search failed: " + "; ".join(errors),
                    response_data=attrs,
                )
            empty = pd.DataFrame(columns=_FEDERATED_RESULT_COLUMNS)
            empty.attrs.update(attrs)
            return empty

        combined = pd.concat(frames, ignore_index=True)
        combined.attrs.update(attrs)
        return combined

    def submit_search(
        self,
        query: str,
        aggregate_expression: bool = True,
        window_size: Optional[int] = None,
        threshold: Optional[float] = None,
        stranded: Optional[bool] = None,
    ) -> str:
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

        from malva_client.storage import save_job_status, save_search

        data = {
            "query": query,
            "wait_for_completion": False,
            "aggregate_expression": aggregate_expression,
        }
        if window_size is not None:
            data["window_size"] = window_size
        if threshold is not None:
            data["threshold"] = threshold
        if stranded is not None:
            data["stranded"] = stranded
        if self.dataset_id:
            data["dataset_id"] = self.dataset_id

        response = self._request("POST", self._endpoint("search"), json=data)
        job_id = response["job_id"]

        # Save to local storage
        search_id = save_search(
            query=query, job_id=job_id, server_url=self.base_url, status="submitted"
        )

        save_job_status(job_id=job_id, status="pending", query=query, server_url=self.base_url)

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
                status=response.get("status", "unknown"),
                server_url=self.base_url,
                results_data=None,
                error_message=response.get("error") if response.get("status") == "error" else None,
            )

            return response

        except MalvaAPIError as e:
            if "404" in str(e) or "not found" in str(e).lower():
                # Job not found on server, check local storage
                from malva_client.storage import get_job_status as get_local_status

                local_status = get_local_status(job_id)
                if local_status:
                    return {
                        "job_id": job_id,
                        "status": local_status["status"],
                        "error": "Job not found on server (may have expired)",
                    }
                else:
                    return {"job_id": job_id, "status": "not_found", "error": "Job not found"}
            raise

    def get_job_results(self, job_id: str) -> "SearchResult":
        """
        Get results of a completed job

        Args:
            job_id: Job ID

        Returns:
            SearchResult object
        """
        from malva_client.storage import get_job_status as get_local_status
        from malva_client.storage import save_job_status

        # Check local cache — always use lazy loading via job_id
        local_status = get_local_status(job_id)
        if local_status and local_status.get("status") == "completed":
            logger.info(f"Job {job_id} found in local cache, using lazy load")
            return SearchResult({"job_id": job_id, "status": "completed"}, self)

        # Otherwise fetch from server
        status_response = self.get_job_status(job_id)

        if status_response["status"] == "completed":
            # Use lazy loading — SearchResult fetches the job results endpoint
            save_job_status(
                job_id=job_id, status="completed", server_url=self.base_url, results_data={}
            )
            return SearchResult({"job_id": job_id, "status": "completed"}, self)

        elif status_response["status"] == "error":
            error_msg = status_response.get("error", "Unknown error")
            raise SearchError(f"Job {job_id} failed: {error_msg}")

        elif status_response["status"] in ["pending", "running"]:
            raise SearchError(
                f"Job {job_id} is still {status_response['status']}. Use get_job_status() to check progress."
            )

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
            if entry.get("job_id"):
                jobs.append(
                    {
                        "job_id": entry["job_id"],
                        "query": entry["query"],
                        "status": entry["status"],
                        "timestamp": entry["timestamp"],
                        "results_count": entry.get("results_count", 0),
                        "error": entry.get("error_message"),
                    }
                )

        return jobs

    def wait_for_job(
        self, job_id: str, poll_interval: int = 2, max_wait: int = 300
    ) -> "SearchResult":
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
            status = status_response.get("status")

            if status == "completed":
                return self.get_job_results(job_id)
            elif status == "error":
                error_msg = status_response.get("error", "Unknown error")
                raise SearchError(f"Job {job_id} failed: {error_msg}")
            elif status in ["pending", "running"]:
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
            response = self._request("DELETE", self._endpoint("search_job", job_id=job_id))

            # Update local storage
            from malva_client.storage import save_job_status

            save_job_status(job_id=job_id, status="cancelled", server_url=self.base_url)

            return response.get("success", True)

        except MalvaAPIError as e:
            if "404" in str(e):
                logger.warning(f"Job {job_id} not found for cancellation")
                return False
            raise

    # ── Coverage Analysis Methods ──

    def get_coverage(
        self,
        chr: str,
        start: int,
        end: int,
        strand: str = "+",
        zoom: int = 1,
        metadata_filters: Optional[Dict[str, Any]] = None,
        poll_interval: int = 2,
        max_wait: int = 300,
        stage_callback: Optional[Callable[[str, Dict[str, Any]], None]] = None,
    ) -> "CoverageResult":
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
            stage_callback: Optional callable receiving coverage stage updates
                as ``callback(stage, payload)``. Stages include search polling,
                result generation, response download, and local JSON parsing.

        Returns:
            CoverageResult object
        """
        # Map +/- notation to server's expected format
        strand_map = {
            "+": "forward",
            "-": "reverse",
            "both": "both",
            "forward": "forward",
            "reverse": "reverse",
        }
        data = {
            "chromosome": chr,
            "start": start,
            "end": end,
            "strand": strand_map.get(strand, "both"),
            "zoom_level": zoom,
        }
        if metadata_filters:
            data["metadata_filters"] = metadata_filters

        self._emit_stage(
            stage_callback, "coverage_search_submitting", region=f"{chr}:{start}-{end}"
        )
        response = self._request("POST", "/api/genome-browser/search", json=data)
        job_id = response.get("job_id")

        if not job_id:
            raise MalvaAPIError("No job_id received from coverage search")
        self._emit_stage(
            stage_callback,
            "coverage_search_submitted",
            job_id=job_id,
            probes=len(response.get("positions") or []),
        )

        region = f"{chr}:{start}-{end}({strand})"
        return self._poll_coverage(
            job_id,
            region,
            positions=response.get("positions"),
            poll_interval=poll_interval,
            max_wait=max_wait,
            stage_callback=stage_callback,
        )

    def get_sequence_coverage(
        self,
        sequence: str,
        sequence_name: str = "",
        metadata_filters: Optional[Dict[str, Any]] = None,
        poll_interval: int = 2,
        max_wait: int = 300,
        stage_callback: Optional[Callable[[str, Dict[str, Any]], None]] = None,
    ) -> "CoverageResult":
        """
        Get coverage for an arbitrary DNA sequence.

        Args:
            sequence: DNA sequence
            sequence_name: Optional name for the sequence
            metadata_filters: Optional metadata filters
            poll_interval: How often to check for completion (seconds)
            max_wait: Maximum time to wait for completion (seconds)
            stage_callback: Optional callable receiving coverage stage updates

        Returns:
            CoverageResult object
        """
        data = {
            "sequence": sequence,
        }
        if sequence_name:
            data["sequence_name"] = sequence_name
        if metadata_filters:
            data["metadata_filters"] = metadata_filters

        self._emit_stage(
            stage_callback,
            "coverage_search_submitting",
            region=sequence_name or f"seq:{len(sequence)}bp",
        )
        response = self._request("POST", "/api/genome-browser/search-sequence", json=data)
        job_id = response.get("job_id")

        if not job_id:
            raise MalvaAPIError("No job_id received from sequence coverage search")
        self._emit_stage(
            stage_callback,
            "coverage_search_submitted",
            job_id=job_id,
            probes=len(response.get("positions") or []),
        )

        region = sequence_name or f"seq:{sequence[:20]}..."
        return self._poll_coverage(
            job_id,
            region,
            positions=response.get("positions"),
            poll_interval=poll_interval,
            max_wait=max_wait,
            stage_callback=stage_callback,
        )

    def _poll_coverage(
        self,
        job_id: str,
        region: str = "",
        positions: Optional[List[int]] = None,
        poll_interval: int = 2,
        max_wait: int = 300,
        stage_callback: Optional[Callable[[str, Dict[str, Any]], None]] = None,
    ) -> "CoverageResult":
        """
        Internal helper to poll for coverage results.

        Newer explorer-backed APIs expose the same /api/coverage/status and
        /api/coverage/results endpoints used by the web Coverage Explorer. Use
        those when submit returned probe positions, and fall back to the older
        /api/genome-browser/coverage endpoint for compatibility.
        """
        start_time = time.time()

        def _fetch_binary_compact_coverage() -> Optional["CoverageResult"]:
            try:
                result = self._request_msgpack_with_transfer_timing(
                    "GET",
                    f"/api/genome-browser/coverage-compact/{job_id}",
                    stage_callback=stage_callback,
                    stage_prefix="coverage_result",
                )
            except MalvaAPIError as exc:
                if self._is_not_found_error(exc):
                    return None
                raise
            coverage_data = result.get("coverage", result.get("coverage_data", result))
            coverage_data["job_id"] = job_id
            coverage_data["region"] = region
            return CoverageResult(coverage_data, self)

        if positions is not None:
            try:
                while time.time() - start_time < max_wait:
                    response = self._request("GET", f"/api/coverage/status/{job_id}")
                    status = response.get("status")
                    if status == "completed":
                        self._emit_stage(stage_callback, "coverage_search_completed", job_id=job_id)
                        compact = _fetch_binary_compact_coverage()
                        if compact is not None:
                            return compact
                        result = self._request_json_with_transfer_timing(
                            "POST",
                            f"/api/coverage/results/{job_id}",
                            json={"positions": positions, "filters": {}},
                            stage_callback=stage_callback,
                            stage_prefix="coverage_result",
                        )
                        coverage_data = result.get("coverage", result.get("coverage_data", result))
                        coverage_data["job_id"] = job_id
                        coverage_data["region"] = region
                        if "probe_data" in result:
                            coverage_data["probe_data"] = result["probe_data"]
                        return CoverageResult(coverage_data, self)
                    if status == "error":
                        error_msg = response.get("error", "Unknown error")
                        raise MalvaAPIError(f"Coverage job failed: {error_msg}")
                    self._emit_stage(
                        stage_callback, "coverage_search_polling", job_id=job_id, status=status
                    )
                    time.sleep(poll_interval)
                raise MalvaAPIError(f"Coverage job {job_id} timed out after {max_wait} seconds")
            except MalvaAPIError as exc:
                if not self._is_not_found_error(exc):
                    raise
                self._emit_stage(stage_callback, "coverage_result_fallback", job_id=job_id)

        while time.time() - start_time < max_wait:
            compact = _fetch_binary_compact_coverage()
            if compact is not None:
                return compact

            response = self._request("GET", f"/api/genome-browser/coverage/{job_id}")
            status = response.get("status")

            if status == "completed":
                coverage_data = response.get("coverage_data", response)
                coverage_data["job_id"] = job_id
                coverage_data["region"] = region
                if "probe_data" in response:
                    coverage_data["probe_data"] = response["probe_data"]
                return CoverageResult(coverage_data, self)
            elif status == "error":
                error_msg = response.get("error", "Unknown error")
                raise MalvaAPIError(f"Coverage job failed: {error_msg}")

            logger.info(f"Coverage job {job_id} status: {status}")
            time.sleep(poll_interval)

        raise MalvaAPIError(f"Coverage job {job_id} timed out after {max_wait} seconds")

    def get_coverage_data(self, job_id: str, **filters) -> "CoverageResult":
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

        response = self._request("GET", f"/api/genome-browser/coverage/{job_id}", params=params)
        coverage_data = response.get("coverage_data", response)
        coverage_data["job_id"] = job_id
        return CoverageResult(coverage_data, self)

    def download_coverage_wig(
        self, job_id: str, output_path: str, poll_interval: int = 2, max_wait: int = 120, **filters
    ) -> str:
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
        url = urljoin(self.base_url, f"/api/genome-browser/coverage-file/{job_id}")

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
        return self._request("GET", f"/api/genome-browser/filter-options/{job_id}")

    # ── Dataset Discovery Methods ──

    def get_datasets(self, provider: Optional[str] = None) -> List[Dict[str, Any]]:
        """
        List the datasets/indexes this account can search.

        Pass one of the returned ``dataset_id`` values to
        :meth:`set_workspace` (or the constructor) to scope subsequent
        searches to that index.
        """
        response = self._request(
            "GET", self._endpoint("datasets"), params=self._provider_params(provider) or None
        )
        if isinstance(response, list):
            return response
        return response.get("datasets", []) if isinstance(response, dict) else []

    def get_datasets_hierarchy(self, provider: Optional[str] = None) -> Dict[str, Any]:
        """
        Get the full dataset → study → sample hierarchy.

        Args:
            provider: Federation provider id. Defaults to the client-level
                provider; an explicit value wins.

        Returns:
            Dictionary with the full hierarchy tree
        """
        return self._request(
            "GET",
            self._endpoint("datasets_hierarchy"),
            params=self._provider_params(provider) or None,
        )

    def get_dataset_studies(
        self,
        dataset_id: str,
        page: int = 1,
        page_size: int = 50,
        provider: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Get studies belonging to a dataset.

        If the dataset-studies endpoint returns no studies for the first page,
        the studies are derived from the dataset hierarchy (the public provider
        currently needs this fallback).

        Args:
            dataset_id: Dataset identifier
            page: Page number (1-based)
            page_size: Number of studies per page
            provider: Federation provider id. Defaults to the client-level
                provider; an explicit value wins.

        Returns:
            Dictionary with studies and pagination info
        """
        params = {"page": page, "page_size": page_size}
        params.update(self._provider_params(provider))
        response = self._request(
            "GET", self._endpoint("dataset_studies", dataset_id=dataset_id), params=params
        )
        if response.get("studies") or page != 1:
            return response

        studies = self._studies_from_hierarchy(provider=provider)
        studies = [study for study in studies if str(study.get("dataset_id")) == str(dataset_id)]
        if not studies:
            return response

        logger.warning(
            "Dataset studies endpoint returned no studies for %s; using the "
            "dataset hierarchy fallback (%d studies).",
            dataset_id,
            len(studies),
        )
        start = (page - 1) * page_size
        return {
            "dataset_id": dataset_id,
            "studies": studies[start : start + page_size],
            "total_count": len(studies),
            "page": page,
            "page_size": page_size,
            "has_more": start + page_size < len(studies),
            "query_time_ms": 0,
            "source": "dataset_hierarchy",
        }

    def get_study_samples(
        self,
        dataset_id: str,
        study_name: str,
        page: int = 1,
        page_size: int = 50,
        provider: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Get samples belonging to a study within a dataset.

        Args:
            dataset_id: Dataset identifier
            study_name: Study name
            page: Page number (1-based)
            page_size: Number of samples per page
            provider: Federation provider id. Defaults to the client-level
                provider; an explicit value wins.

        Returns:
            Dictionary with samples and pagination info

        Raises:
            MalvaAPIError: If the server returns no samples for a study that the
                dataset hierarchy says has samples. This currently happens for
                some public-provider routes and requires a server-side fix; the
                client never fabricates sample rows.
        """
        params = {"page": page, "page_size": page_size}
        params.update(self._provider_params(provider))
        response = self._request(
            "GET",
            self._endpoint("study_samples", dataset_id=dataset_id, study_name=study_name),
            params=params,
        )
        if response.get("samples") or page != 1:
            return response

        hierarchy_studies = self._studies_from_hierarchy(provider=provider)
        hierarchy_study = next(
            (
                study
                for study in hierarchy_studies
                if str(study.get("dataset_id")) == str(dataset_id)
                and str(study.get("study_name") or study.get("name")) == str(study_name)
            ),
            None,
        )
        if hierarchy_study and int(hierarchy_study.get("sample_count") or 0) > 0:
            provider_id = params.get("provider") or "the active provider"
            raise MalvaAPIError(
                f"Server returned no samples for study {study_name!r} in dataset "
                f"{dataset_id!r} ({provider_id}), although the dataset hierarchy "
                f"reports {hierarchy_study['sample_count']} sample(s). The "
                "provider-scoped study-samples route is not returning hierarchy "
                "studies; this requires a server-side fix."
            )
        return response

    def get_sample_details(
        self, sample_uuid: str, provider: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Get full metadata for a single sample.

        Args:
            sample_uuid: UUID of the sample
            provider: Federation provider id. Defaults to the client-level
                provider; an explicit value wins.

        Returns:
            Dictionary with complete sample metadata

        Raises:
            MalvaAPIError: With the server's error body when metadata cannot be
                retrieved (the public provider currently returns HTTP 500 for
                some samples; that needs a server-side fix).
        """
        params = self._provider_params(provider)
        try:
            return self._request(
                "GET",
                self._endpoint("sample_details", sample_uuid=sample_uuid),
                params=params or None,
            )
        except MalvaAPIError as exc:
            if self._is_not_found_error(exc):
                raise
            provider_note = f" from provider {params['provider']!r}" if params else ""
            raise MalvaAPIError(
                f"Failed to retrieve metadata for sample {sample_uuid!r}{provider_note}: {exc}",
                getattr(exc, "response_data", None),
            ) from exc

    def get_filter_values(
        self,
        column: str,
        search: str = "",
        limit: int = 50,
        provider: Optional[str] = None,
    ) -> List[str]:
        """
        Get unique values for a filter column, with optional search.

        Args:
            column: Filter column name (e.g., 'organ', 'species')
            search: Optional text to filter values by
            limit: Maximum number of values to return
            provider: Federation provider id. Defaults to the client-level
                provider; an explicit value wins.

        Returns:
            List of unique values for the column
        """
        params = {"limit": limit}
        params.update(self._provider_params(provider))
        if search:
            params["search"] = search
        response = self._request(
            "GET", self._endpoint("filter_values", column=column), params=params
        )
        return response.get("values", [])

    def get_overview_stats(self, provider: Optional[str] = None) -> Dict[str, Any]:
        """
        Get platform-wide statistics.

        Args:
            provider: Federation provider id. Defaults to the client-level
                provider; an explicit value wins.

        Returns:
            Dictionary with overview statistics (total samples, studies, etc.)
        """
        return self._request(
            "GET", self._endpoint("overview"), params=self._provider_params(provider) or None
        )

    # ── Coexpression Analysis Methods ──

    def get_umap_coordinates(self, dataset_id: str) -> "UMAPCoordinates":
        """
        Get UMAP coordinates for a dataset's metacells.

        Args:
            dataset_id: Dataset identifier

        Returns:
            UMAPCoordinates object with x, y, cluster, and sample data
        """
        response = self._request("GET", f"/api/coexpression/umap/{dataset_id}")
        return UMAPCoordinates(response, self)

    def _coexpression_from_response(
        self,
        response: Dict[str, Any],
        task_endpoint: str,
        poll_interval: int = 2,
        max_wait: int = 300,
    ) -> "CoexpressionResult":
        """Poll async coexpression task responses and return a result object."""
        status = response.get("status")
        task_id = response.get("task_id")
        if status not in {"pending", "running"} or not task_id:
            return CoexpressionResult(response, self)

        start_time = time.time()
        while time.time() - start_time < max_wait:
            task = self._request("GET", f"{task_endpoint.rstrip('/')}/{task_id}")
            task_status = task.get("status")
            if task_status == "completed":
                return CoexpressionResult(task, self)
            if task_status in {"failed", "error"}:
                raise MalvaAPIError(f"Coexpression task failed: {task.get('error', 'unknown')}")
            logger.info(f"Coexpression task {task_id} status: {task_status}")
            time.sleep(poll_interval)
        raise MalvaAPIError(f"Coexpression task {task_id} timed out after {max_wait} seconds")

    @staticmethod
    def _coexpression_task_endpoint(endpoint: str) -> str:
        if endpoint.startswith("/api/coexpression/"):
            return "/api/coexpression/tasks"
        return "/tasks"

    def _wait_coexpression_task(
        self,
        response: Dict[str, Any],
        task_endpoint: str,
        poll_interval: int = 2,
        max_wait: int = 600,
    ) -> Dict[str, Any]:
        """Poll a coexpression task response until completion and return the raw result."""
        status = response.get("status")
        task_id = response.get("task_id")
        if status not in {"pending", "running"} or not task_id:
            return response

        start_time = time.time()
        while time.time() - start_time < max_wait:
            task = self._request("GET", f"{task_endpoint.rstrip('/')}/{task_id}")
            task_status = task.get("status")
            if task_status == "completed":
                return task
            if task_status in {"failed", "error"}:
                raise MalvaAPIError(f"Coexpression task failed: {task.get('error', 'unknown')}")
            logger.info(f"Coexpression task {task_id} status: {task_status}")
            time.sleep(poll_interval)
        raise MalvaAPIError(f"Coexpression task {task_id} timed out after {max_wait} seconds")

    def project_cells_to_metacells(
        self,
        cell_ids: List[int],
        sample_ids: List[int],
        dataset_id: str,
        poll_interval: int = 2,
        max_wait: int = 300,
    ) -> List[str]:
        """Map explicit cells to atlas metacell IDs without a coexpression query.

        This is the efficient group-definition step: send the cell IDs and
        encoded sample IDs, receive only the positive metacell IDs, then pass
        those IDs to :meth:`differential_expression`.
        """
        if len(cell_ids) != len(sample_ids):
            raise ValueError("cell_ids and sample_ids must have the same length")
        if not cell_ids:
            raise ValueError("At least one cell_id is required")

        data = {
            "dataset_id": dataset_id,
            "cell_ids": [int(c) for c in cell_ids],
            "sample_ids": [int(s) for s in sample_ids],
        }
        response, endpoint = self._request_first_available(
            "POST",
            ["/api/coexpression/cell-to-metacell", "/query/cell-to-metacell"],
            json=data,
        )
        task_endpoint = self._coexpression_task_endpoint(endpoint)
        result = self._wait_coexpression_task(
            response, task_endpoint, poll_interval=poll_interval, max_wait=max_wait
        )
        return list(result.get("metacell_ids") or [])

    def differential_expression(
        self,
        metacell_ids_a: List[str],
        metacell_ids_b: List[str],
        dataset_id: str,
        top_n_genes: int = 0,
        min_log2fc: float = 0.0,
        min_pct: float = 0.1,
        poll_interval: int = 2,
        max_wait: int = 600,
    ) -> DifferentialExpressionResult:
        """Run differential expression between two metacell groups.

        ``top_n_genes=0`` returns the full tested gene set, matching the
        Explorer analysis workflow.
        """
        if not metacell_ids_a or not metacell_ids_b:
            raise ValueError("Both metacell_ids_a and metacell_ids_b are required")

        data = {
            "dataset_id": dataset_id,
            "metacell_ids_a": [str(m) for m in metacell_ids_a],
            "metacell_ids_b": [str(m) for m in metacell_ids_b],
            "top_n_genes": int(top_n_genes),
            "min_log2fc": float(min_log2fc),
            "min_pct": float(min_pct),
        }
        response, endpoint = self._request_first_available(
            "POST",
            ["/api/coexpression/differential", "/query/differential"],
            json=data,
        )
        task_endpoint = self._coexpression_task_endpoint(endpoint)
        raw = self._wait_coexpression_task(
            response, task_endpoint, poll_interval=poll_interval, max_wait=max_wait
        )
        return DifferentialExpressionResult(raw, self)

    def differential_expression_from_cells(
        self,
        cell_ids_a: List[int],
        sample_ids_a: List[int],
        cell_ids_b: List[int],
        sample_ids_b: List[int],
        dataset_id: str,
        top_n_genes: int = 0,
        min_log2fc: float = 0.0,
        min_pct: float = 0.1,
        poll_interval: int = 2,
        max_wait: int = 600,
    ) -> DifferentialExpressionResult:
        """Resolve two explicit cell groups to metacells and run DE server-side.

        Use this when the caller already has raw cell IDs and does not need the
        intermediate metacell list; no sample metadata is sent.
        """
        for name, cells, samples in (
            ("A", cell_ids_a, sample_ids_a),
            ("B", cell_ids_b, sample_ids_b),
        ):
            if len(cells) != len(samples):
                raise ValueError(f"cell_ids_{name} and sample_ids_{name} must match")
            if not cells:
                raise ValueError(f"At least one cell_id is required for group {name}")

        data = {
            "dataset_id": dataset_id,
            "cell_ids_a": [int(c) for c in cell_ids_a],
            "sample_ids_a": [int(s) for s in sample_ids_a],
            "cell_ids_b": [int(c) for c in cell_ids_b],
            "sample_ids_b": [int(s) for s in sample_ids_b],
            "top_n_genes": int(top_n_genes),
            "min_log2fc": float(min_log2fc),
            "min_pct": float(min_pct),
        }
        response, endpoint = self._request_first_available(
            "POST",
            [
                "/api/coexpression/differential-from-cells",
                "/query/differential-from-cells",
            ],
            json=data,
        )
        task_endpoint = self._coexpression_task_endpoint(endpoint)
        raw = self._wait_coexpression_task(
            response, task_endpoint, poll_interval=poll_interval, max_wait=max_wait
        )
        return DifferentialExpressionResult(raw, self)

    def project_cells(
        self,
        cell_ids: List[int],
        sample_ids: List[int],
        dataset_id: str,
        poll_interval: int = 2,
        max_wait: int = 300,
        **kwargs,
    ) -> "CoexpressionResult":
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
            "dataset_id": dataset_id,
            "cell_ids": [int(c) for c in cell_ids],
            "sample_ids": [int(s) for s in sample_ids],
        }
        data.update(kwargs)

        response, endpoint = self._request_first_available(
            "POST",
            ["/api/coexpression/query", "/query"],
            json=data,
        )
        return self._coexpression_from_response(
            response,
            task_endpoint=self._coexpression_task_endpoint(endpoint),
            poll_interval=poll_interval,
            max_wait=max_wait,
        )

    def get_coexpression(
        self,
        job_id: str,
        dataset_id: str,
        filter_sample_ids: Optional[List[int]] = None,
        poll_interval: int = 2,
        max_wait: int = 300,
        **kwargs,
    ) -> "CoexpressionResult":
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
            "job_id": job_id,
            "dataset_id": dataset_id,
        }
        if filter_sample_ids:
            data["filter_sample_ids"] = [int(s) for s in filter_sample_ids]
        data.update(kwargs)
        response, endpoint = self._request_first_available(
            "POST",
            ["/api/coexpression/query-by-job", "/query-by-job"],
            json=data,
        )
        return self._coexpression_from_response(
            response,
            task_endpoint=self._coexpression_task_endpoint(endpoint),
            poll_interval=poll_interval,
            max_wait=max_wait,
        )

    def get_coexpression_genes(
        self,
        job_id: str,
        dataset_id: str,
        filter_sample_ids: Optional[List[int]] = None,
        poll_interval: int = 2,
        max_wait: int = 300,
        **kwargs,
    ) -> "CoexpressionResult":
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
            "job_id": job_id,
            "dataset_id": dataset_id,
            "include_go_enrichment": False,
            "include_umap_scores": False,
            "include_cell_type_enrichment": False,
            "include_tissue_breakdown": False,
        }
        if filter_sample_ids:
            data["filter_sample_ids"] = [int(s) for s in filter_sample_ids]
        data.update(kwargs)
        response, endpoint = self._request_first_available(
            "POST",
            ["/api/coexpression/query-by-job", "/query-by-job"],
            json=data,
        )
        return self._coexpression_from_response(
            response,
            task_endpoint=self._coexpression_task_endpoint(endpoint),
            poll_interval=poll_interval,
            max_wait=max_wait,
        )

    def find_metacells(
        self,
        adata,
        groupby: str,
        top_n: int = 50,
        dataset_id: Optional[str] = None,
        *,
        layer: Optional[str] = None,
        feature_ids: Optional[Iterable[str]] = None,
        space: str = "hellinger",
        sample_metadata: Optional[Union[pd.DataFrame, Dict[Any, Any]]] = None,
    ) -> MetacellMatchResult:
        """Find Malva atlas metacells that look most like each local cluster.

        The AnnData object stays local: clusters are averaged into pseudobulk
        profiles and only those few cluster × gene vectors are sent to Malva.

        Args:
            adata: AnnData object with cells in ``obs`` and genes in ``var``.
            groupby: Categorical ``adata.obs`` column defining the clusters.
            top_n: Number of metacell neighbours per cluster (1–1000, default 50).
            dataset_id: Malva dataset / coexpression index identifier.
            layer: Optional ``adata.layers`` key to use instead of ``adata.X``.
            feature_ids: Optional Malva feature IDs aligned to ``adata.var_names``.
                Defaults to ``adata.var_names``.
            space: Server representation used for similarity. ``"hellinger"``
                (default) scores the square root of the panel composition
                (Bhattacharyya cosine), which downweights housekeeping genes and
                generally improves rare-subtype matching. ``"clr"`` keeps the
                centred-log-ratio behaviour of older servers. Servers that do
                not know the field ignore it and return CLR.
            sample_metadata: Optional local sample metadata (DataFrame indexed by
                ``sample_id``, or a mapping ``sample_id -> row``). It is merged
                locally after the response; no metadata is requested from Malva.

        Returns:
            :class:`~malva_client.models.MetacellMatchResult`, a DataFrame-like
            object with ``to_dataframe()``, ``enrich()`` and ``assign_labels()``.
            Its columns are ``cluster``, ``rank``, ``metacell_id``, ``sample_id``,
            ``cell_type``, ``n_cells``, ``similarity``, ``distance``,
            ``enrichment_score`` and ``percentile`` (plus any local sample
            metadata columns). Similarity is computed server-side in CLR
            (centred-log-ratio) space on a fixed variable-gene panel; values are
            normally normalized expression (e.g. ``log1p(CP10K)``), and raw
            counts are also accepted.
        """
        if not dataset_id:
            raise ValueError("dataset_id is required")
        if not 1 <= int(top_n) <= 1000:
            raise ValueError("top_n must be between 1 and 1000")
        space = str(space or "hellinger").strip().lower()
        if space not in {"clr", "hellinger"}:
            raise ValueError("space must be 'clr' or 'hellinger'")
        if groupby not in adata.obs.columns:
            raise KeyError(f"'{groupby}' is not a column in adata.obs")

        if layer is None:
            X = adata.X
        else:
            if layer not in adata.layers:
                raise KeyError(f"'{layer}' is not a layer in adata.layers")
            X = adata.layers[layer]

        obs = adata.obs[groupby]
        if not isinstance(obs.dtype, pd.CategoricalDtype):
            warnings.warn(
                f"adata.obs['{groupby}'] is not categorical; treating each unique "
                "value as a cluster. For Scanpy-style workflows, make it categorical "
                "with adata.obs['{}'] = adata.obs['{}'].astype('category').".format(
                    groupby, groupby
                ),
                UserWarning,
                stacklevel=2,
            )
        categorical = obs if isinstance(obs.dtype, pd.CategoricalDtype) else obs.astype("category")
        codes = categorical.cat.codes.to_numpy(dtype=np.int64)
        labels = [str(label) for label in categorical.cat.categories]

        raw_counts = _looks_like_raw_counts(X)
        if raw_counts:
            matrix_name = "adata.X" if layer is None else f"adata.layers['{layer}']"
            warnings.warn(
                f"{matrix_name} looks like raw counts (non-negative integers). "
                "Malva can score these as proportions, but normalized expression is "
                "recommended: run sc.pp.normalize_total(adata, target_sum=1e4) "
                "followed by sc.pp.log1p(adata). The data is not modified.",
                UserWarning,
                stacklevel=2,
            )

        var_names = np.asarray([str(v) for v in adata.var_names], dtype=object)
        if feature_ids is None:
            gene_ids = var_names
        else:
            gene_ids = np.asarray([str(f) for f in feature_ids], dtype=object)
            if len(gene_ids) != len(var_names):
                raise ValueError(
                    "feature_ids must have one entry per adata.var_names ("
                    f"{len(var_names)} expected, got {len(gene_ids)})"
                )

        means = _cluster_pseudobulk(X, codes, len(labels))
        transform = "counts" if raw_counts else "auto"
        clusters = []
        empty_clusters = []
        for i, label in enumerate(labels):
            row = means.getrow(i)
            if row.nnz == 0:
                empty_clusters.append(label)
                continue
            clusters.append(
                {
                    "name": label,
                    "gene_ids": [str(gene_ids[j]) for j in row.indices],
                    "values": [float(v) for v in row.data],
                    "transform": transform,
                }
            )

        if empty_clusters:
            warnings.warn(
                "Skipping clusters with no non-zero expression: "
                + ", ".join(empty_clusters[:10])
                + ("..." if len(empty_clusters) > 10 else ""),
                UserWarning,
                stacklevel=2,
            )
        if not clusters:
            raise ValueError("No clusters with non-zero expression to compare")

        payload = {
            "dataset_id": dataset_id,
            "clusters": clusters,
            "top_n": int(top_n),
            "metric": "cosine",
            "space": space,
        }
        response, _endpoint = self._request_first_available(
            "POST",
            ["/api/coexpression/clusters-to-metacells", "/query/clusters-to-metacells"],
            json=payload,
        )

        missing_count = 0
        missing_examples: List[str] = []
        for cluster in response.get("clusters", []):
            n_missing = int(cluster.get("n_missing_genes") or 0)
            if n_missing:
                missing_count += n_missing
                missing_examples.extend(cluster.get("missing_genes") or [])
        if missing_count:
            examples = ", ".join(missing_examples[:10])
            warnings.warn(
                f"{missing_count} gene(s) were not found in the Malva dataset and "
                f"were ignored{': ' + examples if examples else ''}",
                UserWarning,
                stacklevel=2,
            )

        columns = [
            "cluster",
            "rank",
            "metacell_id",
            "sample_id",
            "cell_type",
            "n_cells",
            "similarity",
            "distance",
            "enrichment_score",
            "percentile",
        ]
        rows = []
        for cluster in response.get("clusters", []):
            for neighbor in cluster.get("neighbors", []):
                rows.append(
                    {
                        "cluster": cluster.get("name"),
                        "rank": neighbor.get("rank"),
                        "metacell_id": neighbor.get("metacell_id"),
                        "sample_id": neighbor.get("sample_id"),
                        "cell_type": neighbor.get("cell_type"),
                        "n_cells": neighbor.get("n_cells"),
                        "similarity": neighbor.get("similarity"),
                        "distance": neighbor.get("distance"),
                        "enrichment_score": neighbor.get("enrichment_score"),
                        "percentile": neighbor.get("percentile"),
                    }
                )
        frame = pd.DataFrame(rows, columns=columns)

        result = MetacellMatchResult(
            frame,
            client=self,
            dataset_id=dataset_id,
            groupby=groupby,
            cluster_order=labels,
            sample_metadata=sample_metadata,
            representation=response.get("representation"),
        )
        if sample_metadata is not None:
            result.enrich()
        return result

    def compare_clusters_to_metacells(
        self,
        adata,
        groupby: str,
        top_n: int = 50,
        dataset_id: Optional[str] = None,
        *,
        layer: Optional[str] = None,
        feature_ids: Optional[Iterable[str]] = None,
        space: str = "hellinger",
        sample_metadata: Optional[Union[pd.DataFrame, Dict[Any, Any]]] = None,
    ) -> MetacellMatchResult:
        """Deprecated alias for :meth:`find_metacells`."""
        warnings.warn(
            "compare_clusters_to_metacells() is deprecated; use find_metacells() instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.find_metacells(
            adata,
            groupby,
            top_n=top_n,
            dataset_id=dataset_id,
            layer=layer,
            feature_ids=feature_ids,
            space=space,
            sample_metadata=sample_metadata,
        )

    # ── Cell Explorer ─────────────────────────────────────────────────────
    def list_cell_explorer_cohorts(self) -> List[Dict[str, Any]]:
        """List the current user's Cell Explorer cohorts."""
        response = self._request("GET", "/api/cell-explorer/cohorts")
        return response.get("cohorts", [])

    def create_cell_explorer_cohort(
        self,
        name: str,
        *,
        cell_types: Optional[List[str]] = None,
        organs: Optional[List[str]] = None,
        diseases: Optional[List[str]] = None,
        genes: Optional[List[str]] = None,
        dataset_id: Optional[str] = None,
        description: str = "",
        filter_definition: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """Build a per-user metacell cohort and its Zarr store.

        Filters are applied at the same (sample × cell type) metacell level the
        search engine uses.  The cohort is also registered as an analysis cell
        group, so it can be used by the existing analysis widgets.

        Example::

            cohort = client.create_cell_explorer_cohort(
                "IPF lung macrophages",
                cell_types=["alveolar macrophage"],
                organs=["lung"],
                diseases=["idiopathic pulmonary fibrosis"],
                genes=["SPP1", "TREM2"],
            )
        """
        filters: Dict[str, Any] = {}
        if filter_definition:
            filters.update(filter_definition.get("filters") or {})
        if cell_types:
            filters["cell_type"] = list(cell_types)
        if organs:
            filters["specimen_from_organism.organ"] = list(organs)
        if diseases:
            filters["specimen_from_organism.diseases"] = list(diseases)
        payload = {
            "name": name,
            "description": description,
            "dataset_id": dataset_id,
            "filter_definition": {"filters": filters},
            "genes": genes or [],
        }
        return self._request("POST", "/api/cell-explorer/cohorts", json=payload)

    def get_cell_explorer_cohort(self, cohort_id: str) -> Dict[str, Any]:
        """Return one cohort's metadata."""
        return self._request("GET", f"/api/cell-explorer/cohorts/{cohort_id}")

    def preview_cell_explorer_cohort(self, cohort_id: str) -> Dict[str, Any]:
        """Return the cohort metadata plus a bounded obs/var/X preview."""
        return self._request("GET", f"/api/cell-explorer/cohorts/{cohort_id}/preview")

    def delete_cell_explorer_cohort(self, cohort_id: str) -> bool:
        """Delete a cohort, its Zarr store, and its analysis group."""
        self._request("DELETE", f"/api/cell-explorer/cohorts/{cohort_id}")
        return True

    def download_cell_explorer_cohort(
        self,
        cohort_id: str,
        *,
        output_path: Optional[str] = None,
        extract: bool = True,
    ) -> str:
        """Download a cohort Zarr bundle.

        Args:
            cohort_id: Cohort identifier.
            output_path: Destination ``.tar.gz`` path.  Defaults to a temp file.
            extract: Extract the ``.zarr`` directory next to the archive.

        Returns:
            Path to the extracted Zarr directory when ``extract=True``, else the
            archive path.
        """
        response = self._request_raw("GET", f"/api/cell-explorer/cohorts/{cohort_id}/download")
        if output_path:
            archive = Path(output_path)
            if archive.is_dir():
                archive = archive / f"{cohort_id}.tar.gz"
        else:
            handle = tempfile.NamedTemporaryFile(delete=False, suffix=".tar.gz")
            handle.close()
            archive = Path(handle.name)
        archive.write_bytes(response.content)
        if not extract:
            return str(archive)
        with tarfile.open(archive, "r:gz") as bundle:
            bundle.extractall(archive.parent)
        return str(archive.parent / f"{cohort_id}.zarr")

    def load_cell_explorer_cohort(self, cohort_id: str, *, output_path: Optional[str] = None):
        """Download and open a cohort as a :mod:`zarr` group.

        The returned group follows ``malva.cell_explorer.v1``:
        ``obs/*``, ``var/gene``, ``X``, ``layers/raw``, optional ``obsm/X_umap``.
        """
        try:
            import zarr  # noqa: PLC0415 - optional analysis dependency
        except ImportError as exc:
            raise MalvaAPIError(
                "Reading a Cell Explorer cohort requires zarr; install malva-client[zarr]"
            ) from exc
        path = self.download_cell_explorer_cohort(cohort_id, output_path=output_path, extract=True)
        return zarr.open_group(path, mode="r")


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
