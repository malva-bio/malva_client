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
from malva_client.models import SearchResult, CoverageResult, SingleCellResult
from malva_client.config import Config

# Set up logging
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
        
        # Test connection and authentication
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
                    # If quota endpoint doesn't exist, that's ok for connection test
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
            
            # Log response details for debugging
            logger.debug(f"{method} {url} -> {response.status_code}")
            logger.debug(f"Response content type: {response.headers.get('content-type', 'unknown')}")
            
            # Handle different status codes before trying to parse JSON
            if response.status_code == 404:
                content_type = response.headers.get('content-type', '')
                if 'text/html' in content_type:
                    raise MalvaAPIError(f"Endpoint not found: {endpoint}. Check if the Malva API is running at {self.base_url}")
                
                # Try to parse JSON error if it's not HTML
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
            
            # Handle authentication errors
            elif response.status_code in [401, 403]:
                try:
                    error_data = response.json()
                    error_msg = error_data.get('error', error_data.get('detail', 'Authentication failed'))
                except json.JSONDecodeError:
                    error_msg = "Authentication failed"
                raise AuthenticationError(error_msg)
            
            # For other error status codes, try to get JSON error message
            elif not response.ok:
                try:
                    error_data = response.json()
                    error_msg = error_data.get('detail', error_data.get('error', f'HTTP {response.status_code}'))
                    raise MalvaAPIError(f"API error: {error_msg}")
                except json.JSONDecodeError:
                    raise MalvaAPIError(f"HTTP {response.status_code}: {response.text[:200]}")
            
            # Success case - try to parse JSON
            try:
                return response.json()
            except json.JSONDecodeError as e:
                # If we can't parse JSON on a successful response, that's a problem
                content_type = response.headers.get('content-type', 'unknown')
                raise MalvaAPIError(f"Expected JSON response but got {content_type}. Content: {response.text[:200]}")
            
        except requests.exceptions.Timeout:
            raise MalvaAPIError(f"Request timed out after {self.timeout} seconds")
        except requests.exceptions.ConnectionError:
            raise MalvaAPIError(f"Could not connect to {self.base_url}. Check if the server is running.")
        except requests.exceptions.RequestException as e:
            raise MalvaAPIError(f"Request failed: {e}")

    def search(self, query: str, wait_for_completion: bool = True, 
               poll_interval: int = 2, max_wait: int = 300, aggregate_expression: bool = True) -> SearchResult:
        """
        Perform a natural language or gene search
        
        Args:
            query: Natural language query or gene symbol (e.g., "BRCA1 expression in cancer")
            wait_for_completion: Whether to wait for results or return immediately
            poll_interval: How often to check for completion (seconds)
            max_wait: Maximum time to wait for completion (seconds)
            
        Returns:
            SearchResult object containing the results
        """
        # Submit search using the /search endpoint
        data = {
            'query': query,
            'wait_for_completion': wait_for_completion,
            'max_wait': max_wait,
            'poll_interval': poll_interval,
            'aggregate_expression': aggregate_expression
        }
        
        response = self._request('POST', '/search', json=data)
        
        # If wait_for_completion was True and search completed
        if response.get('status') == 'completed':
            logger.info(f"Search completed with job ID: {response['job_id']}")
            return SearchResult(response['results'], self)
        
        # If search is still pending or we didn't wait
        job_id = response['job_id']
        logger.info(f"Search submitted with job ID: {job_id}")
        
        if not wait_for_completion or response.get('status') == 'pending':
            return SearchResult({'job_id': job_id, 'status': 'pending'}, self)
        
        # If there was an error during submissionx
        if not response.get('success', True):
            error_msg = response.get('error', 'Unknown error')
            raise MalvaAPIError(f"Search failed: {error_msg}")
        
        # Fallback: poll manually if the server-side waiting didn't work
        start_time = time.time()
        while time.time() - start_time < max_wait:
            try:
                results_response = self._request('GET', f'/search/{job_id}')
                status = results_response.get('status')
                
                logger.info(f"Search status: {status}")
                
                if status == 'completed':
                    logger.info(f"Search completed with job ID: {response['job_id']}")
                    return SearchResult(response, self)
                elif status == 'error':
                    error_msg = results_response.get('error', 'Unknown error')
                    raise MalvaAPIError(f"Search failed: {error_msg}")
                
            except MalvaAPIError as e:
                if "404" in str(e) or "not found" in str(e).lower():
                    # Job not found yet, continue polling
                    pass
                else:
                    raise
            
            time.sleep(poll_interval)
        
        raise MalvaAPIError(f"Search timed out after {max_wait} seconds")

    def search_cells(self, query: str, wait_for_completion: bool = True,
                    poll_interval: int = 2, max_wait: int = 300) -> SingleCellResult:
        """
        Perform a natural language or gene search and returns individual cells
        
        Args:
            query: Natural language query or gene symbol (e.g., "BRCA1 expression in cancer")
            wait_for_completion: Whether to wait for results or return immediately
            poll_interval: How often to check for completion (seconds)
            max_wait: Maximum time to wait for completion (seconds)
            
        Returns:
            SingleCellResult object containing the individual cell results
        """
        # Submit search using the /search endpoint
        data = {
            'query': query,
            'wait_for_completion': wait_for_completion,
            'max_wait': max_wait,
            'poll_interval': poll_interval,
            'aggregate_expression': False  # This is the key difference
        }
        
        response = self._request('POST', '/search', json=data)
        
        # If wait_for_completion was True and search completed
        if response.get('status') == 'completed':
            logger.info(f"Search completed with job ID: {response['job_id']}")
            return SingleCellResult(response['results'], self)
        
        # If search is still pending or we didn't wait
        job_id = response['job_id']
        logger.info(f"Search submitted with job ID: {job_id}")
        
        if not wait_for_completion or response.get('status') == 'pending':
            return SingleCellResult({'job_id': job_id, 'status': 'pending'}, self)
        
        # If there was an error during submission
        if not response.get('success', True):
            error_msg = response.get('error', 'Unknown error')
            raise MalvaAPIError(f"Search failed: {error_msg}")
        
        # Fallback: poll manually if the server-side waiting didn't work
        start_time = time.time()
        while time.time() - start_time < max_wait:
            try:
                results_response = self._request('GET', f'/search/{job_id}')
                status = results_response.get('status')
                
                logger.info(f"Search status: {status}")
                
                if status == 'completed':
                    logger.info(f"Search completed with job ID: {job_id}")
                    return SingleCellResult(results_response, self)
                elif status == 'error':
                    error_msg = results_response.get('error', 'Unknown error')
                    raise MalvaAPIError(f"Search failed: {error_msg}")
                    
            except MalvaAPIError as e:
                if "404" in str(e) or "not found" in str(e).lower():
                    # Job not found yet, continue polling
                    pass
                else:
                    raise
            
            time.sleep(poll_interval)
        
        raise MalvaAPIError(f"Search timed out after {max_wait} seconds")
    
    def search_sequence(self, sequence: str, **kwargs) -> SearchResult:
        """
        Search for a specific DNA sequence
        
        Args:
            sequence: DNA sequence to search for
            **kwargs: Additional arguments passed to search()
            
        Returns:
            SearchResult object
        """
        # Validate sequence length
        if len(sequence) > 500000:  # 500kb limit
            raise ValueError("Sequence length cannot exceed 500kb (500,000 nucleotides)")
        
        # Validate sequence contains only valid nucleotides
        valid_nucleotides = set('ATGCUN')
        if not all(c.upper() in valid_nucleotides for c in sequence):
            raise ValueError("Sequence contains invalid nucleotides. Only A, T, G, C, U, N are allowed.")
        
        return self.search(sequence, **kwargs)
    
    def search_genes(self, genes: List[str], **kwargs) -> SearchResult:
        """
        Search for multiple genes
        
        Args:
            genes: List of gene symbols
            **kwargs: Additional arguments passed to search()
            
        Returns:
            SearchResult object
        """
        # Join genes into a query
        gene_query = f"expression of {', '.join(genes)}"
        return self.search(gene_query, **kwargs)
    
    def get_coverage(self, chromosome: str, start: int, end: int, 
                    strand: str = 'both', split_by_cell_type: bool = False,
                    preserve_samples: bool = True,
                    wait_for_completion: bool = True, max_wait: int = 300) -> 'CoverageResult':
        """
        Get expression coverage for a genomic region
        
        Args:
            chromosome: Chromosome name (e.g., 'chr1', '1')
            start: Start position
            end: End position  
            strand: Strand to analyze ('positive', 'negative', 'forward', 'reverse')
            split_by_cell_type: Whether to split analysis by cell type
            preserve_samples: Whether to preserve sample-level data (True for API clients)
            wait_for_completion: Whether to wait for completion
            max_wait: Maximum time to wait for completion
            
        Returns:
            CoverageResult object
        """
        # Validate inputs
        if end <= start:
            raise ValueError("End position must be greater than start position")
        
        range_size = end - start
        if range_size > 10_000_000:  # 10Mb limit
            raise ValueError("Range cannot exceed 10Mb")
        
        if strand not in ['forward', 'reverse', 'positive', 'negative']:
            raise ValueError("The strand has to be one of ['forward', 'reverse', 'positive', 'negative']")
        
        if strand == 'negative':
            strand = 'reverse' # standardise name for server
        
        # Submit coverage request
        data = {
            'chromosome': chromosome,
            'start': start,
            'end': end,
            'strand': strand,
            'split_by_cell_type': split_by_cell_type,
            'preserve_samples': preserve_samples
        }
        
        response = self._request('POST', '/api/genome-browser/search', json=data)
        job_id = response['job_id']
        
        logger.info(f"Coverage analysis submitted with job ID: {job_id}")
        
        if not wait_for_completion:
            return CoverageResult({'job_id': job_id, 'status': 'pending'})
        
        # Poll for completion using the RAW endpoint
        start_time = time.time()
        poll_interval = 2
        
        while time.time() - start_time < max_wait:
            try:
                params = {'preserve_samples': 'true'} if preserve_samples else {}
                coverage_response = self._request('GET', f'/api/genome-browser/coverage/{job_id}/raw', params=params)
                
                if coverage_response.get('status') == 'completed':
                    result = CoverageResult(coverage_response)
                    result.set_client(self)  # set client for metadata enrichment
                    return result
                elif coverage_response.get('status') == 'error':
                    error_msg = coverage_response.get('error', 'Unknown error')
                    raise MalvaAPIError(f"Coverage analysis failed: {error_msg}")
                    
            except MalvaAPIError as e:
                if "404" in str(e) or "not found" in str(e).lower():
                    # Still processing
                    pass
                else:
                    raise
            
            time.sleep(poll_interval)
        
        raise MalvaAPIError(f"Coverage analysis timed out after {max_wait} seconds")
    
    def get_sequence_coverage(self, sequence: str, window_size: int = 64,
                            split_by_cell_type: bool = False, preserve_samples: bool = True,
                            wait_for_completion: bool = True, max_wait: int = 300) -> 'CoverageResult':
        """
        Get expression coverage for a DNA sequence
        
        Args:
            sequence: DNA sequence to analyze (A, T, G, C, N only)
            window_size: Size of sliding windows (default: 64)
            split_by_cell_type: Whether to split analysis by cell type
            preserve_samples: Whether to preserve sample-level data
            wait_for_completion: Whether to wait for completion
            max_wait: Maximum time to wait for completion
            
        Returns:
            CoverageResult object with sequence positions instead of genomic coordinates
        """
        # Validate sequence
        if not sequence:
            raise ValueError("Sequence cannot be empty")
        
        sequence = sequence.upper().strip()
        valid_nucleotides = set('ATGCN')
        if not all(c in valid_nucleotides for c in sequence):
            raise ValueError("Sequence contains invalid nucleotides. Only A, T, G, C, N are allowed.")
        
        if len(sequence) < window_size:
            raise ValueError(f"Sequence too short: minimum {window_size} nucleotides required")
        
        if len(sequence) > 50000:  # 50kb limit
            raise ValueError("Sequence too long: maximum 50,000 nucleotides allowed")
        
        # Submit sequence coverage request
        data = {
            'sequence': sequence,
            'window_size': window_size,
            'split_by_cell_type': split_by_cell_type,
            'preserve_samples': preserve_samples
        }
        
        response = self._request('POST', '/api/sequence-coverage/search', json=data)
        job_id = response['job_id']
        
        logger.info(f"Sequence coverage analysis submitted with job ID: {job_id}")
        
        if not wait_for_completion:
            return CoverageResult({'job_id': job_id, 'status': 'pending'})
        
        # Poll for completion using the sequence coverage endpoint
        start_time = time.time()
        poll_interval = 2
        
        while time.time() - start_time < max_wait:
            try:
                params = {'preserve_samples': 'true'} if preserve_samples else {}
                coverage_response = self._request('GET', f'/api/sequence-coverage/coverage/{job_id}/raw', params=params)
                
                if coverage_response.get('status') == 'completed':
                    result = CoverageResult(coverage_response)
                    result.set_client(self)  # set client for metadata enrichment
                    return result
                elif coverage_response.get('status') == 'error':
                    error_msg = coverage_response.get('error', 'Unknown error')
                    raise MalvaAPIError(f"Sequence coverage analysis failed: {error_msg}")
                    
            except MalvaAPIError as e:
                if "404" in str(e) or "not found" in str(e).lower():
                    # Still processing
                    pass
                else:
                    raise
            
            time.sleep(poll_interval)
        
        raise MalvaAPIError(f"Sequence coverage analysis timed out after {max_wait} seconds")

    def download_sequence_coverage_data(self, job_id: str, output_path: Optional[str] = None,
                                    cell_type: str = None, smoothing: int = 0, 
                                    min_expression: float = 0.0) -> Union[str, pd.DataFrame]:
        """
        Download sequence coverage data as TSV file or DataFrame
        
        Args:
            job_id: Job ID from sequence coverage analysis
            output_path: Path to save TSV file. If None, returns as DataFrame
            cell_type: Optional cell type filter
            smoothing: Additional smoothing window size (0 = no additional smoothing)
            min_expression: Minimum expression threshold
            
        Returns:
            File path if saved to disk, or DataFrame if loaded in memory
        """
        # Build parameters for filtering
        params = {}
        if cell_type:
            params['cell_type'] = cell_type
        if smoothing > 0:
            params['smoothing'] = smoothing
        if min_expression > 0:
            params['min_expression'] = min_expression
        
        # Get the TSV data
        url = f'/api/sequence-coverage/coverage-file/{job_id}'
        
        try:
            response = self.session.get(
                urljoin(self.base_url, url.lstrip('/')), 
                params=params,
                timeout=self.timeout
            )
            response.raise_for_status()
            
            if output_path:
                # Save to file
                output_path = Path(output_path)
                output_path.parent.mkdir(parents=True, exist_ok=True)
                
                with open(output_path, 'w') as f:
                    f.write(response.text)
                
                logger.info(f"Sequence coverage data saved to {output_path}")
                return str(output_path)
            else:
                # Load as DataFrame
                try:
                    import io
                    df = pd.read_csv(io.StringIO(response.text), sep='\t')
                    logger.info(f"Loaded sequence coverage data with {len(df)} positions")
                    return df
                except ImportError:
                    raise ImportError("pandas is required to load TSV data into DataFrame. "
                                    "Install with: pip install pandas")
                    
        except requests.exceptions.RequestException as e:
            if "404" in str(e):
                raise MalvaAPIError(f"Sequence coverage job {job_id} not found")
            else:
                raise MalvaAPIError(f"Failed to download sequence coverage data: {e}")
    
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
        
        # Add filters as query parameters
        if filters:
            for key, value in filters.items():
                if isinstance(value, list):
                    params[key] = value
                else:
                    params[key] = [value] if value else []
        
        return self._request('GET', '/samples/search', params=params)
    
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
        response = self._request('GET', '/samples/studies/summary')
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
        return self._request('GET', '/samples/available-filters')


# Convenience functions for common operations
def search_gene(gene: str, base_url: str = "https://malva.mdc-berlin.net") -> SearchResult:
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


def search_sequence(sequence: str, base_url: str = "https://malva.mdc-berlin.net") -> SearchResult:
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

def get_coverage_for_region(chromosome: str, start: int, end: int, 
                           base_url: str = "https://malva.mdc-berlin.net") -> CoverageResult:
    """
    Quick coverage analysis for a genomic region
    
    Args:
        chromosome: Chromosome name
        start: Start position
        end: End position
        base_url: Malva API base URL
        
    Returns:
        CoverageResult object
    """
    client = MalvaClient(base_url)
    return client.get_coverage(chromosome, start, end)

def get_sequence_coverage(sequence: str, window_size: int = 64,
                         base_url: str = "https://malva.mdc-berlin.net") -> 'CoverageResult':
    """
    Quick sequence coverage analysis
    
    Args:
        sequence: DNA sequence to analyze
        window_size: Size of sliding windows (default: 64)
        base_url: Malva API base URL
        
    Returns:
        CoverageResult object
    """
    client = MalvaClient(base_url)
    return client.get_sequence_coverage(sequence, window_size=window_size)