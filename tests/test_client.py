"""Tests for malva_client.client."""

import json
import time
from pathlib import Path
from unittest.mock import patch, MagicMock, PropertyMock

import pytest
import responses
from responses import matchers

from malva_client.client import MalvaClient
from malva_client.models import SearchResult, SingleCellResult, CoverageResult, CoexpressionResult, UMAPCoordinates
from malva_client.exceptions import (
    MalvaAPIError,
    AuthenticationError,
    QuotaExceededError,
    SearchError,
)


BASE = "https://test.malva.bio"


# ===========================================================================
# Constructor & Connection
# ===========================================================================

class TestConstructor:
    @responses.activate
    def test_init_with_explicit_params(self):
        responses.add(responses.GET, f"{BASE}/health", json={"status": "ok"})
        responses.add(responses.GET, f"{BASE}/api/quota-status",
                      json={"account_type": "test"})

        client = MalvaClient(base_url=BASE, api_token="tok-123")
        assert client.base_url == BASE
        assert client.api_token == "tok-123"

    def test_init_loads_config_when_no_url(self):
        mock_config = MagicMock()
        mock_config.server_url = "https://cfg.malva.bio"
        mock_config.api_token = "cfg-tok"
        mock_config.verify_ssl = True

        with patch('malva_client.client.Config.load', return_value=mock_config), \
             patch.object(MalvaClient, '_test_connection'):
            client = MalvaClient()
            assert client.base_url == "https://cfg.malva.bio"

    def test_init_raises_without_url(self):
        mock_config = MagicMock()
        mock_config.server_url = None
        mock_config.api_token = None
        mock_config.verify_ssl = True

        with patch('malva_client.client.Config.load', return_value=mock_config):
            with pytest.raises(ValueError, match="base_url must be provided"):
                MalvaClient()

    @responses.activate
    def test_init_connection_failure(self):
        responses.add(responses.GET, f"{BASE}/health",
                      json={"error": "down"}, status=500)

        with pytest.raises(MalvaAPIError):
            MalvaClient(base_url=BASE, api_token="tok")


# ===========================================================================
# Search Operations
# ===========================================================================

class TestSearch:
    @responses.activate
    def test_search_completed_immediately(self, mock_client):
        responses.add(
            responses.POST, f"{BASE}/search",
            json={
                'status': 'completed',
                'job_id': 'j1',
                'results': {
                    'BRCA1': {
                        'expression_data': {
                            'samples': [1],
                            'cell_types': ['neuron'],
                            'columns': ['sample_idx', 'cell_type_idx', 'norm_expr', 'raw_expr', 'cell_count'],
                            'data': [[0, 0, 0.9, 1.0, 100]],
                        }
                    }
                },
            },
        )
        result = mock_client.search("BRCA1")
        assert isinstance(result, SearchResult)

    @responses.activate
    def test_search_polls_until_complete(self, mock_client):
        responses.add(
            responses.POST, f"{BASE}/search",
            json={'job_id': 'j1', 'status': 'running'},
        )
        responses.add(
            responses.GET, f"{BASE}/search/j1",
            json={'status': 'pending'},
        )
        responses.add(
            responses.GET, f"{BASE}/search/j1",
            json={'status': 'completed', 'results': {}},
        )

        result = mock_client.search("BRCA1", poll_interval=0, max_wait=10)
        assert isinstance(result, SearchResult)

    @responses.activate
    def test_search_timeout(self, mock_client):
        responses.add(
            responses.POST, f"{BASE}/search",
            json={'job_id': 'j1', 'status': 'running'},
        )
        # Always return pending
        responses.add(
            responses.GET, f"{BASE}/search/j1",
            json={'status': 'pending'},
        )

        with pytest.raises(MalvaAPIError, match="timed out"):
            mock_client.search("BRCA1", poll_interval=0, max_wait=0.01)

    @responses.activate
    def test_search_server_error(self, mock_client):
        responses.add(
            responses.POST, f"{BASE}/search",
            json={'job_id': 'j1', 'status': 'running'},
        )
        responses.add(
            responses.GET, f"{BASE}/search/j1",
            json={'status': 'error', 'error': 'internal failure'},
        )

        with pytest.raises(MalvaAPIError, match="internal failure"):
            mock_client.search("BRCA1", poll_interval=0, max_wait=10)

    @responses.activate
    def test_search_with_parameters(self, mock_client):
        def check_body(request):
            body = json.loads(request.body)
            assert body['window_size'] == 24
            assert body['threshold'] == 0.5
            assert body['stranded'] is True
            return (200, {}, json.dumps({
                'status': 'completed', 'job_id': 'j1',
                'results': {},
            }))

        responses.add_callback(responses.POST, f"{BASE}/search", callback=check_body)
        mock_client.search("BRCA1", window_size=24, threshold=0.5, stranded=True)

    def test_search_invalid_window_size(self, mock_client):
        with pytest.raises(ValueError, match="window_size"):
            mock_client.search("BRCA1", window_size=0)

    def test_search_invalid_threshold(self, mock_client):
        with pytest.raises(ValueError, match="threshold"):
            mock_client.search("BRCA1", threshold=1.5)

    @responses.activate
    def test_search_cells(self, mock_client):
        responses.add(
            responses.POST, f"{BASE}/search",
            json={
                'status': 'completed',
                'job_id': 'j1',
                'results': {
                    'seq_1': {
                        'cell': [1, 2],
                        'expression': [0.5, 0.8],
                        'sample': [10, 10],
                        'ncells': 2,
                    }
                },
            },
        )
        result = mock_client.search_cells("BRCA1")
        assert isinstance(result, SingleCellResult)

    def test_search_sequence_validation(self, mock_client):
        with pytest.raises(ValueError, match="invalid nucleotides"):
            mock_client.search_sequence("ATGCXYZ")

    def test_search_sequence_too_long(self, mock_client):
        with pytest.raises(ValueError, match="500kb"):
            mock_client.search_sequence("A" * 500_001)

    @responses.activate
    def test_search_sequence_valid(self, mock_client):
        responses.add(
            responses.POST, f"{BASE}/search",
            json={'status': 'completed', 'job_id': 'j1', 'results': {}},
        )
        result = mock_client.search_sequence("ATGCN")
        assert isinstance(result, SearchResult)

    @responses.activate
    def test_search_returns_search_result_type(self, mock_client):
        responses.add(
            responses.POST, f"{BASE}/search",
            json={'status': 'completed', 'job_id': 'j1', 'results': {}},
        )
        result = mock_client.search("BRCA1")
        assert isinstance(result, SearchResult)


# ===========================================================================
# Batch Search
# ===========================================================================

class TestBatchSearch:
    @responses.activate
    def test_search_sequences_fasta_payload(self, mock_client):
        def check_fasta(request):
            body = json.loads(request.body)
            assert 'fasta_sequences' in body
            assert '>seq_1' in body['fasta_sequences']
            assert '>seq_2' in body['fasta_sequences']
            return (200, {}, json.dumps({
                'status': 'completed', 'job_id': 'j1', 'results': {},
            }))

        responses.add_callback(responses.POST, f"{BASE}/search", callback=check_fasta)
        mock_client.search_sequences(["ATGC", "GCTA"])

    @responses.activate
    def test_search_sequences_completed(self, mock_client):
        responses.add(
            responses.POST, f"{BASE}/search",
            json={'status': 'completed', 'job_id': 'j1', 'results': {}},
        )
        result = mock_client.search_sequences(["ATGC", "GCTA"])
        assert isinstance(result, SearchResult)

    @responses.activate
    def test_search_sequences_async(self, mock_client):
        responses.add(
            responses.POST, f"{BASE}/search",
            json={'job_id': 'j1', 'status': 'pending'},
        )
        result = mock_client.search_sequences(["ATGC"], wait_for_completion=False)
        assert isinstance(result, SearchResult)
        assert result.status == 'pending'

    def test_search_sequences_empty_list(self, mock_client):
        with pytest.raises(ValueError, match="At least one sequence"):
            mock_client.search_sequences([])

    def test_search_sequences_invalid_nucleotides(self, mock_client):
        with pytest.raises(ValueError, match="invalid nucleotides"):
            mock_client.search_sequences(["ATGCXYZ"])

    def test_search_sequences_exceeds_100kb(self, mock_client):
        with pytest.raises(ValueError, match="100,000"):
            mock_client.search_sequences(["A" * 100_001])

    @responses.activate
    def test_search_genes_completed(self, mock_client):
        def check_query(request):
            body = json.loads(request.body)
            assert 'BRCA1' in body['query']
            assert 'TP53' in body['query']
            return (200, {}, json.dumps({
                'status': 'completed', 'job_id': 'j1', 'results': {},
            }))

        responses.add_callback(responses.POST, f"{BASE}/search", callback=check_query)
        result = mock_client.search_genes(["BRCA1", "TP53"])
        assert isinstance(result, SearchResult)

    def test_search_genes_too_many(self, mock_client):
        genes = [f"GENE{i}" for i in range(11)]
        with pytest.raises(ValueError, match="Too many genes"):
            mock_client.search_genes(genes)


# ===========================================================================
# Job Management
# ===========================================================================

class TestJobManagement:
    @responses.activate
    def test_submit_search(self, mock_client):
        responses.add(
            responses.POST, f"{BASE}/search",
            json={'job_id': 'j1', 'status': 'pending'},
        )
        with patch('malva_client.storage.get_storage') as mock_gs:
            mock_storage = MagicMock()
            mock_storage.save_search.return_value = 1
            mock_gs.return_value = mock_storage
            job_id = mock_client.submit_search("BRCA1")
            assert job_id == 'j1'

    @responses.activate
    def test_get_job_status(self, mock_client):
        responses.add(
            responses.GET, f"{BASE}/search/j1",
            json={'status': 'running', 'job_id': 'j1'},
        )
        with patch('malva_client.storage.get_storage') as mock_gs:
            mock_gs.return_value = MagicMock()
            result = mock_client.get_job_status("j1")
            assert result['status'] == 'running'

    @responses.activate
    def test_get_job_results_completed(self, mock_client):
        responses.add(
            responses.GET, f"{BASE}/search/j1",
            json={'status': 'completed', 'job_id': 'j1', 'results': {}},
        )
        with patch('malva_client.storage.get_storage') as mock_gs:
            mock_storage = MagicMock()
            mock_storage.get_job_status.return_value = None
            mock_gs.return_value = mock_storage
            result = mock_client.get_job_results("j1")
            assert isinstance(result, SearchResult)

    @responses.activate
    def test_get_job_results_still_pending(self, mock_client):
        responses.add(
            responses.GET, f"{BASE}/search/j1",
            json={'status': 'pending', 'job_id': 'j1'},
        )
        with patch('malva_client.storage.get_storage') as mock_gs:
            mock_storage = MagicMock()
            mock_storage.get_job_status.return_value = None
            mock_gs.return_value = mock_storage
            with pytest.raises(SearchError, match="still pending"):
                mock_client.get_job_results("j1")

    @responses.activate
    def test_wait_for_job(self, mock_client):
        responses.add(
            responses.GET, f"{BASE}/search/j1",
            json={'status': 'pending', 'job_id': 'j1'},
        )
        responses.add(
            responses.GET, f"{BASE}/search/j1",
            json={'status': 'completed', 'job_id': 'j1', 'results': {}},
        )
        with patch('malva_client.storage.get_storage') as mock_gs:
            mock_storage = MagicMock()
            mock_storage.get_job_status.return_value = None
            mock_gs.return_value = mock_storage
            result = mock_client.wait_for_job("j1", poll_interval=0, max_wait=10)
            assert isinstance(result, SearchResult)

    @responses.activate
    def test_wait_for_job_timeout(self, mock_client):
        responses.add(
            responses.GET, f"{BASE}/search/j1",
            json={'status': 'pending', 'job_id': 'j1'},
        )
        with patch('malva_client.storage.get_storage') as mock_gs:
            mock_gs.return_value = MagicMock()
            with pytest.raises(MalvaAPIError, match="timed out"):
                mock_client.wait_for_job("j1", poll_interval=0, max_wait=0.01)

    @responses.activate
    def test_cancel_job(self, mock_client):
        responses.add(
            responses.DELETE, f"{BASE}/search/j1",
            json={'success': True},
        )
        with patch('malva_client.storage.get_storage') as mock_gs:
            mock_gs.return_value = MagicMock()
            assert mock_client.cancel_job("j1") is True

    @responses.activate
    def test_list_jobs(self, mock_client):
        with patch('malva_client.storage.get_storage') as mock_gs:
            mock_storage = MagicMock()
            mock_storage.get_search_history.return_value = [
                {'job_id': 'j1', 'query': 'BRCA1', 'status': 'completed',
                 'timestamp': '2024-01-01', 'results_count': 5, 'error_message': None},
            ]
            mock_gs.return_value = mock_storage
            jobs = mock_client.list_jobs()
            assert len(jobs) == 1
            assert jobs[0]['job_id'] == 'j1'


# ===========================================================================
# Coverage
# ===========================================================================

class TestCoverage:
    @responses.activate
    def test_get_coverage(self, mock_client):
        responses.add(
            responses.POST, f"{BASE}/api/genome-browser/search",
            json={'job_id': 'cov1'},
        )
        responses.add(
            responses.GET, f"{BASE}/api/genome-browser/coverage/cov1",
            json={
                'status': 'completed',
                'coverage_data': {
                    'positions': [100, 200],
                    'cell_types': ['neuron'],
                    'coverage_matrix': [[5, 10]],
                },
            },
        )
        result = mock_client.get_coverage('chr1', 100, 200, poll_interval=0)
        assert isinstance(result, CoverageResult)

    @responses.activate
    def test_get_sequence_coverage(self, mock_client):
        responses.add(
            responses.POST, f"{BASE}/api/genome-browser/search-sequence",
            json={'job_id': 'cov2'},
        )
        responses.add(
            responses.GET, f"{BASE}/api/genome-browser/coverage/cov2",
            json={
                'status': 'completed',
                'coverage_data': {
                    'positions': [0, 1],
                    'cell_types': ['neuron'],
                    'coverage_matrix': [[3, 7]],
                },
            },
        )
        result = mock_client.get_sequence_coverage("ATGC", poll_interval=0)
        assert isinstance(result, CoverageResult)

    @responses.activate
    def test_get_coverage_data(self, mock_client):
        responses.add(
            responses.GET, f"{BASE}/api/genome-browser/coverage/cov3",
            json={
                'status': 'completed',
                'coverage_data': {
                    'positions': [0],
                    'cell_types': ['astrocyte'],
                    'coverage_matrix': [[5]],
                },
            },
        )
        result = mock_client.get_coverage_data("cov3")
        assert isinstance(result, CoverageResult)

    @responses.activate
    def test_get_coverage_filter_options(self, mock_client):
        responses.add(
            responses.GET, f"{BASE}/api/genome-browser/filter-options/cov4",
            json={'organ': ['brain', 'liver']},
        )
        result = mock_client.get_coverage_filter_options("cov4")
        assert 'organ' in result

    @responses.activate
    def test_download_coverage_wig(self, mock_client, tmp_path):
        output = tmp_path / "out.wig"
        responses.add(
            responses.GET, f"{BASE}/api/genome-browser/coverage-file/cov5",
            body="track type=wiggle_0\nchr1 100 200 5\n",
            status=200,
        )
        result = mock_client.download_coverage_wig("cov5", str(output), poll_interval=0)
        assert Path(result).exists()
        assert "wiggle" in Path(result).read_text()


# ===========================================================================
# Dataset Discovery
# ===========================================================================

class TestDatasetDiscovery:
    @responses.activate
    def test_get_datasets_hierarchy(self, mock_client):
        responses.add(
            responses.GET, f"{BASE}/samples/api/datasets/hierarchy",
            json={'datasets': [{'id': 'd1'}]},
        )
        result = mock_client.get_datasets_hierarchy()
        assert 'datasets' in result

    @responses.activate
    def test_get_overview_stats(self, mock_client):
        responses.add(
            responses.GET, f"{BASE}/samples/api/overview",
            json={'total_samples': 1000, 'total_studies': 50},
        )
        result = mock_client.get_overview_stats()
        assert result['total_samples'] == 1000

    @responses.activate
    def test_get_sample_metadata(self, mock_client):
        responses.add(
            responses.POST, f"{BASE}/api/sample-metadata",
            json={'metadata': {'101': {'organ': 'brain'}, '102': {'organ': 'liver'}}},
        )
        result = mock_client.get_sample_metadata([101, 102])
        assert 101 in result
        assert result[101]['organ'] == 'brain'

    @responses.activate
    def test_get_studies(self, mock_client):
        responses.add(
            responses.GET, f"{BASE}/samples/api/studies",
            json={'studies': [{'name': 'Study1', 'count': 10}]},
        )
        import pandas as pd
        result = mock_client.get_studies()
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 1


# ===========================================================================
# Error Handling
# ===========================================================================

# ===========================================================================
# Coexpression
# ===========================================================================

class TestCoexpressionMethods:
    @responses.activate
    def test_get_umap_coordinates(self, mock_client):
        responses.add(
            responses.GET, f"{BASE}/api/coexpression/umap/human_cortex",
            json={
                'x': [1.0, 2.0], 'y': [3.0, 4.0],
                'metacell_ids': [1, 2],
                'clusters': ['c1', 'c2'],
            },
        )
        result = mock_client.get_umap_coordinates("human_cortex")
        assert isinstance(result, UMAPCoordinates)
        assert len(result) == 2

    @responses.activate
    def test_get_coexpression(self, mock_client):
        def check_body(request):
            body = json.loads(request.body)
            assert body['job_id'] == 'j1'
            assert body['dataset_id'] == 'human_cortex'
            # Full query should NOT have include_* flags set to False
            assert 'include_go_enrichment' not in body
            return (200, {}, json.dumps({
                'dataset_id': 'human_cortex',
                'correlated_genes': [{'gene': 'FOX', 'correlation': 0.9, 'p_value': 1e-5}],
                'n_query_cells': 100,
                'n_mapped_metacells': 10,
            }))

        responses.add_callback(
            responses.POST, f"{BASE}/api/coexpression/query-by-job",
            callback=check_body,
        )
        result = mock_client.get_coexpression("j1", "human_cortex")
        assert isinstance(result, CoexpressionResult)
        assert result.dataset_id == 'human_cortex'
        assert len(result.correlated_genes) == 1

    @responses.activate
    def test_get_coexpression_genes(self, mock_client):
        def check_body(request):
            body = json.loads(request.body)
            assert body['job_id'] == 'j1'
            assert body['dataset_id'] == 'human_cortex'
            assert body['include_go_enrichment'] is False
            assert body['include_umap_scores'] is False
            assert body['include_cell_type_enrichment'] is False
            assert body['include_tissue_breakdown'] is False
            return (200, {}, json.dumps({
                'dataset_id': 'human_cortex',
                'correlated_genes': [{'gene': 'FOX', 'correlation': 0.9, 'p_value': 1e-5}],
                'n_query_cells': 100,
                'n_mapped_metacells': 10,
            }))

        responses.add_callback(
            responses.POST, f"{BASE}/api/coexpression/query-by-job",
            callback=check_body,
        )
        result = mock_client.get_coexpression_genes("j1", "human_cortex")
        assert isinstance(result, CoexpressionResult)
        assert len(result.correlated_genes) == 1


# ===========================================================================
# Error Handling
# ===========================================================================

class TestErrorHandling:
    @responses.activate
    def test_request_401_raises_auth_error(self, mock_client):
        responses.add(
            responses.GET, f"{BASE}/health",
            json={'error': 'unauthorized'}, status=401,
        )
        with pytest.raises(AuthenticationError):
            mock_client._request('GET', '/health')

    @responses.activate
    def test_request_429_raises_quota_error(self, mock_client):
        responses.add(
            responses.GET, f"{BASE}/health",
            json={'error': 'quota exceeded', 'account_type': 'free'}, status=429,
        )
        with pytest.raises(QuotaExceededError):
            mock_client._request('GET', '/health')

    @responses.activate
    def test_request_404_raises_api_error(self, mock_client):
        responses.add(
            responses.GET, f"{BASE}/nonexistent",
            json={'detail': 'not found'}, status=404,
        )
        with pytest.raises(MalvaAPIError, match="not found"):
            mock_client._request('GET', '/nonexistent')

    @responses.activate
    def test_request_connection_error(self, mock_client):
        import requests as req
        responses.add(
            responses.GET, f"{BASE}/health",
            body=req.exceptions.ConnectionError("connection refused"),
        )
        with pytest.raises(MalvaAPIError, match="connect"):
            mock_client._request('GET', '/health')
