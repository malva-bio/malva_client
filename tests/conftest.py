"""Shared fixtures for malva_client tests."""

import pytest
from unittest.mock import patch, MagicMock

from malva_client.client import MalvaClient


@pytest.fixture
def mock_client():
    """Create a MalvaClient with _test_connection patched out."""
    with patch.object(MalvaClient, '_test_connection'):
        client = MalvaClient(
            base_url="https://test.malva.bio",
            api_token="test-token-123",
        )
    return client


@pytest.fixture
def sample_search_response():
    """Realistic completed search response in the compact format."""
    return {
        'job_id': 'abc-123',
        'status': 'completed',
        'results': {
            'BRCA1': {
                'expression_data': {
                    'samples': [101, 102, 103],
                    'cell_types': ['neuron', 'astrocyte', 'oligodendrocyte'],
                    'columns': ['sample_idx', 'cell_type_idx', 'norm_expr', 'raw_expr', 'cell_count'],
                    'data': [
                        [0, 0, 0.85, 1.2, 150],
                        [0, 1, 0.42, 0.6, 80],
                        [1, 0, 0.91, 1.4, 200],
                        [1, 2, 0.33, 0.5, 60],
                        [2, 1, 0.55, 0.8, 120],
                    ],
                }
            }
        },
    }


@pytest.fixture
def sample_pending_response():
    """Pending search response."""
    return {
        'job_id': 'pending-456',
        'status': 'pending',
    }


@pytest.fixture
def sample_cell_response():
    """Realistic SingleCellResult data."""
    return {
        'job_id': 'cell-789',
        'status': 'completed',
        'results': {
            'seq_1': {
                'cell': [1, 2, 3, 4, 5],
                'expression': [0.5, 0.8, 0.0, 1.2, 0.3],
                'sample': [10, 10, 20, 20, 30],
                'ncells': 5,
                'sequence_length': 100,
            }
        },
    }


@pytest.fixture
def sample_coverage_response():
    """Realistic CoverageResult data."""
    return {
        'job_id': 'cov-101',
        'region': 'chr1:1000-2000(+)',
        'positions': [1000, 1100, 1200, 1300, 1400],
        'cell_types': ['neuron', 'astrocyte'],
        'coverage_matrix': [
            [10, 20, 30, 25, 15],
            [5, 10, 15, 12, 8],
        ],
        'total_windows': 5,
    }


@pytest.fixture
def sample_metadata():
    """Sample metadata dict keyed by sample_id."""
    return {
        101: {
            'specimen_from_organism.organ': 'brain',
            'genus_species': 'Homo sapiens',
            'project.project_core.project_short_name': 'BrainStudy',
        },
        102: {
            'specimen_from_organism.organ': 'liver',
            'genus_species': 'Homo sapiens',
            'project.project_core.project_short_name': 'LiverStudy',
        },
        103: {
            'specimen_from_organism.organ': 'brain',
            'genus_species': 'Mus musculus',
            'project.project_core.project_short_name': 'MouseBrain',
        },
    }
