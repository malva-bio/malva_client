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


@pytest.fixture
def sample_coexpression_response():
    """Realistic CoexpressionResult data."""
    return {
        'dataset_id': 'human_cortex',
        'correlated_genes': [
            {'gene': 'FOXP3', 'correlation': 0.95, 'p_value': 1e-10},
            {'gene': 'IL2RA', 'correlation': 0.88, 'p_value': 1e-8},
            {'gene': 'CTLA4', 'correlation': 0.82, 'p_value': 1e-7},
            {'gene': 'IKZF2', 'correlation': 0.75, 'p_value': 1e-6},
            {'gene': 'TNFRSF18', 'correlation': 0.70, 'p_value': 1e-5},
        ],
        'umap_scores': {
            'metacell_ids': [1, 2, 3, 4],
            'x': [1.0, 2.0, 3.0, 4.0],
            'y': [5.0, 6.0, 7.0, 8.0],
            'positive_fraction': [0.8, 0.1, 0.5, 0.3],
        },
        'go_enrichment': [
            {'go_id': 'GO:0001', 'name': 'immune regulation', 'fdr': 1e-5},
            {'go_id': 'GO:0002', 'name': 'T cell activation', 'fdr': 1e-3},
            {'go_id': 'GO:0003', 'name': 'cytokine signaling', 'fdr': 0.01},
        ],
        'cell_type_enrichment': [
            {'cell_type': 'T regulatory', 'enrichment_score': 0.9, 'p_value': 1e-8},
            {'cell_type': 'CD4+ T cell', 'enrichment_score': 0.6, 'p_value': 1e-4},
        ],
        'tissue_breakdown': [
            {'tissue': 'blood', 'fraction': 0.6},
            {'tissue': 'lymph node', 'fraction': 0.3},
            {'tissue': 'spleen', 'fraction': 0.1},
        ],
        'n_query_cells': 5000,
        'n_mapped_metacells': 120,
    }


@pytest.fixture
def sample_umap_compact():
    """Realistic compact UMAP coordinate data."""
    return {
        'x': [1.0, 2.0, 3.0, 4.0, 5.0],
        'y': [10.0, 20.0, 30.0, 40.0, 50.0],
        'metacell_id': [101, 102, 103, 104, 105],
        'n_cells': [50, 80, 120, 60, 90],
        'sample': ['s1', 's1', 's2', 's2', 's3'],
        'cluster': ['c1', 'c1', 'c2', 'c2', 'c3'],
    }
