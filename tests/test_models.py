"""Tests for malva_client.models."""

import pandas as pd
import numpy as np
import pytest
from unittest.mock import MagicMock

from malva_client.models import SearchResult, MalvaDataFrame, SingleCellResult, CoverageResult


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_mock_client():
    return MagicMock()


def _make_expr_df():
    """A small expression DataFrame for MalvaDataFrame tests."""
    return pd.DataFrame({
        'sample_id': [1, 1, 2, 2, 3],
        'cell_type': ['neuron', 'astrocyte', 'neuron', 'oligodendrocyte', 'astrocyte'],
        'norm_expr': [0.8, 0.4, 0.9, 0.3, 0.5],
        'cell_count': [150, 80, 200, 60, 120],
        'gene_sequence': ['BRCA1'] * 5,
    })


# ===========================================================================
# SearchResult
# ===========================================================================

class TestSearchResultFromResponse:
    def test_from_completed_response(self, sample_search_response):
        client = _make_mock_client()
        sr = SearchResult(sample_search_response, client)
        assert len(sr) > 0
        assert 'norm_expr' in sr.df.columns
        assert 'cell_type' in sr.df.columns
        assert 'sample_id' in sr.df.columns

    def test_from_empty_response(self):
        sr = SearchResult({}, _make_mock_client())
        assert len(sr) == 0

    def test_from_none_results(self):
        sr = SearchResult({'results': {}}, _make_mock_client())
        assert len(sr) == 0

    def test_from_pending_response(self, sample_pending_response):
        sr = SearchResult(sample_pending_response, _make_mock_client())
        assert sr.job_id == 'pending-456'
        assert sr.status == 'pending'

    def test_job_id_property(self, sample_search_response):
        sr = SearchResult(sample_search_response, _make_mock_client())
        assert sr.job_id == 'abc-123'

    def test_status_property(self, sample_search_response):
        sr = SearchResult(sample_search_response, _make_mock_client())
        assert sr.status == 'completed'

    def test_total_cells(self, sample_search_response):
        sr = SearchResult(sample_search_response, _make_mock_client())
        # Sum of cell_count column: 150+80+200+60+120 = 610
        assert sr.total_cells == 610

    def test_multiple_gene_sequences(self):
        data = {
            'results': {
                'BRCA1': {
                    'expression_data': {
                        'samples': [1],
                        'cell_types': ['neuron'],
                        'columns': ['sample_idx', 'cell_type_idx', 'norm_expr', 'raw_expr', 'cell_count'],
                        'data': [[0, 0, 0.9, 1.0, 100]],
                    }
                },
                'TP53': {
                    'expression_data': {
                        'samples': [1],
                        'cell_types': ['neuron'],
                        'columns': ['sample_idx', 'cell_type_idx', 'norm_expr', 'raw_expr', 'cell_count'],
                        'data': [[0, 0, 0.5, 0.7, 50]],
                    }
                },
            }
        }
        sr = SearchResult(data, _make_mock_client())
        assert sr.df['gene_sequence'].nunique() == 2

    def test_old_format_backward_compat(self):
        data = {
            'results': {
                'BRCA1': {
                    'cell': [1, 2, 3],
                    'expression': [0.5, 0.8, 0.3],
                    'sample': [10, 10, 20],
                }
            }
        }
        sr = SearchResult(data, _make_mock_client())
        assert len(sr) == 3
        assert 'norm_expr' in sr.df.columns

    def test_returns_search_result_type(self, sample_search_response):
        sr = SearchResult(sample_search_response, _make_mock_client())
        assert isinstance(sr, SearchResult)
        assert isinstance(sr, MalvaDataFrame)


class TestSearchResultRepr:
    def test_repr_does_not_crash(self, sample_search_response):
        sr = SearchResult(sample_search_response, _make_mock_client())
        text = repr(sr)
        assert isinstance(text, str)
        assert len(text) > 0

    def test_repr_html(self, sample_search_response):
        sr = SearchResult(sample_search_response, _make_mock_client())
        html = sr._repr_html_()
        assert isinstance(html, str)
        assert '<' in html  # contains HTML tags


# ===========================================================================
# MalvaDataFrame
# ===========================================================================

class TestMalvaDataFrameFilterBy:
    def test_filter_by_single_field(self):
        df = _make_expr_df()
        mdf = MalvaDataFrame(df, _make_mock_client())
        filtered = mdf.filter_by(cell_type='neuron')
        assert all(filtered.df['cell_type'] == 'neuron')

    def test_filter_by_multiple_fields(self):
        df = _make_expr_df()
        mdf = MalvaDataFrame(df, _make_mock_client())
        filtered = mdf.filter_by(cell_type='neuron', sample_id=1)
        assert len(filtered) == 1


class TestMalvaDataFrameAggregateBy:
    def test_aggregate_by_mean(self):
        df = _make_expr_df()
        mdf = MalvaDataFrame(df, _make_mock_client())
        result = mdf.aggregate_by('cell_type', agg_func='mean')
        assert 'mean_norm_expr' in result.columns
        assert len(result) == 3  # neuron, astrocyte, oligodendrocyte

    def test_aggregate_by_multiple_columns(self):
        df = _make_expr_df()
        mdf = MalvaDataFrame(df, _make_mock_client())
        result = mdf.aggregate_by(['sample_id', 'cell_type'])
        assert 'mean_norm_expr' in result.columns

    def test_aggregate_by_different_funcs(self):
        df = _make_expr_df()
        mdf = MalvaDataFrame(df, _make_mock_client())
        for func in ['median', 'sum', 'count', 'std', 'min', 'max']:
            result = mdf.aggregate_by('cell_type', agg_func=func)
            assert f'{func}_norm_expr' in result.columns

    def test_aggregate_unsupported_func_raises(self):
        df = _make_expr_df()
        mdf = MalvaDataFrame(df, _make_mock_client())
        with pytest.raises(ValueError, match="Unsupported"):
            mdf.aggregate_by('cell_type', agg_func='bogus')

    def test_aggregate_missing_expr_column_raises(self):
        df = pd.DataFrame({'cell_type': ['a'], 'other': [1]})
        mdf = MalvaDataFrame(df, _make_mock_client())
        with pytest.raises(ValueError, match="not found"):
            mdf.aggregate_by('cell_type')


class TestMalvaDataFrameBasics:
    def test_to_pandas(self):
        df = _make_expr_df()
        mdf = MalvaDataFrame(df, _make_mock_client())
        result = mdf.to_pandas()
        assert isinstance(result, pd.DataFrame)
        # Should be a copy
        result['norm_expr'] = 999
        assert mdf.df['norm_expr'].iloc[0] != 999

    def test_len(self):
        df = _make_expr_df()
        mdf = MalvaDataFrame(df, _make_mock_client())
        assert len(mdf) == 5

    def test_getitem(self):
        df = _make_expr_df()
        mdf = MalvaDataFrame(df, _make_mock_client())
        col = mdf['norm_expr']
        assert isinstance(col, pd.Series)
        assert len(col) == 5

    def test_available_fields(self):
        df = _make_expr_df()
        mdf = MalvaDataFrame(df, _make_mock_client())
        fields = mdf.available_fields()
        assert isinstance(fields, dict)
        assert 'expression' in fields or 'basic' in fields

    def test_unique_values(self):
        df = _make_expr_df()
        mdf = MalvaDataFrame(df, _make_mock_client())
        vals = mdf.unique_values('cell_type')
        assert 'neuron' in vals
        assert 'astrocyte' in vals

    def test_unique_values_with_limit(self):
        df = _make_expr_df()
        mdf = MalvaDataFrame(df, _make_mock_client())
        vals = mdf.unique_values('cell_type', limit=1)
        assert len(vals) == 1

    def test_value_counts(self):
        df = _make_expr_df()
        mdf = MalvaDataFrame(df, _make_mock_client())
        vc = mdf.value_counts('cell_type')
        assert isinstance(vc, pd.Series)
        # neuron=2, astrocyte=2, oligodendrocyte=1
        assert vc.sum() == 5

    def test_value_counts_with_limit(self):
        df = _make_expr_df()
        mdf = MalvaDataFrame(df, _make_mock_client())
        vc = mdf.value_counts('cell_type', limit=1)
        assert len(vc) == 1


# ===========================================================================
# SingleCellResult
# ===========================================================================

class TestSingleCellResult:
    def test_from_completed_response(self, sample_cell_response):
        scr = SingleCellResult(sample_cell_response)
        assert scr.is_completed
        assert scr.cell_count == 5
        assert scr.sample_count == 3  # samples 10, 20, 30

    def test_to_dataframe(self, sample_cell_response):
        scr = SingleCellResult(sample_cell_response)
        df = scr.to_dataframe()
        assert list(df.columns) == ['cell_id', 'expression', 'sample_id']
        assert len(df) == 5

    def test_filter_by_expression(self, sample_cell_response):
        scr = SingleCellResult(sample_cell_response)
        filtered = scr.filter_by_expression(min_expression=0.5)
        assert filtered.cell_count == 3  # 0.5, 0.8, 1.2

    def test_get_expression_stats(self, sample_cell_response):
        scr = SingleCellResult(sample_cell_response)
        stats = scr.get_expression_stats()
        assert 'mean' in stats
        assert 'median' in stats
        assert 'std' in stats
        assert 'min' in stats
        assert 'max' in stats
        assert stats['total_cells'] == 5
        assert stats['min'] == 0.0
        assert stats['max'] == 1.2

    def test_aggregate_by_sample(self, sample_cell_response):
        scr = SingleCellResult(sample_cell_response)
        agg = scr.aggregate_by_sample()
        assert 'sample_id' in agg.columns
        assert 'mean_expression' in agg.columns
        assert len(agg) == 3  # 3 unique samples

    def test_pending_result(self):
        scr = SingleCellResult({'job_id': 'p1', 'status': 'pending'})
        assert scr.is_pending
        assert not scr.has_results
        assert scr.cell_count == 0

    def test_repr_completed(self, sample_cell_response):
        scr = SingleCellResult(sample_cell_response)
        text = repr(scr)
        assert 'cells=5' in text

    def test_repr_pending(self):
        scr = SingleCellResult({'job_id': 'p1', 'status': 'pending'})
        text = repr(scr)
        assert 'pending' in text

    def test_len(self, sample_cell_response):
        scr = SingleCellResult(sample_cell_response)
        assert len(scr) == 5


# ===========================================================================
# CoverageResult
# ===========================================================================

class TestCoverageResult:
    def test_from_response(self, sample_coverage_response):
        cr = CoverageResult(sample_coverage_response)
        assert cr.job_id == 'cov-101'
        assert cr.region == 'chr1:1000-2000(+)'
        assert len(cr.cell_types) == 2
        assert len(cr.positions) == 5
        assert cr.total_windows == 5

    def test_to_dataframe(self, sample_coverage_response):
        cr = CoverageResult(sample_coverage_response)
        df = cr.to_dataframe()
        assert list(df.columns) == ['neuron', 'astrocyte']
        assert len(df) == 5
        assert df.index.name == 'position'

    def test_properties(self, sample_coverage_response):
        cr = CoverageResult(sample_coverage_response)
        assert cr.job_id == 'cov-101'
        assert cr.region == 'chr1:1000-2000(+)'
        assert cr.cell_types == ['neuron', 'astrocyte']
        assert cr.total_windows == 5

    def test_len(self, sample_coverage_response):
        cr = CoverageResult(sample_coverage_response)
        assert len(cr) == 5

    def test_repr(self, sample_coverage_response):
        cr = CoverageResult(sample_coverage_response)
        text = repr(cr)
        assert 'cov-101' in text
        assert 'positions=5' in text

    def test_empty_coverage(self):
        cr = CoverageResult({'job_id': 'empty'})
        df = cr.to_dataframe()
        assert df.empty
