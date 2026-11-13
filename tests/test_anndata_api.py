"""Offline tests for the AnnData-friendly client API."""

import warnings
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from malva_client.models import (
    CellExpressionStream,
    MetacellMatchResult,
    SearchResult,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _search_response():
    columns = [
        "sample_idx",
        "cell_type_idx",
        "norm_expr",
        "kpt_expr",
        "cell_count",
        "fraction_positive",
    ]
    expression = {
        "samples": [2, 1],
        "cell_types": ["b", "a"],
        "columns": columns,
        "data": [
            [0, 0, 0.5, 1.0, 10, 0.5],
            [1, 1, 0.7, 2.0, 20, 0.7],
        ],
    }
    return {
        "job_id": "job-1",
        "status": "completed",
        "results": {
            "G2": {"expression_data": dict(expression)},
            "G1": {"expression_data": dict(expression)},
        },
    }


def _stream_frames():
    return [
        pd.DataFrame(
            {
                "sample_id": [2, 2],
                "cell_id": [20, 20],
                "store_index": [0, 0],
                "gene": ["G2", "G1"],
                "value": [0.5, 1.0],
            }
        ),
        pd.DataFrame(
            {
                "sample_id": [1],
                "cell_id": [10],
                "store_index": [0],
                "gene": ["G2"],
                "value": [0.7],
            }
        ),
    ]


def _metacell_response():
    return {
        "dataset_id": "ds1",
        "metric": "cosine",
        "top_n": 50,
        "n_metacells": 100,
        "n_clusters": 2,
        "clusters": [
            {
                "name": "a",
                "neighbors": [
                    {
                        "metacell_id": "m1",
                        "sample_id": 1,
                        "sample_uuid": "uuid-1",
                        "cell_type": "T",
                        "n_cells": 10,
                        "similarity": 0.9,
                        "distance": 0.1,
                        "enrichment_score": 1.0,
                        "rank": 1,
                    },
                    {
                        "metacell_id": "m2",
                        "sample_id": 2,
                        "sample_uuid": "uuid-2",
                        "cell_type": "T",
                        "n_cells": 12,
                        "similarity": 0.4,
                        "distance": 0.6,
                        "enrichment_score": 0.5,
                        "rank": 2,
                    },
                ],
            },
            {
                "name": "b",
                "neighbors": [
                    {
                        "metacell_id": "m3",
                        "sample_id": 1,
                        "sample_uuid": "uuid-1",
                        "cell_type": "B",
                        "n_cells": 8,
                        "similarity": 0.8,
                        "distance": 0.2,
                        "enrichment_score": 2.0,
                        "rank": 1,
                    },
                    {
                        "metacell_id": "m4",
                        "sample_id": 3,
                        "sample_uuid": "uuid-3",
                        "cell_type": "T",
                        "n_cells": 7,
                        "similarity": 0.2,
                        "distance": 0.8,
                        "enrichment_score": 0.1,
                        "rank": 2,
                    },
                ],
            },
        ],
    }


def _capture_response(monkeypatch, client, response):
    captured = {}

    def fake(method, endpoints, **kwargs):
        captured["json"] = kwargs["json"]
        return response, endpoints[0]

    monkeypatch.setattr(client, "_request_first_available", fake)
    return captured


# ---------------------------------------------------------------------------
# SearchResult
# ---------------------------------------------------------------------------


class TestSearchResultAnnData:
    def test_to_anndata_fixed_layout(self, mock_client):
        result = SearchResult(_search_response(), mock_client)

        adata = result.to_anndata()

        assert adata.shape == (2, 2)
        assert isinstance(adata.X, sparse.csr_matrix)
        assert adata.X.dtype == np.float32
        assert adata.obs_names.tolist() == ["1:a", "2:b"]
        assert adata.var_names.tolist() == ["G1", "G2"]
        assert "sample_id" in adata.obs.columns
        assert "cell_type" in adata.obs.columns
        assert adata.uns["malva"]["layout"] == "sample_x_cell_type_by_gene"
        assert adata.uns["malva"]["job_id"] == "job-1"

    def test_deterministic_var_order(self, mock_client):
        first = SearchResult(_search_response(), mock_client).to_anndata()
        second = SearchResult(_search_response(), mock_client).to_anndata()
        assert first.var_names.tolist() == second.var_names.tolist()

    def test_empty_result_returns_empty_anndata(self, mock_client):
        adata = SearchResult({}, mock_client).to_anndata()
        assert adata.shape == (0, 0)
        assert isinstance(adata.X, sparse.csr_matrix)

    def test_enrich_fetches_only_missing_samples(self, mock_client, monkeypatch):
        result = SearchResult(_search_response(), mock_client)
        result._sample_metadata = {1: {"organ": "brain"}}
        calls = []

        def fake_metadata(sample_ids):
            calls.append(list(sample_ids))
            return {2: {"organ": "lung"}}

        monkeypatch.setattr(mock_client, "get_sample_metadata", fake_metadata)

        returned = result.enrich()

        assert returned is result
        assert calls == [[2]]
        assert set(result.df["organ"].tolist()) == {"brain", "lung"}

    def test_enrich_is_idempotent(self, mock_client, monkeypatch):
        result = SearchResult(_search_response(), mock_client)
        result._sample_metadata = {1: {"organ": "brain"}, 2: {"organ": "lung"}}
        calls = []
        monkeypatch.setattr(
            mock_client,
            "get_sample_metadata",
            lambda sample_ids: calls.append(list(sample_ids)) or {},
        )

        result.enrich()
        result.enrich()

        assert calls == []
        assert result.df.columns.tolist().count("organ") == 1

    def test_enrich_with_metadata_is_alias(self, mock_client):
        result = SearchResult(_search_response(), mock_client)
        result._sample_metadata = {1: {"organ": "brain"}, 2: {"organ": "lung"}}
        assert result.enrich_with_metadata() is result

    def test_to_pandas_alias_matches_to_dataframe(self, mock_client):
        result = SearchResult(_search_response(), mock_client)
        pd.testing.assert_frame_equal(result.to_dataframe(), result.to_pandas())

    def test_write_h5ad(self, mock_client, tmp_path):
        result = SearchResult(_search_response(), mock_client)
        path = tmp_path / "search.h5ad"
        assert result.write_h5ad(path) is result
        adata = ad.read_h5ad(path, backed="r")
        try:
            assert adata.shape == (2, 2)
        finally:
            adata.file.close()


# ---------------------------------------------------------------------------
# CellExpressionStream
# ---------------------------------------------------------------------------


class TestCellExpressionStream:
    def test_iteration_still_yields_dataframes(self):
        stream = CellExpressionStream(_stream_frames())
        batches = list(stream)
        assert len(batches) == 2
        assert all(isinstance(batch, pd.DataFrame) for batch in batches)
        assert batches[0].columns.tolist() == [
            "sample_id",
            "cell_id",
            "store_index",
            "gene",
            "value",
        ]

    def test_to_anndata_cells_by_genes(self):
        stream = CellExpressionStream(_stream_frames(), dataset_id="ds1")

        adata = stream.to_anndata()

        assert adata.shape == (2, 2)
        assert adata.obs_names.tolist() == ["1:10", "2:20"]
        assert adata.var_names.tolist() == ["G1", "G2"]
        assert isinstance(adata.X, sparse.csr_matrix)
        assert adata.X.dtype == np.float32
        assert adata.uns["malva"]["source"] == "cell_expression_stream"
        # Duplicate (cell, gene) entries are summed, not duplicated as observations.
        np.testing.assert_allclose(adata.X.toarray(), [[0.0, 0.7], [1.0, 0.5]])

    def test_obs_names_are_unique(self):
        frame = pd.DataFrame(
            {
                "sample_id": [1, 1, 1],
                "cell_id": [10, 10, 11],
                "store_index": [0, 0, 0],
                "gene": ["G1", "G2", "G1"],
                "value": [1.0, 2.0, 3.0],
            }
        )
        adata = CellExpressionStream([frame]).to_anndata()
        assert adata.obs_names.is_unique
        assert adata.obs_names.tolist() == ["1:10", "1:11"]

    def test_iter_anndata_supports_concat(self):
        frames = [
            pd.DataFrame(
                {
                    "sample_id": [1],
                    "cell_id": [10],
                    "store_index": [0],
                    "gene": ["G2"],
                    "value": [1.0],
                }
            ),
            pd.DataFrame(
                {
                    "sample_id": [2],
                    "cell_id": [20],
                    "store_index": [0],
                    "gene": ["G1"],
                    "value": [2.0],
                }
            ),
        ]
        parts = list(CellExpressionStream(frames).iter_anndata())
        assert len(parts) == 2
        assert all(part.var_names.tolist() == sorted(part.var_names) for part in parts)

        combined = ad.concat(parts, join="outer", fill_value=0)
        assert combined.shape == (2, 2)
        assert combined.obs_names.is_unique
        assert combined.var_names.tolist() == ["G1", "G2"]

    def test_write_h5ad_is_readable_backed(self, tmp_path):
        stream = CellExpressionStream(_stream_frames(), dataset_id="ds1")
        path = tmp_path / "cells.h5ad"

        returned = stream.write_h5ad(path)

        assert returned is stream
        adata = ad.read_h5ad(path, backed="r")
        try:
            assert adata.shape == (2, 2)
            assert adata.obs_names.tolist() == ["1:10", "2:20"]
            assert adata.var_names.tolist() == ["G1", "G2"]
            matrix: Any = adata.X
            assert matrix is not None
            block = matrix[:]
            assert sparse.issparse(block)
            np.testing.assert_allclose(block.toarray(), [[0.0, 0.7], [1.0, 0.5]])
        finally:
            adata.file.close()

    def test_write_h5ad_empty_stream(self, tmp_path):
        path = tmp_path / "empty.h5ad"
        CellExpressionStream([]).write_h5ad(path)
        adata = ad.read_h5ad(path)
        assert adata.shape == (0, 0)

    def test_write_h5ad_preserves_sample_metadata(self, tmp_path):
        stream = CellExpressionStream(
            _stream_frames(),
            sample_metadata=pd.DataFrame(
                {"organ": ["brain"]}, index=pd.Index([1], name="sample_id")
            ),
        )
        path = tmp_path / "cells_meta.h5ad"

        stream.write_h5ad(path)

        adata = ad.read_h5ad(path)
        obs: Any = adata.obs
        assert obs.loc["1:10", "organ"] == "brain"
        assert pd.isna(obs.loc["2:20", "organ"])
        assert adata.obs["sample_id"].dtype == np.int64

    def test_stream_enrich_merges_local_metadata(self):
        stream = CellExpressionStream(
            _stream_frames(),
            sample_metadata=pd.DataFrame(
                {"organ": ["brain"]}, index=pd.Index([1], name="sample_id")
            ),
        )
        adata = stream.enrich().to_anndata()
        obs: Any = adata.obs
        assert obs.loc["1:10", "organ"] == "brain"
        assert pd.isna(obs.loc["2:20", "organ"])

    def test_client_stream_method_returns_wrapper(self, mock_client, monkeypatch):
        monkeypatch.setattr(
            mock_client,
            "_iter_cell_expression_batches",
            lambda **kwargs: iter(_stream_frames()),
        )

        stream = mock_client.stream_cell_expression(
            sample_ids=[1, 2], gene_ids=["G1", "G2"], dataset_id="ds1"
        )

        assert isinstance(stream, CellExpressionStream)
        assert [len(batch) for batch in stream] == [2, 1]


# ---------------------------------------------------------------------------
# Metacell matches
# ---------------------------------------------------------------------------


class TestMetacellMatchResult:
    def test_find_metacells_returns_result_object(self, mock_client, monkeypatch):
        _capture_response(monkeypatch, mock_client, _metacell_response())
        adata = ad.AnnData(
            X=np.array([[1.0, 0.0], [0.0, 2.0]]),
            obs=pd.DataFrame({"cluster": pd.Categorical(["a", "b"])}, index=pd.Index(["c1", "c2"])),
            var=pd.DataFrame(index=pd.Index(["G1", "G2"])),
        )

        result = mock_client.find_metacells(adata, groupby="cluster", dataset_id="ds1")

        assert isinstance(result, MetacellMatchResult)
        frame = result.to_dataframe()
        assert frame.columns.tolist() == [
            "cluster",
            "rank",
            "metacell_id",
            "sample_id",
            "sample_uuid",
            "cell_type",
            "n_cells",
            "similarity",
            "distance",
            "enrichment_score",
        ]
        pd.testing.assert_frame_equal(frame, result.to_pandas())
        assert list(result.columns)[0] == "cluster"

    def test_enrich_adds_local_sample_metadata(self, mock_client, monkeypatch):
        _capture_response(monkeypatch, mock_client, _metacell_response())
        adata = ad.AnnData(
            X=np.ones((2, 2)),
            obs=pd.DataFrame({"cluster": pd.Categorical(["a", "b"])}, index=pd.Index(["c1", "c2"])),
            var=pd.DataFrame(index=pd.Index(["G1", "G2"])),
        )
        metadata = pd.DataFrame({"organ": ["brain"]}, index=pd.Index([1], name="sample_id"))

        result = mock_client.find_metacells(
            adata, groupby="cluster", dataset_id="ds1", sample_metadata=metadata
        )
        result.enrich()

        frame = result.to_dataframe()
        assert set(frame.loc[frame["sample_id"] == 1, "organ"].tolist()) == {"brain"}
        assert pd.isna(frame.loc[frame["sample_id"] == 2, "organ"]).all()

    def test_assign_labels_and_write_obs(self, mock_client, monkeypatch):
        _capture_response(monkeypatch, mock_client, _metacell_response())
        adata = ad.AnnData(
            X=np.ones((2, 2)),
            obs=pd.DataFrame({"cluster": pd.Categorical(["a", "b"])}, index=pd.Index(["c1", "c2"])),
            var=pd.DataFrame(index=pd.Index(["G1", "G2"])),
        )
        result = mock_client.find_metacells(adata, groupby="cluster", dataset_id="ds1")

        labels = result.assign_labels(
            adata, column="malva_cell_type", confidence_column="malva_confidence"
        )

        assert labels.set_index("cluster").loc["a", "cell_type"] == "T"
        assert labels.set_index("cluster").loc["b", "cell_type"] == "B"
        assert labels.set_index("cluster").loc["a", "confidence"] == pytest.approx(1.0)
        assert labels.set_index("cluster").loc["b", "confidence"] == pytest.approx(0.8)
        assert adata.obs["malva_cell_type"].tolist() == ["T", "B"]
        np.testing.assert_allclose(adata.obs["malva_confidence"].tolist(), [1.0, 0.8])

    def test_assign_labels_projects_arbitrary_metadata_column(self, mock_client, monkeypatch):
        response = _metacell_response()
        # Make the sample-3 hit the strongest neighbour of cluster b so the
        # projected disease labels differ between clusters.
        response["clusters"][1]["neighbors"][0]["similarity"] = 0.2
        response["clusters"][1]["neighbors"][0]["distance"] = 0.8
        response["clusters"][1]["neighbors"][1]["similarity"] = 0.8
        response["clusters"][1]["neighbors"][1]["distance"] = 0.2
        _capture_response(monkeypatch, mock_client, response)
        adata = ad.AnnData(
            X=np.ones((2, 2)),
            obs=pd.DataFrame({"cluster": pd.Categorical(["a", "b"])}, index=pd.Index(["c1", "c2"])),
            var=pd.DataFrame(index=pd.Index(["G1", "G2"])),
        )
        metadata = pd.DataFrame(
            {"disease": ["normal", "normal", "cancer"]},
            index=pd.Index([1, 2, 3], name="sample_id"),
        )

        result = mock_client.find_metacells(
            adata, groupby="cluster", dataset_id="ds1", sample_metadata=metadata
        )
        labels = result.assign_labels(adata, column="malva_disease", label_column="disease")

        assert labels.columns.tolist() == [
            "cluster",
            "disease",
            "confidence",
            "n_metacells",
            "weight",
        ]
        by_cluster = labels.set_index("cluster")
        assert by_cluster.loc["a", "disease"] == "normal"
        assert by_cluster.loc["a", "confidence"] == pytest.approx(1.0)
        assert by_cluster.loc["b", "disease"] == "cancer"
        assert by_cluster.loc["b", "confidence"] == pytest.approx(0.8)
        assert result.labels_ is labels
        assert adata.obs["malva_disease"].tolist() == ["normal", "cancer"]

    def test_assign_labels_rejects_missing_label_column(self, mock_client, monkeypatch):
        _capture_response(monkeypatch, mock_client, _metacell_response())
        adata = ad.AnnData(
            X=np.ones((2, 2)),
            obs=pd.DataFrame({"cluster": pd.Categorical(["a", "b"])}, index=pd.Index(["c1", "c2"])),
            var=pd.DataFrame(index=pd.Index(["G1", "G2"])),
        )
        result = mock_client.find_metacells(adata, groupby="cluster", dataset_id="ds1")

        with pytest.raises(ValueError, match="no 'batch' column"):
            result.assign_labels(label_column="batch", write_to_obs=False)

    def test_assign_labels_without_adata(self, mock_client, monkeypatch):
        _capture_response(monkeypatch, mock_client, _metacell_response())
        adata = ad.AnnData(
            X=np.ones((2, 2)),
            obs=pd.DataFrame({"cluster": pd.Categorical(["a", "b"])}, index=pd.Index(["c1", "c2"])),
            var=pd.DataFrame(index=pd.Index(["G1", "G2"])),
        )
        result = mock_client.find_metacells(adata, groupby="cluster", dataset_id="ds1")

        labels = result.assign_labels()

        assert labels.columns.tolist() == [
            "cluster",
            "cell_type",
            "confidence",
            "n_metacells",
            "weight",
        ]
        assert result.labels_ is labels

    def test_deprecated_alias_still_works(self, mock_client, monkeypatch):
        _capture_response(monkeypatch, mock_client, _metacell_response())
        adata = ad.AnnData(
            X=np.ones((2, 2)),
            obs=pd.DataFrame({"cluster": pd.Categorical(["a", "b"])}, index=pd.Index(["c1", "c2"])),
            var=pd.DataFrame(index=pd.Index(["G1", "G2"])),
        )

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            frame = mock_client.compare_clusters_to_metacells(
                adata, groupby="cluster", dataset_id="ds1"
            )

        assert any(issubclass(item.category, DeprecationWarning) for item in caught)
        assert isinstance(frame, MetacellMatchResult)
        assert "cluster" in frame.columns
