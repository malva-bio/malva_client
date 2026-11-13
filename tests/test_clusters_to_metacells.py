"""Offline tests for cluster pseudobulk → Malva metacell comparison."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

BASE = "https://test.malva.bio"


def _adata(X, labels):
    return ad.AnnData(
        X=X,
        obs=pd.DataFrame(
            {"cluster": labels}, index=pd.Index([f"cell{i}" for i in range(len(labels))])
        ),
        var=pd.DataFrame(index=pd.Index(["G1", "G2", "G3"])),
    )


def _response(clusters=None, missing=0, missing_genes=None, scores=None):
    clusters = clusters or ["a", "b"]
    return {
        "dataset_id": "ds1",
        "metric": "cosine",
        "top_n": 50,
        "n_metacells": 100,
        "n_clusters": len(clusters),
        "clusters": [
            {
                "name": name,
                "n_query_genes": 3,
                "n_matched_genes": 3 - missing,
                "n_missing_genes": missing,
                "missing_genes": missing_genes or [],
                "cell_type_scores": (scores or {}).get(name, []),
                "neighbors": [
                    {
                        "metacell_id": f"1__{name}",
                        "sample_id": 1,
                        "sample_uuid": f"uuid-{name}",
                        "cell_type": name,
                        "n_cells": 10,
                        "similarity": 0.9,
                        "distance": 0.1,
                        "enrichment_score": 1.5,
                        "rank": 1,
                    }
                ],
            }
            for name in clusters
        ],
        "query_time_ms": 2.0,
    }


def _capture(monkeypatch, client, response):
    captured = {}

    def fake(method, endpoints, **kwargs):
        captured["method"] = method
        captured["endpoints"] = endpoints
        captured["json"] = kwargs["json"]
        return response, endpoints[0]

    monkeypatch.setattr(client, "_request_first_available", fake)
    return captured


class TestAggregation:
    def test_sparse_cluster_means_are_sent_as_sparse_profiles(self, mock_client, monkeypatch):
        X = sparse.csr_matrix(np.array([[1.5, 0.0, 2.25], [1.5, 0.0, 0.0], [0.0, 3.5, 0.0]]))
        adata = _adata(X, pd.Categorical(["a", "a", "b"]))
        captured = _capture(monkeypatch, mock_client, _response())

        frame = mock_client.compare_clusters_to_metacells(
            adata, "cluster", top_n=5, dataset_id="ds1"
        )

        payload = captured["json"]
        assert payload["dataset_id"] == "ds1"
        assert payload["top_n"] == 5
        assert payload["metric"] == "cosine"
        assert "space" not in payload
        clusters = {c["name"]: c for c in payload["clusters"]}
        assert set(clusters) == {"a", "b"}
        assert clusters["a"]["gene_ids"] == ["G1", "G3"]
        assert clusters["a"]["values"] == pytest.approx([1.5, 1.125])
        assert clusters["b"]["gene_ids"] == ["G2"]
        assert clusters["b"]["values"] == pytest.approx([3.5])
        # only a few cluster × gene values, not all cells
        assert sum(len(c["values"]) for c in payload["clusters"]) == 3

        assert list(frame.columns) == [
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
        assert frame["cluster"].tolist() == ["a", "b"]
        assert frame["rank"].tolist() == [1, 1]
        assert frame["sample_uuid"].tolist() == ["uuid-a", "uuid-b"]

    def test_find_metacells_prefers_canonical_endpoint(self, mock_client, monkeypatch):
        adata = _adata(sparse.csr_matrix(np.eye(3)), pd.Categorical(["a", "b", "c"]))
        captured = _capture(monkeypatch, mock_client, _response(["a", "b", "c"]))

        mock_client.find_metacells(adata, "cluster", dataset_id="ds1")

        assert captured["endpoints"] == [
            "/api/coexpression/metacells/nearest",
            "/query/metacells/nearest",
            "/api/coexpression/clusters-to-metacells",
            "/query/clusters-to-metacells",
        ]

    def test_dense_input_is_supported(self, mock_client, monkeypatch):
        adata = _adata(np.array([[1.5, 0.0, 2.25], [1.5, 0.0, 0.0]]), pd.Categorical(["a", "a"]))
        captured = _capture(monkeypatch, mock_client, _response(clusters=["a"]))

        mock_client.compare_clusters_to_metacells(adata, "cluster", dataset_id="ds1")

        assert captured["json"]["clusters"][0]["gene_ids"] == ["G1", "G3"]

    def test_empty_cluster_is_skipped_with_warning(self, mock_client, monkeypatch):
        X = sparse.csr_matrix(np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]]))
        adata = _adata(X, pd.Categorical(["empty", "a"]))
        captured = _capture(monkeypatch, mock_client, _response(clusters=["a"]))

        with pytest.warns(UserWarning, match="Skipping clusters with no non-zero"):
            mock_client.compare_clusters_to_metacells(adata, "cluster", dataset_id="ds1")

        assert [c["name"] for c in captured["json"]["clusters"]] == ["a"]

    def test_all_empty_clusters_raises(self, mock_client, monkeypatch):
        adata = _adata(np.zeros((2, 3)), pd.Categorical(["a", "a"]))
        _capture(monkeypatch, mock_client, _response())

        with pytest.warns(UserWarning):
            with pytest.raises(ValueError, match="No clusters with non-zero expression"):
                mock_client.compare_clusters_to_metacells(adata, "cluster", dataset_id="ds1")


class TestWarnings:
    def test_non_categorical_groupby_warns(self, mock_client, monkeypatch):
        adata = _adata(np.array([[1.5, 0.0, 0.0], [1.5, 0.0, 0.0]]), ["a", "a"])
        _capture(monkeypatch, mock_client, _response(clusters=["a"]))

        with pytest.warns(UserWarning, match="is not categorical"):
            mock_client.compare_clusters_to_metacells(adata, "cluster", dataset_id="ds1")

    def test_raw_count_x_warns_and_data_is_unchanged(self, mock_client, monkeypatch):
        X = np.array([[1.0, 0.0, 2.0], [4.0, 0.0, 0.0]])
        adata = _adata(X, pd.Categorical(["a", "a"]))
        before = adata.X.copy()
        _capture(monkeypatch, mock_client, _response(clusters=["a"]))

        with pytest.warns(UserWarning, match="looks like raw counts"):
            mock_client.compare_clusters_to_metacells(adata, "cluster", dataset_id="ds1")

        assert np.array_equal(adata.X, before)

    def test_missing_atlas_genes_warn(self, mock_client, monkeypatch):
        adata = _adata(np.array([[1.5, 0.0, 2.25], [1.5, 0.0, 0.0]]), pd.Categorical(["a", "a"]))
        _capture(
            monkeypatch,
            mock_client,
            _response(clusters=["a"], missing=2, missing_genes=["X", "Y"]),
        )

        with pytest.warns(UserWarning, match="2 gene.*not found"):
            mock_client.compare_clusters_to_metacells(adata, "cluster", dataset_id="ds1")


class TestLocalMetadata:
    def test_sample_metadata_is_merged_locally(self, mock_client, monkeypatch):
        adata = _adata(np.array([[1.5, 0.0, 2.25], [1.5, 0.0, 0.0]]), pd.Categorical(["a", "a"]))
        _capture(monkeypatch, mock_client, _response(clusters=["a"]))
        metadata = pd.DataFrame(
            {"tissue": ["brain"], "cell_type": ["should_not_duplicate"]},
            index=pd.Index([1], name="sample_id"),
        )

        frame = mock_client.compare_clusters_to_metacells(
            adata, "cluster", dataset_id="ds1", sample_metadata=metadata
        )

        assert frame["tissue"].tolist() == ["brain"]
        assert frame["cell_type"].tolist() == ["a"]
        assert list(frame.columns).count("cell_type") == 1

    def test_dict_metadata_is_merged(self, mock_client, monkeypatch):
        adata = _adata(np.array([[1.5, 0.0, 2.25], [1.5, 0.0, 0.0]]), pd.Categorical(["a", "a"]))
        _capture(monkeypatch, mock_client, _response(clusters=["a"]))

        frame = mock_client.compare_clusters_to_metacells(
            adata,
            "cluster",
            dataset_id="ds1",
            sample_metadata={1: {"organ": "brain"}},
        )

        assert frame["organ"].tolist() == ["brain"]


class TestValidation:
    def test_top_n_is_capped(self, mock_client):
        adata = _adata(np.full((2, 3), 1.5), pd.Categorical(["a", "a"]))
        with pytest.raises(ValueError, match="top_n must be between 1 and 1000"):
            mock_client.compare_clusters_to_metacells(
                adata, "cluster", top_n=1001, dataset_id="ds1"
            )

    def test_dataset_id_is_required(self, mock_client):
        adata = _adata(np.full((2, 3), 1.5), pd.Categorical(["a", "a"]))
        with pytest.raises(ValueError, match="dataset_id is required"):
            mock_client.compare_clusters_to_metacells(adata, "cluster")

    def test_unknown_groupby_raises(self, mock_client):
        adata = _adata(np.full((2, 3), 1.5), pd.Categorical(["a", "a"]))
        with pytest.raises(KeyError, match="not a column"):
            mock_client.compare_clusters_to_metacells(adata, "nope", dataset_id="ds1")

    def test_feature_ids_must_align(self, mock_client):
        adata = _adata(np.full((2, 3), 1.5), pd.Categorical(["a", "a"]))
        with pytest.raises(ValueError, match="feature_ids must have one entry"):
            mock_client.compare_clusters_to_metacells(
                adata, "cluster", dataset_id="ds1", feature_ids=["G1"]
            )

    def test_hellinger_method_is_forwarded(self, mock_client, monkeypatch):
        adata = _adata(np.array([[1.5, 0.0, 2.5], [1.5, 0.0, 0.0]]), pd.Categorical(["a", "a"]))
        captured = _capture(monkeypatch, mock_client, _response(clusters=["a"]))
        mock_client.find_metacells(adata, "cluster", dataset_id="ds1", method="hellinger")
        payload = captured["json"]
        assert payload["space"] == "hellinger"
        assert payload["clusters"][0]["transform"] == "auto"

    def test_panel_lda_method_is_forwarded(self, mock_client, monkeypatch):
        adata = _adata(np.array([[1.0, 0.0, 2.0]]), pd.Categorical(["a"]))
        captured = _capture(monkeypatch, mock_client, _response(clusters=["a"]))
        mock_client.find_metacells(adata, "cluster", dataset_id="ds1", method="panel_lda")
        assert captured["json"]["space"] == "panel_lda"

    def test_invalid_method_raises(self, mock_client):
        adata = _adata(np.array([[1.0, 0.0, 2.0]]), pd.Categorical(["a"]))
        with pytest.raises(ValueError, match="method must be"):
            mock_client.find_metacells(adata, "cluster", dataset_id="ds1", method="cosine2")

    def test_server_scores_drive_assign_labels_by_default(self, mock_client, monkeypatch):
        adata = _adata(np.array([[1.0, 0.0, 2.0], [1.0, 0.0, 0.0]]), pd.Categorical(["a", "a"]))
        scores = {
            "a": [
                {"cell_type": "T", "score": 0.2, "n_metacells": 2, "rank": 2},
                {"cell_type": "B", "score": 0.8, "n_metacells": 2, "rank": 1},
            ]
        }
        _capture(monkeypatch, mock_client, _response(clusters=["a"], scores=scores))
        result = mock_client.find_metacells(adata, "cluster", dataset_id="ds1", method="hellinger")
        labels = result.assign_labels(write_to_obs=False)
        assert labels.iloc[0]["cell_type"] == "B"
        assert labels.iloc[0]["confidence"] == pytest.approx(0.8 / 1.0)

    def test_explicit_similarity_weight_keeps_legacy_aggregation(self, mock_client, monkeypatch):
        adata = _adata(np.array([[1.0, 0.0, 2.0], [1.0, 0.0, 0.0]]), pd.Categorical(["a", "a"]))
        scores = {
            "a": [
                {"cell_type": "T", "score": 0.2, "n_metacells": 2, "rank": 2},
                {"cell_type": "B", "score": 0.8, "n_metacells": 2, "rank": 1},
            ]
        }
        _capture(monkeypatch, mock_client, _response(clusters=["a"], scores=scores))
        result = mock_client.find_metacells(adata, "cluster", dataset_id="ds1", method="hellinger")
        labels = result.assign_labels(weight="similarity", write_to_obs=False)
        assert labels.iloc[0]["cell_type"] == "a"
