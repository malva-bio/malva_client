"""Opt-in Leiden grouping: request payloads and SearchResult metadata."""

import pytest

from malva_client.client import MalvaClient
from malva_client.models import SearchResult


def _client(monkeypatch, responses=None):
    client = MalvaClient.__new__(MalvaClient)
    client.base_url = "http://explorer.example"
    client.dataset_id = None
    client.provider_id = None
    client._endpoint_style = "legacy"
    calls = []

    def fake_request(method, endpoint, **kwargs):
        calls.append((method, endpoint, kwargs))
        if responses is not None:
            return responses(endpoint, kwargs)
        return {"job_id": "job-1", "status": "completed"}

    monkeypatch.setattr(client, "_request", fake_request)
    return client, calls


def _legacy_results() -> dict:
    return {
        "FOXP3": {
            "expression_data": {
                "samples": ["s1"],
                "cell_types": ["T cell"],
                "columns": [
                    "sample_idx",
                    "cell_type_idx",
                    "norm_expr",
                    "raw_expr",
                    "cell_count",
                ],
                "data": [[0, 0, 0.9, 1.1, 10]],
            }
        }
    }


def test_search_without_group_by_keeps_the_legacy_payload(monkeypatch):
    client, calls = _client(monkeypatch)
    client.search("FOXP3")

    payload = calls[0][2]["json"]
    assert payload == {
        "query": "FOXP3",
        "wait_for_completion": False,
        "aggregate_expression": True,
        "min_kmer_presence": 0,
        "max_kmer_presence": 100000,
    }


def test_search_passes_leiden_only_when_requested(monkeypatch):
    client, calls = _client(monkeypatch)
    client.search("FOXP3", group_by="leiden")

    assert calls[0][2]["json"]["group_by"] == "leiden"


def test_search_genes_and_federated_forward_group_by(monkeypatch):
    client, calls = _client(monkeypatch)
    client.search_genes(["FOXP3"], group_by="leiden")
    assert calls[0][2]["json"]["group_by"] == "leiden"

    def responses(endpoint, kwargs):
        if endpoint == "/api/workspace":
            return {"provider_ids": ["p1"]}
        if endpoint == "/api/expression/submit":
            assert kwargs["json"]["group_by"] == "leiden"
            return {
                "job_id": "fed-1",
                "status": "completed",
                "shards": [
                    {
                        "provider_id": "p1",
                        "provider_name": "P1",
                        "search_job_id": "shard-1",
                        "status": "completed",
                        "error": None,
                    }
                ],
            }
        if endpoint == "/api/expression/results/shard-1":
            return {
                "status": "completed",
                "genes": _legacy_results(),
                "grouping": {
                    "requested": "leiden",
                    "applied": "leiden",
                    "axes": ["cell_type", "leiden"],
                    "annotation": {"version": "v1", "method": "marker"},
                    "fallback_reason": None,
                },
                "cluster_annotations": {"s1": {"0": {"cell_type": "T cell", "n_cells": 10}}},
            }
        raise AssertionError(endpoint)

    client, _ = _client(monkeypatch, responses)
    frame = client.search_federated("FOXP3", provider_ids=["p1"], group_by="leiden")
    assert frame.attrs["grouping"]["applied"] == "leiden"
    assert frame.attrs["cluster_annotations"]["p1"]["s1"]["0"]["cell_type"] == "T cell"


def test_invalid_group_by_is_rejected_before_submission(monkeypatch):
    client, calls = _client(monkeypatch)
    with pytest.raises(ValueError, match="group_by"):
        client.search("FOXP3", group_by="clusters")
    assert calls == []


def test_search_result_exposes_grouping_and_cluster_annotations(monkeypatch):
    client, _ = _client(monkeypatch)
    raw = {
        "job_id": "job-1",
        "status": "completed",
        "results": {
            "FOXP3": {
                "expression_data": {
                    "samples": ["s1"],
                    "clusters": ["0", "1"],
                    "columns": [
                        "sample_idx",
                        "cell_type_idx",
                        "norm_expr",
                        "raw_expr",
                        "cell_count",
                    ],
                    "data": [[0, 0, 0.9, 1.1, 10], [0, 1, 0.2, 0.3, 5]],
                }
            }
        },
        "grouping": {
            "requested": "leiden",
            "applied": "leiden",
            "axes": ["cell_type", "leiden"],
            "annotation": {"version": "v1", "method": "marker"},
            "fallback_reason": None,
        },
        "cluster_annotations": {
            "s1": {
                "0": {
                    "cell_type": "T cell",
                    "cl_id": "cl-0",
                    "annotation_method": "marker",
                    "annotation_version": "v1",
                    "n_cells": 10,
                },
                "1": {
                    "cell_type": "B cell",
                    "cl_id": "cl-1",
                    "annotation_method": "marker",
                    "annotation_version": "v1",
                    "n_cells": 5,
                },
            }
        },
    }
    result = SearchResult(raw, client)

    assert result.grouping["requested"] == "leiden"
    assert result.grouping["applied"] == "leiden"
    assert result.grouping["annotation"]["version"] == "v1"

    annotations = result.cluster_annotations
    assert list(annotations.columns) == [
        "sample_id",
        "cluster_id",
        "cell_type",
        "cl_id",
        "annotation_method",
        "annotation_version",
        "n_cells",
    ]
    assert set(annotations["cell_type"]) == {"T cell", "B cell"}

    frame = result.df
    assert "cluster_id" in frame.columns
    assert frame.loc[frame["cluster_id"] == "0", "cell_type"].iloc[0] == "T cell"
    assert frame.loc[frame["cluster_id"] == "1", "cell_type"].iloc[0] == "B cell"


def test_legacy_result_reports_cell_type_and_an_empty_annotation_table(monkeypatch):
    client, _ = _client(monkeypatch)
    result = SearchResult(
        {"job_id": "job-1", "status": "completed", "results": _legacy_results()}, client
    )

    assert result.grouping["applied"] == "cell_type"
    assert result.grouping["requested"] == "cell_type"
    annotations = result.cluster_annotations
    assert annotations.empty
    assert "cell_type" in annotations.columns
    assert "cluster_id" not in result.df.columns
