"""Canonical endpoint routing and automatic legacy fallback."""

import json
from unittest.mock import patch
from urllib.parse import parse_qs, urlparse

import msgpack
import numpy as np
import pandas as pd
import pytest
import responses

from malva_client.client import MalvaClient
from malva_client.exceptions import AuthenticationError, MalvaAPIError
from malva_client.models import SearchResult

BASE = "https://test.malva.bio"

CANONICAL_TEMPLATES = {
    "health": "/api/health",
    "quota": "/api/account/quota",
    "datasets": "/api/datasets",
    "samples": "/api/samples",
    "sample_details": "/api/samples/sample-1",
    "sample_metadata": "/api/samples/metadata",
    "studies": "/api/metadata/studies",
    "overview": "/api/metadata/overview",
    "available_filters": "/api/metadata/available-filters",
    "filter_values": "/api/metadata/filters/organ/values",
    "datasets_hierarchy": "/api/metadata/datasets/hierarchy",
    "dataset_studies": "/api/metadata/datasets/ds-1/studies",
    "study_samples": "/api/metadata/datasets/ds-1/studies/study-a/samples",
    "search": "/api/search",
    "search_job": "/api/search/job-1",
    "search_job_status": "/api/search/job-1/status",
    "search_job_results": "/api/search/job-1/results",
    "search_job_cells": "/api/search/job-1/cells",
    "search_job_cell_ids": "/api/search/job-1/cell-ids",
    "search_job_global_cell_ids": "/api/search/job-1/global-cell-ids",
    "search_federated": "/api/search",
    "search_federated_job": "/api/search/federated/job-1",
    "search_shard_results": "/api/search/job-1/results",
    "cells_by_metadata": "/api/cells/by-metadata",
    "cell_expression_export": "/api/exports/cell-expression",
    "cell_expression_export_job": "/api/exports/cell-expression/export-1",
    "cell_expression_export_download": "/api/exports/cell-expression/export-1/download",
    "all_cells_export": "/api/exports/all-cells",
    "all_cells_export_job": "/api/exports/all-cells/export-1",
    "all_cells_export_download": "/api/exports/all-cells/export-1/download",
    "download_token": "/api/downloads/token-1",
}

LEGACY_TEMPLATES = {
    "health": "/health",
    "quota": "/api/quota-status",
    "datasets": "/datasets",
    "samples": "/samples/api/search",
    "sample_details": "/samples/api/sample/sample-1",
    "sample_metadata": "/api/sample-metadata",
    "studies": "/samples/api/studies",
    "overview": "/samples/api/overview",
    "available_filters": "/samples/api/filters",
    "filter_values": "/samples/api/filters/organ/values",
    "datasets_hierarchy": "/samples/api/datasets/hierarchy",
    "dataset_studies": "/samples/api/datasets/ds-1/studies",
    "study_samples": "/samples/api/datasets/ds-1/studies/study-a/samples",
    "search": "/search",
    "search_job": "/search/job-1",
    "search_job_status": "/search/job-1/status",
    "search_job_results": "/api/expression-data/job-1/results",
    "search_job_cells": "/api/expression/results/job-1/cells",
    "search_job_cell_ids": "/search/job-1/cell-ids",
    "search_job_global_cell_ids": "/search/job-1/global-cell-ids",
    "search_federated": "/api/expression/submit",
    "search_federated_job": "/api/expression/federated/job-1",
    "search_shard_results": "/api/expression/results/job-1",
    "cells_by_metadata": "/api/cells/by-filter",
    "cell_expression_export": "/exports/cell-expression",
    "cell_expression_export_job": "/exports/cell-expression/export-1",
    "cell_expression_export_download": "/exports/cell-expression/export-1/download",
    "all_cells_export": "/exports/all-cells",
    "all_cells_export_job": "/exports/all-cells/export-1",
    "all_cells_export_download": "/exports/all-cells/export-1/download",
    "download_token": "/api/download/token-1",
}

TEMPLATE_VALUES = {
    "job_id": "job-1",
    "export_id": "export-1",
    "sample_uuid": "sample-1",
    "column": "organ",
    "dataset_id": "ds-1",
    "study_name": "study-a",
    "token": "token-1",
}


def _client(monkeypatch, style):
    if style is None:
        monkeypatch.delenv("MALVA_ENDPOINT_STYLE", raising=False)
    else:
        monkeypatch.setenv("MALVA_ENDPOINT_STYLE", style)
    with patch.object(MalvaClient, "_test_connection"):
        return MalvaClient(base_url=BASE, api_token="tok")


def _assert_request(request, method, path, expected_query=None):
    parsed = urlparse(request.url)
    assert request.method == method
    assert parsed.path == path
    if expected_query is not None:
        assert parse_qs(parsed.query) == expected_query


def _json_callback(method, path, payload, expected_query=None):
    def callback(request):
        _assert_request(request, method, path, expected_query)
        return 200, {"content-type": "application/json"}, json.dumps(payload)

    return callback


def _msgpack_callback(method, path, payload, expected_query=None):
    def callback(request):
        _assert_request(request, method, path, expected_query)
        return 200, {"content-type": "application/msgpack"}, msgpack.packb(payload)

    return callback


def _gene_payload(sample_id=1):
    return {
        "GENE": {
            "expression_data": {
                "samples": [sample_id],
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


@pytest.mark.parametrize("operation,path", sorted(CANONICAL_TEMPLATES.items()))
def test_canonical_endpoint_template(monkeypatch, operation, path):
    client = _client(monkeypatch, "canonical")
    assert client._endpoint(operation, **TEMPLATE_VALUES) == path


@pytest.mark.parametrize("operation,path", sorted(LEGACY_TEMPLATES.items()))
def test_legacy_endpoint_template(monkeypatch, operation, path):
    client = _client(monkeypatch, "legacy")
    assert client._endpoint(operation, **TEMPLATE_VALUES) == path


@responses.activate
def test_canonical_dataset_and_metadata_requests(monkeypatch):
    client = _client(monkeypatch, "canonical")

    responses.add_callback(
        responses.GET,
        f"{BASE}/api/datasets",
        _json_callback(
            "GET",
            "/api/datasets",
            {"datasets": [{"dataset_id": "ds-1"}]},
            {"provider": ["p1"]},
        ),
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/api/metadata/studies",
        _json_callback(
            "GET",
            "/api/metadata/studies",
            {"studies": [{"study_name": "study-a"}]},
            {"provider": ["p1"]},
        ),
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/api/metadata/overview",
        _json_callback("GET", "/api/metadata/overview", {"total_samples": 1}, {"provider": ["p1"]}),
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/api/metadata/available-filters",
        _json_callback(
            "GET",
            "/api/metadata/available-filters",
            {"filters": {}},
            {"provider": ["p1"]},
        ),
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/api/metadata/datasets/hierarchy",
        _json_callback(
            "GET",
            "/api/metadata/datasets/hierarchy",
            {"datasets": []},
            {"provider": ["p1"]},
        ),
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/api/metadata/datasets/ds-1/studies",
        _json_callback(
            "GET",
            "/api/metadata/datasets/ds-1/studies",
            {"studies": [{"study_name": "study-a"}]},
            {"provider": ["p1"], "page": ["1"], "page_size": ["50"]},
        ),
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/api/metadata/datasets/ds-1/studies/study-a/samples",
        _json_callback(
            "GET",
            "/api/metadata/datasets/ds-1/studies/study-a/samples",
            {"samples": [{"sample_id": 1}]},
            {"provider": ["p1"], "page": ["1"], "page_size": ["50"]},
        ),
    )

    assert client.get_datasets(provider="p1") == [{"dataset_id": "ds-1"}]
    assert client.get_studies(provider="p1").to_dict("records") == [{"study_name": "study-a"}]
    assert client.get_overview_stats(provider="p1")["total_samples"] == 1
    assert client.get_available_filters(provider="p1") == {"filters": {}}
    assert client.get_datasets_hierarchy(provider="p1") == {"datasets": []}
    assert client.get_dataset_studies("ds-1", provider="p1")["studies"]
    assert client.get_study_samples("ds-1", "study-a", provider="p1")["samples"]


@responses.activate
def test_canonical_sample_requests(monkeypatch):
    client = _client(monkeypatch, "canonical")

    responses.add_callback(
        responses.GET,
        f"{BASE}/api/samples",
        _json_callback(
            "GET",
            "/api/samples",
            {"samples": [{"sample_id": 1}], "has_more": False},
            {
                "page": ["2"],
                "page_size": ["10"],
                "q": ["brain"],
                "provider": ["p1"],
                "organ": ["brain"],
            },
        ),
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/api/samples/sample-1",
        _json_callback(
            "GET",
            "/api/samples/sample-1",
            {"sample_uuid": "sample-1"},
            {"provider": ["p1"]},
        ),
    )
    responses.add_callback(
        responses.POST,
        f"{BASE}/api/samples/metadata",
        _json_callback(
            "POST",
            "/api/samples/metadata",
            {"metadata": {"1": {"organ": "brain"}}},
            {"provider": ["p1"]},
        ),
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/api/metadata/filters/organ/values",
        _json_callback(
            "GET",
            "/api/metadata/filters/organ/values",
            {"values": ["brain"]},
            {"provider": ["p1"], "limit": ["50"], "search": ["bra"]},
        ),
    )

    response = client.get_samples(
        page=2,
        page_size=10,
        filters={"organ": ["brain"]},
        search_query="brain",
        provider="p1",
    )
    assert response["samples"] == [{"sample_id": 1}]
    assert client.get_sample_details("sample-1", provider="p1")["sample_uuid"] == "sample-1"
    assert client.get_sample_metadata([1], provider="p1") == {1: {"organ": "brain"}}
    assert client.get_filter_values("organ", search="bra", provider="p1") == ["brain"]


@responses.activate
def test_canonical_search_and_job_requests(monkeypatch):
    client = _client(monkeypatch, "canonical")

    responses.add_callback(
        responses.POST,
        f"{BASE}/api/search",
        _json_callback("POST", "/api/search", {"job_id": "job-1", "status": "completed"}),
    )
    responses.add_callback(
        responses.POST,
        f"{BASE}/api/search",
        _json_callback("POST", "/api/search", {"job_id": "job-2", "status": "completed"}),
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/api/search/job-1/status",
        _json_callback("GET", "/api/search/job-1/status", {"status": "completed"}),
    )
    responses.add_callback(
        responses.DELETE,
        f"{BASE}/api/search/job-1",
        _json_callback("DELETE", "/api/search/job-1", {"success": True}),
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/api/search/job-1/results",
        _json_callback("GET", "/api/search/job-1/results", {"results": _gene_payload()}),
    )

    assert client.search("BRCA1").job_id == "job-1"
    assert client.search_genes(["BRCA1"]).job_id == "job-2"
    assert client.get_job_status("job-1")["status"] == "completed"
    assert client.cancel_job("job-1") is True

    result = SearchResult({"job_id": "job-1"}, client)
    assert "GENE" in result.df["gene_sequence"].tolist()


@responses.activate
def test_canonical_federated_requests(monkeypatch):
    client = _client(monkeypatch, "canonical")
    shards = [
        {
            "provider_id": "p1",
            "provider_name": "Provider One",
            "search_job_id": "shard-1",
            "status": "completed",
            "error": None,
        }
    ]

    responses.add_callback(
        responses.GET,
        f"{BASE}/api/workspace",
        _json_callback("GET", "/api/workspace", {"provider_ids": ["p1"]}),
    )

    def federated_submit(request):
        _assert_request(request, "POST", "/api/search")
        body = json.loads(request.body)
        assert body["provider_ids"] == ["p1"]
        assert body["query"] == "FOXP3"
        return (
            200,
            {"content-type": "application/json"},
            json.dumps({"job_id": "fed-1", "status": "running", "shards": shards}),
        )

    responses.add_callback(responses.POST, f"{BASE}/api/search", federated_submit)
    responses.add_callback(
        responses.GET,
        f"{BASE}/api/search/federated/fed-1",
        _json_callback(
            "GET",
            "/api/search/federated/fed-1",
            {"job_id": "fed-1", "status": "completed", "shards": shards},
        ),
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/api/search/shard-1/results",
        _json_callback(
            "GET",
            "/api/search/shard-1/results",
            {"status": "completed", "genes": _gene_payload()},
            {"provider": ["p1"]},
        ),
    )

    frame = client.search_federated("FOXP3", provider_ids=["p1"], poll_interval=0, max_wait=5)
    assert frame["provider_id"].tolist() == ["p1"]
    assert frame.attrs["job_id"] == "fed-1"


@responses.activate
def test_canonical_cells_by_metadata_request(monkeypatch):
    client = _client(monkeypatch, "canonical")

    responses.add_callback(
        responses.POST,
        f"{BASE}/api/cells/by-metadata",
        _msgpack_callback(
            "POST",
            "/api/cells/by-metadata",
            {
                "s": np.array([1], dtype=np.int64).tobytes(),
                "c": np.array([10], dtype=np.int64).tobytes(),
            },
            {"provider": ["p1"]},
        ),
    )

    cells = client.get_cells_by_metadata(
        sample_ids=[1],
        include_cell_metadata=False,
        provider="p1",
    )
    assert cells.to_dict("records") == [{"sample_id": 1, "cell_id": 10}]


@responses.activate
def test_canonical_all_cells_export_requests(monkeypatch):
    client = _client(monkeypatch, "canonical")

    responses.add_callback(
        responses.POST,
        f"{BASE}/api/exports/all-cells",
        _json_callback(
            "POST",
            "/api/exports/all-cells",
            {"export_id": "export-1", "status": "pending"},
        ),
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/api/exports/all-cells/export-1",
        _json_callback(
            "GET",
            "/api/exports/all-cells/export-1",
            {"export_id": "export-1", "status": "completed"},
        ),
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/api/exports/all-cells/export-1/download",
        _msgpack_callback(
            "GET",
            "/api/exports/all-cells/export-1/download",
            {
                "s": np.array([1], dtype=np.int64).tobytes(),
                "c": np.array([10], dtype=np.int64).tobytes(),
            },
        ),
    )

    cells = client.get_cells_by_metadata(
        include_all_database_cells=True,
        include_cell_metadata=False,
        poll_interval=0,
    )
    assert cells.to_dict("records") == [{"sample_id": 1, "cell_id": 10}]


@responses.activate
def test_canonical_retrieve_cells_fetchers(monkeypatch):
    client = _client(monkeypatch, "canonical")
    job_id = "job-1"

    responses.add_callback(
        responses.POST,
        f"{BASE}/api/exports/cell-expression",
        _json_callback(
            "POST",
            "/api/exports/cell-expression",
            {"export_id": "export-1", "status": "pending"},
        ),
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/api/exports/cell-expression/export-1",
        _json_callback(
            "GET",
            "/api/exports/cell-expression/export-1",
            {"export_id": "export-1", "status": "completed"},
        ),
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/api/exports/cell-expression/export-1/download",
        _msgpack_callback(
            "GET",
            "/api/exports/cell-expression/export-1/download",
            {
                "format": "cell_expression_values_v1",
                "features": [
                    {
                        "feature_index": 1,
                        "job_id": job_id,
                        "feature": "GENE",
                        "label": "GENE",
                    }
                ],
                "s": np.array([1], dtype=np.int64).tobytes(),
                "c": np.array([10], dtype=np.int64).tobytes(),
                "fi": np.array([1], dtype=np.int32).tobytes(),
                "v": np.array([2.5], dtype=np.float32).tobytes(),
                "n": 1,
            },
        ),
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/api/search/{job_id}",
        _json_callback(
            "GET",
            f"/api/search/{job_id}",
            {
                "status": "completed",
                "results": {"GENE": {"cell": [10], "sample": [1], "expression": [2.5]}},
            },
        ),
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/api/search/{job_id}/cell-ids",
        _json_callback(
            "GET",
            f"/api/search/{job_id}/cell-ids",
            {
                "results": {
                    "GENE": {"cell": [10], "sample": [1], "expression": [2.5]},
                }
            },
        ),
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/api/search/{job_id}/global-cell-ids",
        _msgpack_callback(
            "GET",
            f"/api/search/{job_id}/global-cell-ids",
            {
                "s": np.array([1], dtype=np.int64).tobytes(),
                "c": np.array([10], dtype=np.int64).tobytes(),
            },
        ),
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/api/search/{job_id}/cells",
        _json_callback(
            "GET",
            f"/api/search/{job_id}/cells",
            {"samples": {"1": {"cell_ids": [10]}}},
        ),
    )

    values = client._retrieve_cells_from_values_export_endpoint(
        job_id,
        ["GENE"],
        None,
        False,
        False,
        poll_interval=0,
    )
    assert values.to_dataframe()["value"].tolist() == [2.5]

    full = client._retrieve_cells_from_full_search_endpoint(job_id, ["GENE"], None, False, False)
    assert full.to_dataframe()["value"].tolist() == [2.5]

    cell_ids = client._retrieve_cells_from_cell_ids_endpoint(job_id, ["GENE"], None, False, False)
    assert cell_ids.to_dataframe()["value"].tolist() == [2.5]

    global_ids = client._retrieve_cells_from_global_ids_endpoint(
        job_id, ["GENE"], None, False, False
    )
    assert global_ids.to_dataframe()["value"].tolist() == [1.0]

    explorer = client._retrieve_cells_from_explorer_endpoint(job_id, ["GENE"], None, False, False)
    assert explorer.to_dataframe()["value"].tolist() == [1.0]


@responses.activate
def test_legacy_requests_still_use_historical_paths(monkeypatch):
    client = _client(monkeypatch, "legacy")

    responses.add_callback(
        responses.GET,
        f"{BASE}/samples/api/search",
        _json_callback(
            "GET",
            "/samples/api/search",
            {"samples": [], "has_more": False},
            {"page": ["1"], "page_size": ["50"]},
        ),
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/search/job-1/cell-ids",
        _json_callback(
            "GET",
            "/search/job-1/cell-ids",
            {"results": {"GENE": {"cell": [10], "sample": [1], "expression": [2.0]}}},
        ),
    )

    assert client.get_samples()["samples"] == []
    cells = client._retrieve_cells_from_cell_ids_endpoint("job-1", ["GENE"], None, False, False)
    assert cells.to_dataframe()["value"].tolist() == [2.0]


@responses.activate
def test_auto_detection_canonical_is_cached(monkeypatch):
    client = _client(monkeypatch, None)

    responses.add(
        responses.GET,
        f"{BASE}/api/account/quota",
        json={"account_type": "pro"},
    )
    responses.add(
        responses.GET,
        f"{BASE}/api/samples",
        json={"samples": [], "has_more": False},
    )
    responses.add(
        responses.GET,
        f"{BASE}/api/samples",
        json={"samples": [], "has_more": False},
    )

    assert client.endpoint_style == "canonical"
    client.get_samples()
    client.get_samples()

    probe_calls = [
        call for call in responses.calls if urlparse(call.request.url).path == "/api/account/quota"
    ]
    assert len(probe_calls) == 1
    assert all(urlparse(call.request.url).path == "/api/samples" for call in responses.calls[1:])

    client.reset_endpoint_style()
    assert client.endpoint_style == "canonical"
    probe_calls = [
        call for call in responses.calls if urlparse(call.request.url).path == "/api/account/quota"
    ]
    assert len(probe_calls) == 2


@responses.activate
def test_auto_detection_falls_back_to_legacy(monkeypatch):
    client = _client(monkeypatch, None)

    responses.add(
        responses.GET,
        f"{BASE}/api/account/quota",
        json={"detail": "not found"},
        status=404,
    )
    responses.add(
        responses.GET,
        f"{BASE}/samples/api/search",
        json={"samples": [], "has_more": False},
    )
    responses.add(
        responses.GET,
        f"{BASE}/samples/api/search",
        json={"samples": [], "has_more": False},
    )

    assert client.endpoint_style == "legacy"
    client.get_samples()
    client.get_samples()

    probe_calls = [
        call for call in responses.calls if urlparse(call.request.url).path == "/api/account/quota"
    ]
    assert len(probe_calls) == 1
    assert all(
        urlparse(call.request.url).path == "/samples/api/search" for call in responses.calls[1:]
    )


@responses.activate
def test_detection_401_raises_authentication_error(monkeypatch):
    client = _client(monkeypatch, None)
    responses.add(
        responses.GET,
        f"{BASE}/api/account/quota",
        json={"detail": "unauthorized"},
        status=401,
    )

    with pytest.raises(AuthenticationError):
        client.endpoint_style


def test_explicit_style_skips_probe(monkeypatch):
    client = _client(monkeypatch, "canonical")
    assert client.endpoint_style == "canonical"


def test_invalid_endpoint_style_override(monkeypatch):
    monkeypatch.setenv("MALVA_ENDPOINT_STYLE", "new-style")
    with pytest.raises(ValueError, match="MALVA_ENDPOINT_STYLE"):
        MalvaClient(base_url=BASE, api_token="tok")


def _federated_frame(shards):
    frame = pd.DataFrame({"provider_id": [], "sample_id": []})
    frame.attrs["shards"] = shards
    return frame


@responses.activate
def test_retrieve_cells_federated_frame_uses_completed_shard(monkeypatch):
    client = _client(monkeypatch, "canonical")
    frame = _federated_frame(
        [
            {
                "provider_id": "p1",
                "provider_name": "Provider One",
                "status": "failed",
                "search_job_id": None,
            },
            {
                "provider_id": "p2",
                "provider_name": "Provider Two",
                "status": "completed",
                "search_job_id": "job-b",
            },
        ]
    )

    responses.add(
        responses.POST,
        f"{BASE}/api/exports/cell-expression",
        json={"detail": "not found"},
        status=404,
    )
    responses.add(
        responses.GET,
        f"{BASE}/api/search/job-b",
        json={"detail": "not found"},
        status=404,
    )
    responses.add(
        responses.GET,
        f"{BASE}/api/search/job-b/cell-ids",
        json={"results": {"GENE": {"cell": [10], "sample": [1], "expression": [2.0]}}},
    )

    cells = client.retrieve_cells(
        frame,
        features="GENE",
        include_sample_metadata=False,
    )

    assert cells.to_dataframe()["value"].tolist() == [2.0]
    paths = [str(urlparse(call.request.url).path) for call in responses.calls]
    assert "/api/search/job-b/cell-ids" in paths
    assert not any("job-a" in path for path in paths)


@responses.activate
def test_retrieve_cells_federated_frame_honours_provider(monkeypatch):
    client = _client(monkeypatch, "canonical")
    frame = _federated_frame(
        [
            {
                "provider_id": "p1",
                "provider_name": "Provider One",
                "status": "completed",
                "search_job_id": "job-a",
            },
            {
                "provider_id": "p2",
                "provider_name": "Provider Two",
                "status": "completed",
                "search_job_id": "job-b",
            },
        ]
    )

    responses.add(
        responses.POST,
        f"{BASE}/api/exports/cell-expression",
        json={"detail": "not found"},
        status=404,
    )
    responses.add(
        responses.GET,
        f"{BASE}/api/search/job-b",
        json={"detail": "not found"},
        status=404,
    )
    responses.add(
        responses.GET,
        f"{BASE}/api/search/job-b/cell-ids",
        json={"results": {"GENE": {"cell": [10], "sample": [1], "expression": [2.0]}}},
    )

    cells = client.retrieve_cells(
        frame,
        features="GENE",
        include_sample_metadata=False,
        provider="Provider Two",
    )

    assert cells.to_dataframe()["value"].tolist() == [2.0]
    paths = [str(urlparse(call.request.url).path) for call in responses.calls]
    assert "/api/search/job-b/cell-ids" in paths
    assert not any("job-a" in path for path in paths)


def test_retrieve_cells_federated_frame_without_search_job_errors(monkeypatch):
    client = _client(monkeypatch, "canonical")
    frame = _federated_frame(
        [
            {
                "provider_id": "p1",
                "provider_name": "Provider One",
                "status": "failed",
                "search_job_id": None,
            }
        ]
    )

    with pytest.raises(MalvaAPIError, match="no completed shard with a search_job_id"):
        client.retrieve_cells(frame)


def test_retrieve_cells_federated_frame_without_shards_errors(monkeypatch):
    client = _client(monkeypatch, "canonical")

    with pytest.raises(MalvaAPIError, match="no shard metadata"):
        client.retrieve_cells(pd.DataFrame({"sample_id": [1]}))


def _not_found_callback(method, path, expected_query):
    def callback(request):
        _assert_request(request, method, path, expected_query)
        return 404, {"content-type": "application/json"}, json.dumps({"detail": "not found"})

    return callback


@responses.activate
def test_federated_frame_sends_provider_to_every_search_fetcher(monkeypatch):
    client = _client(monkeypatch, "canonical")
    frame = _federated_frame(
        [
            {
                "provider_id": "p1",
                "provider_name": "Provider One",
                "status": "completed",
                "search_job_id": "job-a",
            },
            {
                "provider_id": "p2",
                "provider_name": "Provider Two",
                "status": "completed",
                "search_job_id": "job-b",
            },
        ]
    )

    responses.add_callback(
        responses.GET,
        f"{BASE}/api/search/job-b",
        _not_found_callback("GET", "/api/search/job-b", {"max_cells": ["0"], "provider": ["p2"]}),
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/api/search/job-b/cell-ids",
        _not_found_callback("GET", "/api/search/job-b/cell-ids", {"provider": ["p2"]}),
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/api/search/job-b/global-cell-ids",
        _not_found_callback(
            "GET",
            "/api/search/job-b/global-cell-ids",
            {"provider": ["p2"], "sample_ids": ["1"]},
        ),
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/api/search/job-b/cells",
        _json_callback(
            "GET",
            "/api/search/job-b/cells",
            {"samples": {"1": {"cell_ids": [10]}}},
            {"provider": ["p2"], "sample_ids": ["1"]},
        ),
    )

    cells = client.retrieve_cells(
        frame,
        features="GENE",
        sample_ids=[1],
        provider="Provider Two",
        include_sample_metadata=False,
    )

    assert cells.source == "expression_cells"
    assert cells.to_dataframe()["value"].tolist() == [1.0]
    paths = [str(urlparse(call.request.url).path) for call in responses.calls]
    assert "/api/exports/cell-expression" not in paths
    assert not any("job-a" in path for path in paths)


@responses.activate
def test_federated_frame_without_provider_sends_no_provider_param(monkeypatch):
    client = _client(monkeypatch, "canonical")
    frame = _federated_frame(
        [
            {
                "provider_id": "p1",
                "provider_name": "Provider One",
                "status": "failed",
                "search_job_id": None,
            },
            {
                "provider_id": "p2",
                "provider_name": "Provider Two",
                "status": "completed",
                "search_job_id": "job-b",
            },
        ]
    )

    responses.add(
        responses.POST,
        f"{BASE}/api/exports/cell-expression",
        json={"detail": "not found"},
        status=404,
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/api/search/job-b",
        _not_found_callback("GET", "/api/search/job-b", {"max_cells": ["0"]}),
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/api/search/job-b/cell-ids",
        _json_callback(
            "GET",
            "/api/search/job-b/cell-ids",
            {"results": {"GENE": {"cell": [10], "sample": [1], "expression": [2.0]}}},
            {},
        ),
    )

    cells = client.retrieve_cells(
        frame,
        features="GENE",
        include_sample_metadata=False,
    )

    assert cells.source == "search_cell_ids"
    assert cells.to_dataframe()["value"].tolist() == [2.0]
    paths = [str(urlparse(call.request.url).path) for call in responses.calls]
    assert "/api/exports/cell-expression" in paths
    assert not any("job-a" in path for path in paths)


def test_federated_provider_error_names_the_selected_provider(monkeypatch):
    client = _client(monkeypatch, "canonical")
    frame = _federated_frame(
        [
            {
                "provider_id": "p1",
                "provider_name": "Provider One",
                "status": "failed",
                "search_job_id": None,
            }
        ]
    )

    with pytest.raises(MalvaAPIError, match=r"Provider One.*no completed search job"):
        client.retrieve_cells(frame, provider="Provider One")


@responses.activate
def test_federated_cell_metadata_enrichment_sends_provider(monkeypatch):
    client = _client(monkeypatch, "canonical")
    frame = _federated_frame(
        [
            {
                "provider_id": "p2",
                "provider_name": "Provider Two",
                "status": "completed",
                "search_job_id": "job-b",
            }
        ]
    )

    responses.add_callback(
        responses.GET,
        f"{BASE}/api/search/job-b",
        _not_found_callback("GET", "/api/search/job-b", {"max_cells": ["0"], "provider": ["p2"]}),
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/api/search/job-b/cell-ids",
        _json_callback(
            "GET",
            "/api/search/job-b/cell-ids",
            {"results": {"GENE": {"cell": [10], "sample": [1], "expression": [2.0]}}},
            {"provider": ["p2"]},
        ),
    )
    responses.add_callback(
        responses.POST,
        f"{BASE}/api/samples/metadata",
        _json_callback(
            "POST",
            "/api/samples/metadata",
            {"metadata": {"1": {"organ": "brain"}}},
            {"provider": ["p2"]},
        ),
    )

    cells = client.retrieve_cells(
        frame,
        features="GENE",
        sample_ids=[1],
        provider="p2",
    )

    assert cells.source == "search_cell_ids"
    assert cells.sample_metadata.to_dict("records") == [{"sample_id": 1, "organ": "brain"}]
    paths = [str(urlparse(call.request.url).path) for call in responses.calls]
    assert "/api/exports/cell-expression" not in paths
