"""Federated Explorer search provenance, partial failures and validation."""

import pytest

from malva_client.client import _FEDERATED_RESULT_COLUMNS, MalvaClient
from malva_client.exceptions import MalvaAPIError


def _genes(sample_id: int) -> dict:
    return {
        "FOXP3": {
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


def _client(monkeypatch, responses):
    client = MalvaClient.__new__(MalvaClient)
    client.base_url = "http://explorer.example"
    client.dataset_id = None
    client.provider_id = None
    client._endpoint_style = "legacy"
    calls = []

    def fake_request(method, endpoint, **kwargs):
        calls.append((method, endpoint, kwargs))
        return responses(endpoint, kwargs)

    monkeypatch.setattr(client, "_request", fake_request)
    return client, calls


def _shards(*, second_status="completed", second_error=None, second_job="remote-1"):
    return [
        {
            "provider_id": "malva-public-index",
            "provider_name": "Malva public index",
            "search_job_id": "pub-1",
            "status": "completed",
            "error": None,
        },
        {
            "provider_id": "singleron-aws",
            "provider_name": "Singleron",
            "search_job_id": second_job,
            "status": second_status,
            "error": second_error,
        },
    ]


def test_search_federated_preserves_provenance_and_surfaces_failed_shard(monkeypatch):
    shards = _shards(second_status="error", second_error="403 Forbidden", second_job=None)

    def responses(endpoint, kwargs):
        if endpoint == "/api/expression/submit":
            return {"job_id": "fed-1", "status": "failed", "shards": shards}
        if endpoint == "/api/expression/results/pub-1":
            assert kwargs["params"] == {"provider": "malva-public-index"}
            return {"status": "completed", "genes": _genes(101)}
        raise AssertionError(endpoint)

    client, _ = _client(monkeypatch, responses)

    with pytest.warns(RuntimeWarning, match="Singleron: 403 Forbidden"):
        frame = client.search_federated("FOXP3", dataset_id="main")

    assert frame.columns[:2].tolist() == ["provider", "provider_id"]
    assert frame["provider"].tolist() == ["Malva public index"]
    assert frame.attrs["job_id"] == "fed-1"
    assert len(frame.attrs["shards"]) == 2
    assert frame.attrs["errors"] == ["Singleron: 403 Forbidden"]
    assert frame.attrs["shard_errors"][0]["provider_id"] == "singleron-aws"
    statuses = {item["provider_id"]: item["status"] for item in frame.attrs["providers"]}
    assert statuses == {"malva-public-index": "completed", "singleron-aws": "error"}


def test_search_federated_strict_raises_for_partial_failure(monkeypatch):
    shards = _shards(second_status="error", second_error="403 Forbidden", second_job=None)

    def responses(endpoint, kwargs):
        if endpoint == "/api/expression/submit":
            return {"job_id": "fed-1", "status": "failed", "shards": shards}
        if endpoint == "/api/expression/results/pub-1":
            return {"status": "completed", "genes": _genes(101)}
        raise AssertionError(endpoint)

    client, _ = _client(monkeypatch, responses)

    with pytest.warns(RuntimeWarning), pytest.raises(MalvaAPIError, match="403 Forbidden"):
        client.search_federated("FOXP3", strict=True)


def test_search_federated_empty_success_keeps_schema_and_provenance(monkeypatch):
    shards = _shards()

    def responses(endpoint, kwargs):
        if endpoint == "/api/expression/submit":
            return {"job_id": "fed-empty", "status": "completed", "shards": shards}
        if endpoint in {
            "/api/expression/results/pub-1",
            "/api/expression/results/remote-1",
        }:
            return {"status": "completed", "genes": {}}
        raise AssertionError(endpoint)

    client, _ = _client(monkeypatch, responses)
    frame = client.search_federated("FOXP3")

    assert frame.empty
    assert frame.columns.tolist() == _FEDERATED_RESULT_COLUMNS
    assert frame.attrs["job_id"] == "fed-empty"
    assert frame.attrs["errors"] == []
    assert len(frame.attrs["shards"]) == 2
    assert {item["provider_id"] for item in frame.attrs["providers"]} == {
        "malva-public-index",
        "singleron-aws",
    }


def test_search_federated_rejects_unknown_provider_before_submit(monkeypatch):
    def responses(endpoint, kwargs):
        if endpoint == "/api/workspace":
            return {
                "providers": [
                    {"id": "malva-public-index", "name": "Malva public index"},
                ],
                "provider_ids": ["malva-public-index"],
            }
        raise AssertionError(endpoint)

    client, calls = _client(monkeypatch, responses)

    with pytest.raises(MalvaAPIError, match="Unknown federation provider"):
        client.search_federated("FOXP3", provider_ids=["does-not-exist"])

    assert [call[1] for call in calls] == ["/api/workspace"]


def test_search_federated_reports_zero_shard_running_job(monkeypatch):
    def responses(endpoint, kwargs):
        if endpoint == "/api/expression/submit":
            return {"job_id": "fed-zero", "status": "running", "shards": []}
        raise AssertionError(endpoint)

    client, _ = _client(monkeypatch, responses)

    with pytest.raises(MalvaAPIError, match="returned no shards"):
        client.search_federated("FOXP3", max_wait=0)


def test_search_federated_reports_all_failed_shards(monkeypatch):
    shards = [
        {
            "provider_id": "singleron-aws",
            "provider_name": "Singleron",
            "search_job_id": None,
            "status": "error",
            "error": "403 Forbidden",
        }
    ]

    def responses(endpoint, kwargs):
        if endpoint == "/api/expression/submit":
            return {"job_id": "fed-2", "status": "failed", "shards": shards}
        raise AssertionError(endpoint)

    client, _ = _client(monkeypatch, responses)

    with (
        pytest.warns(RuntimeWarning),
        pytest.raises(MalvaAPIError, match="Singleron: 403 Forbidden"),
    ):
        client.search_federated("FOXP3")
