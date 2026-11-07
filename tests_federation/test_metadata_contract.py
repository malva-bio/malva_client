"""Offline contract tests for provider-scoped Explorer metadata.

The Explorer metadata proxies route with ``?provider=<id>``.  These tests pin
the wire shape the legacy client must preserve: provider tags, local numeric
sample ids, explicit failures for unknown providers, and the legacy response
shape when no provider is requested.
"""

import json
from urllib.parse import parse_qs, urlparse

import pytest
import responses
from responses import matchers

from malva_client.client import MalvaClient
from malva_client.exceptions import MalvaAPIError

pytestmark = pytest.mark.federation

BASE = "https://federation.test"
PROVIDER_A = "malva-public-index"
PROVIDER_B = "singleron-aws"


def _query(request):
    return parse_qs(urlparse(request.url).query)


def _json_body(request):
    body = request.body
    if isinstance(body, bytes):
        body = body.decode()
    return json.loads(body or "{}")


def _provider_payload_complete(payload):
    """A provider payload is complete only when it explicitly says ``ok``."""
    return payload.get("ok") is True


@responses.activate
def test_db_stats_are_provider_tagged(client):
    responses.add(
        responses.GET,
        f"{BASE}/api/metadata/db-stats",
        json={
            "provider_id": PROVIDER_A,
            "provider_name": "Malva public index",
            "ok": True,
            "total_cells": 1234,
            "total_samples": 12,
        },
        match=[matchers.query_param_matcher({"provider": PROVIDER_A})],
    )

    payload = client._request(
        "GET",
        "/api/metadata/db-stats",
        params={"provider": PROVIDER_A},
    )

    assert payload["provider_id"] == PROVIDER_A
    assert payload["provider_name"] == "Malva public index"
    assert payload["ok"] is True
    assert payload["total_cells"] == 1234
    assert payload["total_samples"] == 12
    assert _query(responses.calls[0].request) == {"provider": [PROVIDER_A]}


@responses.activate
def test_celltype_sample_counts_keep_provider_and_local_ids(client):
    responses.add(
        responses.GET,
        f"{BASE}/api/metadata/celltype-sample-counts",
        json={
            "provider_id": PROVIDER_B,
            "provider_name": "Singleron",
            "ok": True,
            "counts": [[101, "T cell", 5], [102, "B cell", 2], [101, "B cell", 1]],
        },
        match=[matchers.query_param_matcher({"provider": PROVIDER_B})],
    )

    payload = client._request(
        "GET",
        "/api/metadata/celltype-sample-counts",
        params={"provider": PROVIDER_B},
    )

    assert payload["provider_id"] == PROVIDER_B
    assert payload["ok"] is True
    assert payload["counts"]

    for sample_id, cell_type, count in payload["counts"]:
        assert isinstance(sample_id, int)
        assert sample_id > 0
        assert isinstance(cell_type, str) and cell_type
        assert isinstance(count, int) and count >= 0

    assert {row[0] for row in payload["counts"]} == {101, 102}
    assert _query(responses.calls[0].request) == {"provider": [PROVIDER_B]}


@responses.activate
def test_group_payloads_carry_top_level_and_per_group_provider_tags(client):
    responses.add(
        responses.POST,
        f"{BASE}/api/metadata/groups",
        json={
            "provider_id": PROVIDER_A,
            "provider_name": "Malva public index",
            "columns": ["cell_type"],
            "groups": [
                {
                    "values": ["T cell"],
                    "n_samples": 3,
                    "provider_id": PROVIDER_A,
                    "provider_name": "Malva public index",
                }
            ],
            "total": 1,
        },
        match=[matchers.query_param_matcher({"provider": PROVIDER_A})],
    )

    payload = client._request(
        "POST",
        "/api/metadata/groups",
        params={"provider": PROVIDER_A},
        json={"columns": ["cell_type"]},
    )

    assert payload["provider_id"] == PROVIDER_A
    assert payload["provider_name"] == "Malva public index"
    assert payload["groups"][0]["provider_id"] == PROVIDER_A
    assert payload["groups"][0]["provider_name"] == "Malva public index"
    assert _json_body(responses.calls[0].request) == {"columns": ["cell_type"]}
    assert _query(responses.calls[0].request) == {"provider": [PROVIDER_A]}


@responses.activate
def test_sample_metadata_posts_local_numeric_ids_to_provider(client):
    responses.add(
        responses.POST,
        f"{BASE}/api/sample-metadata",
        json={
            "metadata": {
                "1": {"organ": "brain", "provider": PROVIDER_B},
                "2": {"organ": "liver", "provider": PROVIDER_B},
            }
        },
        match=[matchers.query_param_matcher({"provider": PROVIDER_B})],
    )

    payload = client._request(
        "POST",
        "/api/sample-metadata",
        params={"provider": PROVIDER_B},
        json={"ids": [1, 2]},
    )

    # Provider-local ids stay numeric on the wire; the Explorer must not be
    # asked to translate them into globally unique ids.
    assert _json_body(responses.calls[0].request) == {"ids": [1, 2]}
    assert set(payload["metadata"]) == {"1", "2"}
    assert {int(key) for key in payload["metadata"]} == {1, 2}
    assert _query(responses.calls[0].request) == {"provider": [PROVIDER_B]}


@responses.activate
def test_studies_hierarchy_and_drilldown_accept_provider_query(client):
    responses.add(
        responses.GET,
        f"{BASE}/samples/api/studies",
        json=["study-a", "study-b"],
        match=[matchers.query_param_matcher({"provider": PROVIDER_A})],
    )
    responses.add(
        responses.GET,
        f"{BASE}/samples/api/datasets/hierarchy",
        json={
            "datasets": [
                {
                    "dataset_id": "main",
                    "provider_id": PROVIDER_A,
                    "provider_name": "Malva public index",
                    "studies": [{"study": "study-a"}],
                }
            ]
        },
        match=[matchers.query_param_matcher({"provider": PROVIDER_A})],
    )
    responses.add(
        responses.GET,
        f"{BASE}/samples/api/datasets/main/studies",
        json={"studies": [{"study": "study-a"}], "total": 1},
        match=[
            matchers.query_param_matcher({"provider": PROVIDER_A, "page": "1", "page_size": "5"})
        ],
    )
    responses.add(
        responses.GET,
        f"{BASE}/samples/api/datasets/main/studies/study-a/samples",
        json={"samples": [{"sample_id": 7}], "total": 1},
        match=[
            matchers.query_param_matcher({"provider": PROVIDER_A, "page": "1", "page_size": "5"})
        ],
    )

    studies = client._request("GET", "/samples/api/studies", params={"provider": PROVIDER_A})
    hierarchy = client._request(
        "GET",
        "/samples/api/datasets/hierarchy",
        params={"provider": PROVIDER_A},
    )
    dataset_studies = client._request(
        "GET",
        "/samples/api/datasets/main/studies",
        params={"provider": PROVIDER_A, "page": 1, "page_size": 5},
    )
    samples = client._request(
        "GET",
        "/samples/api/datasets/main/studies/study-a/samples",
        params={"provider": PROVIDER_A, "page": 1, "page_size": 5},
    )

    assert studies == ["study-a", "study-b"]
    assert hierarchy["datasets"][0]["provider_id"] == PROVIDER_A
    assert dataset_studies["studies"][0]["study"] == "study-a"
    assert samples["samples"][0]["sample_id"] == 7

    for call in responses.calls:
        assert _query(call.request)["provider"] == [PROVIDER_A]


@responses.activate
def test_workspace_lists_metadata_providers_and_tags_datasets(client):
    responses.add(
        responses.GET,
        f"{BASE}/api/workspace",
        json={
            "metadata_providers": [
                {
                    "id": PROVIDER_A,
                    "name": "Malva public index",
                    "services": {"metadata": "http://public-metadata.test"},
                },
                {
                    "id": PROVIDER_B,
                    "name": "Singleron",
                    "services": {"metadata": "http://singleron-metadata.test"},
                },
            ],
            "providers": [
                {
                    "id": PROVIDER_A,
                    "name": "Malva public index",
                    "datasets": [
                        {
                            "dataset_id": "public-main",
                            "provider_id": PROVIDER_A,
                            "provider_name": "Malva public index",
                        }
                    ],
                },
                {
                    "id": PROVIDER_B,
                    "name": "Singleron",
                    "datasets": [
                        {
                            "dataset_id": "singleron-main",
                            "provider_id": PROVIDER_B,
                            "provider_name": "Singleron",
                        }
                    ],
                },
            ],
            "provider_ids": [PROVIDER_A, PROVIDER_B],
        },
        match=[matchers.query_param_matcher({})],
    )

    workspace = client._request("GET", "/api/workspace")

    providers = workspace["metadata_providers"]
    assert {provider["id"] for provider in providers} == {PROVIDER_A, PROVIDER_B}
    assert all(provider["services"]["metadata"] for provider in providers)

    datasets = [dataset for group in workspace["providers"] for dataset in group["datasets"]]
    assert {dataset["provider_id"] for dataset in datasets} == {PROVIDER_A, PROVIDER_B}
    assert all(dataset["provider_name"] for dataset in datasets)


@responses.activate
def test_search_only_provider_is_an_error_and_does_not_fall_back(client):
    responses.add(
        responses.GET,
        f"{BASE}/api/metadata/celltype-sample-counts",
        status=404,
        json={"error": 'provider "search-only" has no metadata service'},
        match=[matchers.query_param_matcher({"provider": "search-only"})],
    )
    # If the client silently dropped the provider param, this legacy response
    # would be returned instead.
    responses.add(
        responses.GET,
        f"{BASE}/api/metadata/celltype-sample-counts",
        json={"counts": [[1, "primary fallback", 1]]},
        match=[matchers.query_param_matcher({})],
    )

    with pytest.raises(MalvaAPIError, match="metadata"):
        client._request(
            "GET",
            "/api/metadata/celltype-sample-counts",
            params={"provider": "search-only"},
        )

    assert len(responses.calls) == 1
    assert _query(responses.calls[0].request) == {"provider": ["search-only"]}


@responses.activate
def test_unknown_provider_is_an_error_for_group_and_sample_metadata(client):
    responses.add(
        responses.POST,
        f"{BASE}/api/metadata/groups",
        status=404,
        json={"error": 'unknown provider "missing"'},
        match=[matchers.query_param_matcher({"provider": "missing"})],
    )
    responses.add(
        responses.POST,
        f"{BASE}/api/sample-metadata",
        status=404,
        json={"error": 'unknown provider "missing"'},
        match=[matchers.query_param_matcher({"provider": "missing"})],
    )

    with pytest.raises(MalvaAPIError, match="unknown provider"):
        client._request(
            "POST",
            "/api/metadata/groups",
            params={"provider": "missing"},
            json={"columns": ["cell_type"]},
        )
    with pytest.raises(MalvaAPIError, match="unknown provider"):
        client._request(
            "POST",
            "/api/sample-metadata",
            params={"provider": "missing"},
            json={"ids": [1]},
        )


@responses.activate
def test_soft_provider_failure_is_not_a_complete_payload(client):
    responses.add(
        responses.GET,
        f"{BASE}/api/metadata/celltype-sample-counts",
        json={
            "provider_id": PROVIDER_B,
            "provider_name": "Singleron",
            "ok": False,
            "counts": [],
            "error": "metadata service unavailable",
        },
        match=[matchers.query_param_matcher({"provider": PROVIDER_B})],
    )

    payload = client._request(
        "GET",
        "/api/metadata/celltype-sample-counts",
        params={"provider": PROVIDER_B},
    )

    assert payload["ok"] is False
    assert payload["provider_id"] == PROVIDER_B
    assert payload["counts"] == []
    assert payload["error"]
    assert not _provider_payload_complete(payload)


@responses.activate
def test_no_provider_preserves_the_legacy_request_shape(client):
    responses.add(
        responses.GET,
        f"{BASE}/api/metadata/db-stats",
        json={"total_cells": 10, "total_samples": 2},
        match=[matchers.query_param_matcher({})],
    )
    responses.add(
        responses.GET,
        f"{BASE}/api/metadata/celltype-sample-counts",
        json={"counts": [[1, "T cell", 4]]},
        match=[matchers.query_param_matcher({})],
    )
    responses.add(
        responses.POST,
        f"{BASE}/api/sample-metadata",
        json={"metadata": {"1": {"organ": "brain"}}},
        match=[matchers.query_param_matcher({})],
    )

    stats = client._request("GET", "/api/metadata/db-stats")
    counts = client._request("GET", "/api/metadata/celltype-sample-counts")
    metadata = client.get_sample_metadata([1])

    assert stats == {"total_cells": 10, "total_samples": 2}
    assert counts == {"counts": [[1, "T cell", 4]]}
    assert metadata == {1: {"organ": "brain"}}

    assert all(_query(call.request) == {} for call in responses.calls)


@pytest.mark.skipif(
    not hasattr(MalvaClient, "search_federated"),
    reason="search_federated is not present in this checkout",
)
@responses.activate
def test_federated_search_keeps_provider_columns(client):
    shards = [
        {
            "provider_id": PROVIDER_A,
            "provider_name": "Malva public index",
            "search_job_id": "job-a",
            "status": "completed",
            "error": None,
        },
        {
            "provider_id": PROVIDER_B,
            "provider_name": "Singleron",
            "search_job_id": "job-b",
            "status": "completed",
            "error": None,
        },
    ]

    def genes(sample_id):
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

    responses.add(
        responses.POST,
        f"{BASE}/api/expression/submit",
        json={"job_id": "fed-1", "status": "running", "shards": shards},
    )
    responses.add(
        responses.GET,
        f"{BASE}/api/expression/federated/fed-1",
        json={"job_id": "fed-1", "status": "completed", "shards": shards},
    )
    responses.add(
        responses.GET,
        f"{BASE}/api/expression/results/job-a",
        json={"status": "completed", "genes": genes(101)},
        match=[matchers.query_param_matcher({"provider": PROVIDER_A})],
    )
    responses.add(
        responses.GET,
        f"{BASE}/api/expression/results/job-b",
        json={"status": "completed", "genes": genes(1)},
        match=[matchers.query_param_matcher({"provider": PROVIDER_B})],
    )

    frame = client.search_federated("FOXP3")

    assert frame.columns[:2].tolist() == ["provider", "provider_id"]
    assert frame["provider"].tolist() == ["Malva public index", "Singleron"]
    assert frame["provider_id"].tolist() == [PROVIDER_A, PROVIDER_B]
    assert frame["sample_id"].tolist() == [101, 1]
