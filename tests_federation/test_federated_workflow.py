"""Offline end-to-end shape test for a federated Explorer workflow.

A federated search fans out to several providers and keeps each shard's origin.
Metadata is provider-local, so the join key must be namespaced
(``provider::sample``) or samples that share a numeric id collide.
"""

from urllib.parse import parse_qs, urlparse

import pytest
import responses
from responses import matchers

from malva_client.client import MalvaClient

pytestmark = pytest.mark.federation

BASE = "https://federation.test"
PROVIDER_A = "malva-public-index"
PROVIDER_B = "singleron-aws"
LOCAL_SAMPLE_ID = 7


def _query(request):
    return parse_qs(urlparse(request.url).query)


def _genes(sample_id):
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


def _samples_from_gene(gene_data):
    expression = gene_data.get("expression_data") or {}
    return expression.get("samples") or gene_data.get("samples") or []


def _rows_from_explorer(client):
    """Build result rows with the same provider contract as search_federated."""
    submitted = client._request("POST", "/api/expression/submit", json={"query": "FOXP3"})
    job_id = submitted["job_id"]
    state = submitted
    if not all(shard.get("status") == "completed" for shard in state.get("shards") or []):
        state = client._request("GET", f"/api/expression/federated/{job_id}")

    rows = []
    for shard in state.get("shards") or []:
        if shard.get("status") != "completed" or not shard.get("search_job_id"):
            continue
        provider_id = shard["provider_id"]
        provider_name = shard.get("provider_name") or provider_id
        data = client._request(
            "GET",
            f"/api/expression/results/{shard['search_job_id']}",
            params={"provider": provider_id},
        )
        genes = data.get("genes") or data.get("results") or {}
        for gene_data in genes.values():
            for sample_id in _samples_from_gene(gene_data):
                rows.append(
                    {
                        "provider": provider_name,
                        "provider_id": provider_id,
                        "sample_id": int(sample_id),
                    }
                )
    return rows


def _register_federated_responses():
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
    # Both providers return the same provider-local sample id on purpose.
    responses.add(
        responses.GET,
        f"{BASE}/api/expression/results/job-a",
        json={"status": "completed", "genes": _genes(LOCAL_SAMPLE_ID)},
        match=[matchers.query_param_matcher({"provider": PROVIDER_A})],
    )
    responses.add(
        responses.GET,
        f"{BASE}/api/expression/results/job-b",
        json={"status": "completed", "genes": _genes(LOCAL_SAMPLE_ID)},
        match=[matchers.query_param_matcher({"provider": PROVIDER_B})],
    )
    responses.add(
        responses.POST,
        f"{BASE}/api/sample-metadata",
        json={"metadata": {str(LOCAL_SAMPLE_ID): {"organ": "brain"}}},
        match=[matchers.query_param_matcher({"provider": PROVIDER_A})],
    )
    responses.add(
        responses.POST,
        f"{BASE}/api/sample-metadata",
        json={"metadata": {str(LOCAL_SAMPLE_ID): {"organ": "liver"}}},
        match=[matchers.query_param_matcher({"provider": PROVIDER_B})],
    )


@responses.activate
def test_federated_rows_join_provider_scoped_metadata_without_collisions(client):
    _register_federated_responses()

    if hasattr(client, "search_federated"):
        frame = client.search_federated("FOXP3")
        rows = [
            {
                "provider": record["provider"],
                "provider_id": record["provider_id"],
                "sample_id": int(record["sample_id"]),
            }
            for record in frame.to_dict("records")
        ]
    else:
        rows = _rows_from_explorer(client)

    assert rows
    assert {row["provider_id"] for row in rows} == {PROVIDER_A, PROVIDER_B}
    assert {row["provider"] for row in rows} == {"Malva public index", "Singleron"}
    assert {row["sample_id"] for row in rows} == {LOCAL_SAMPLE_ID}

    joined = {}
    for provider_id in (PROVIDER_A, PROVIDER_B):
        local_ids = sorted({row["sample_id"] for row in rows if row["provider_id"] == provider_id})
        payload = client._request(
            "POST",
            "/api/sample-metadata",
            params={"provider": provider_id},
            json={"ids": local_ids},
        )
        for sample_id, values in payload["metadata"].items():
            joined[f"{provider_id}::{int(sample_id)}"] = values

    assert set(joined) == {
        f"{PROVIDER_A}::{LOCAL_SAMPLE_ID}",
        f"{PROVIDER_B}::{LOCAL_SAMPLE_ID}",
    }
    assert joined[f"{PROVIDER_A}::{LOCAL_SAMPLE_ID}"]["organ"] == "brain"
    assert joined[f"{PROVIDER_B}::{LOCAL_SAMPLE_ID}"]["organ"] == "liver"
    assert len(joined) == len(rows)

    provider_queries = {
        _query(call.request)["provider"][0]
        for call in responses.calls
        if "/api/expression/results/" in call.request.url
        or "/api/sample-metadata" in call.request.url
    }
    assert provider_queries == {PROVIDER_A, PROVIDER_B}


@pytest.mark.skipif(
    not hasattr(MalvaClient, "search_federated"),
    reason="search_federated is not present in this checkout",
)
@responses.activate
def test_search_federated_helper_returns_provider_columns(client):
    _register_federated_responses()

    frame = client.search_federated("FOXP3")

    assert frame.columns[:2].tolist() == ["provider", "provider_id"]
    assert set(frame["provider_id"]) == {PROVIDER_A, PROVIDER_B}
    assert set(frame["sample_id"]) == {LOCAL_SAMPLE_ID}
