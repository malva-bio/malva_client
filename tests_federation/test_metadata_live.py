"""Live tests against a federation-enabled Explorer.

These tests are opt-in and require::

    MALVA_FEDERATION_API_URL
    MALVA_FEDERATION_API_TOKEN
    MALVA_FEDERATION_GENE   # optional, defaults to FOXP3

They skip when the credentials are absent, and also when the Explorer exposes
fewer than two providers with a metadata service. Assertions cover contract
shape only; exact counts and study names are intentionally not pinned.
"""

import os
from dataclasses import dataclass
from urllib.parse import quote

import pytest

from malva_client.client import MalvaClient

pytestmark = [pytest.mark.federation, pytest.mark.federation_live]

DEFAULT_GENE = "FOXP3"


@dataclass(frozen=True)
class FederationSession:
    client: MalvaClient
    workspace: dict
    providers: list
    provider_ids: list


def _metadata_providers(workspace):
    providers = workspace.get("metadata_providers") or []
    return [
        provider
        for provider in providers
        if isinstance(provider, dict)
        and provider.get("id")
        and (provider.get("services") or {}).get("metadata")
    ]


def _datasets_for_provider(session, provider_id):
    if provider_id not in session.provider_ids:
        return []
    for group in session.workspace.get("providers") or []:
        if isinstance(group, dict) and group.get("id") == provider_id:
            return [row for row in group.get("datasets") or [] if isinstance(row, dict)]
    return []


def _dataset_id(row):
    if not isinstance(row, dict):
        return None
    value = row.get("dataset_id") or row.get("id") or row.get("dataset")
    return str(value) if value else None


def _dataset_rows(payload):
    if isinstance(payload, list):
        return [row for row in payload if isinstance(row, dict)]
    if isinstance(payload, dict):
        for key in ("datasets", "hierarchy", "data", "items"):
            rows = payload.get(key)
            if isinstance(rows, list):
                return [row for row in rows if isinstance(row, dict)]
    return []


def _study_rows(payload):
    if isinstance(payload, list):
        return [row for row in payload if isinstance(row, (dict, str))]
    if isinstance(payload, dict):
        for key in ("studies", "data", "items"):
            rows = payload.get(key)
            if isinstance(rows, list):
                return [row for row in rows if isinstance(row, (dict, str))]
    return []


def _study_name(row):
    if isinstance(row, str):
        return row
    if isinstance(row, dict):
        value = row.get("study") or row.get("study_name") or row.get("name") or row.get("id")
        return str(value) if value else None
    return None


def _sample_ids_from_counts(payload, limit=3):
    sample_ids = []
    for entry in payload.get("counts") or []:
        if not isinstance(entry, (list, tuple)) or len(entry) < 3:
            continue
        try:
            sample_id = int(entry[0])
        except (TypeError, ValueError):
            continue
        if sample_id not in sample_ids:
            sample_ids.append(sample_id)
        if len(sample_ids) >= limit:
            break
    return sample_ids


@pytest.fixture(scope="session")
def federation_session():
    api_url = os.environ.get("MALVA_FEDERATION_API_URL")
    api_token = os.environ.get("MALVA_FEDERATION_API_TOKEN")
    if not api_url or not api_token:
        pytest.skip(
            "MALVA_FEDERATION_API_URL and MALVA_FEDERATION_API_TOKEN are required "
            "for live federation tests"
        )

    client = MalvaClient(base_url=api_url, api_token=api_token, timeout=120)
    workspace = client._request("GET", "/api/workspace")
    providers = _metadata_providers(workspace)
    if len(providers) < 2:
        pytest.skip(
            "federation is not enabled: fewer than two providers expose "
            f"services.metadata (found {len(providers)})"
        )

    provider_ids = [provider["id"] for provider in providers]
    return FederationSession(
        client=client,
        workspace=workspace,
        providers=providers,
        provider_ids=provider_ids,
    )


@pytest.fixture(scope="session")
def provider_counts(federation_session):
    payloads = {}
    for provider_id in federation_session.provider_ids:
        payloads[provider_id] = federation_session.client._request(
            "GET",
            "/api/metadata/celltype-sample-counts",
            params={"provider": provider_id},
        )
    return payloads


@pytest.fixture(scope="session")
def sample_ids_by_provider(provider_counts):
    return {
        provider_id: _sample_ids_from_counts(payload)
        for provider_id, payload in provider_counts.items()
    }


def test_workspace_lists_multiple_metadata_providers(federation_session):
    assert len(federation_session.provider_ids) >= 2
    assert len(set(federation_session.provider_ids)) == len(federation_session.provider_ids)
    for provider in federation_session.providers:
        assert provider["services"]["metadata"]

    for provider_id in federation_session.provider_ids:
        for dataset in _datasets_for_provider(federation_session, provider_id):
            assert dataset.get("provider_id") == provider_id
            assert dataset.get("provider_name")


def test_db_stats_are_provider_scoped(federation_session):
    seen = set()
    for provider_id in federation_session.provider_ids:
        payload = federation_session.client._request(
            "GET",
            "/api/metadata/db-stats",
            params={"provider": provider_id},
        )
        assert payload["provider_id"] == provider_id
        assert payload["provider_name"]
        assert payload["ok"] is True
        assert int(payload["total_cells"]) >= 0
        assert int(payload["total_samples"]) >= 0
        seen.add(provider_id)

    assert seen == set(federation_session.provider_ids)


def test_celltype_sample_counts_are_provider_scoped(federation_session, provider_counts):
    for provider_id, payload in provider_counts.items():
        assert payload["provider_id"] == provider_id
        assert payload["provider_name"]
        assert payload["ok"] is True
        assert isinstance(payload.get("counts"), list)

        for entry in (payload.get("counts") or [])[:50]:
            assert isinstance(entry, (list, tuple))
            assert len(entry) >= 3
            assert isinstance(entry[0], int) and entry[0] > 0
            assert isinstance(entry[1], str) and entry[1]
            assert isinstance(entry[2], int) and entry[2] >= 0


def test_group_payloads_are_provider_tagged(federation_session):
    for provider_id in federation_session.provider_ids:
        payload = federation_session.client._request(
            "POST",
            "/api/metadata/groups",
            params={"provider": provider_id},
            json={"columns": ["cell_type"], "limit": 10},
        )
        assert payload["provider_id"] == provider_id
        assert payload["provider_name"]
        for group in payload.get("groups") or []:
            assert group["provider_id"] == provider_id
            assert group["provider_name"]


def test_sample_metadata_uses_provider_local_numeric_ids(
    federation_session, sample_ids_by_provider
):
    checked = 0
    for provider_id, local_ids in sample_ids_by_provider.items():
        if not local_ids:
            continue
        payload = federation_session.client._request(
            "POST",
            "/api/sample-metadata",
            params={"provider": provider_id},
            json={"ids": local_ids},
        )
        metadata = payload.get("metadata") or {}
        assert metadata, f"{provider_id} returned no sample metadata"

        metadata_ids = set()
        for key in metadata:
            metadata_ids.add(int(key))
        assert metadata_ids <= set(local_ids)
        checked += 1

    if checked == 0:
        pytest.skip("no provider exposed sample counts to check metadata")


def test_studies_hierarchy_and_drilldown_are_provider_scoped(federation_session):
    checked_hierarchy = 0
    checked_studies = 0
    checked_samples = 0

    for provider_id in federation_session.provider_ids:
        studies = federation_session.client._request(
            "GET",
            "/samples/api/studies",
            params={"provider": provider_id},
        )
        assert isinstance(studies, (list, dict))

        hierarchy = federation_session.client._request(
            "GET",
            "/samples/api/datasets/hierarchy",
            params={"provider": provider_id},
        )
        assert isinstance(hierarchy, (dict, list))
        checked_hierarchy += 1

        dataset_ids = [
            dataset_id
            for dataset_id in (
                _dataset_id(row) for row in _datasets_for_provider(federation_session, provider_id)
            )
            if dataset_id
        ]
        dataset_ids.extend(
            dataset_id
            for dataset_id in (_dataset_id(row) for row in _dataset_rows(hierarchy))
            if dataset_id
        )
        if not dataset_ids:
            continue

        dataset_id = dataset_ids[0]
        dataset_studies = federation_session.client._request(
            "GET",
            f"/samples/api/datasets/{quote(dataset_id, safe='')}/studies",
            params={"provider": provider_id, "page": 1, "page_size": 5},
        )
        assert isinstance(dataset_studies, (dict, list))
        checked_studies += 1

        study_rows = _study_rows(dataset_studies)
        study_name = _study_name(study_rows[0]) if study_rows else None
        if not study_name:
            continue

        samples = federation_session.client._request(
            "GET",
            (
                f"/samples/api/datasets/{quote(dataset_id, safe='')}"
                f"/studies/{quote(study_name, safe='')}/samples"
            ),
            params={"provider": provider_id, "page": 1, "page_size": 5},
        )
        assert isinstance(samples, (dict, list))
        checked_samples += 1

    assert checked_hierarchy == len(federation_session.provider_ids)
    assert checked_studies > 0
    assert checked_samples > 0


def test_namespaced_sample_keys_do_not_collide(sample_ids_by_provider):
    namespaced = {
        f"{provider_id}::{sample_id}"
        for provider_id, sample_ids in sample_ids_by_provider.items()
        for sample_id in sample_ids
    }
    expected = sum(len(set(sample_ids)) for sample_ids in sample_ids_by_provider.values())
    if not namespaced:
        pytest.skip("no provider exposed sample counts to namespace")

    assert len(namespaced) == expected
    assert all("::" in key for key in namespaced)

    local_ids = {
        sample_id for sample_ids in sample_ids_by_provider.values() for sample_id in sample_ids
    }
    for local_id in local_ids:
        keys = {
            f"{provider_id}::{local_id}"
            for provider_id, sample_ids in sample_ids_by_provider.items()
            if local_id in sample_ids
        }
        assert len(keys) == sum(
            1 for sample_ids in sample_ids_by_provider.values() if local_id in sample_ids
        )


def _live_federated_rows(client, gene, provider_ids):
    """Exercise the federated endpoints when search_federated is unavailable."""
    submitted = client._request(
        "POST",
        "/api/expression/submit",
        json={"query": gene, "provider_ids": provider_ids},
    )
    job_id = submitted["job_id"]
    state = submitted
    terminal = {"completed", "failed", "error", "canceled", "interrupted"}
    shards = submitted.get("shards") or []
    if not (shards and all(shard.get("status") in terminal for shard in shards)):
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
            expression = gene_data.get("expression_data") or {}
            samples = expression.get("samples") or gene_data.get("samples") or []
            for sample_id in samples:
                rows.append(
                    {
                        "provider": provider_name,
                        "provider_id": provider_id,
                        "sample_id": int(sample_id),
                    }
                )
    return rows


def test_federated_search_carries_provider_columns(federation_session):
    gene = os.environ.get("MALVA_FEDERATION_GENE") or DEFAULT_GENE

    if hasattr(federation_session.client, "search_federated"):
        frame = federation_session.client.search_federated(
            gene,
            provider_ids=federation_session.provider_ids,
        )
        assert frame.columns[:2].tolist() == ["provider", "provider_id"]
        if not frame.empty:
            assert set(frame["provider_id"].dropna()) <= set(federation_session.provider_ids)
            assert frame["provider"].notna().all()
        return

    rows = _live_federated_rows(
        federation_session.client,
        gene,
        federation_session.provider_ids,
    )
    assert rows
    assert {row["provider_id"] for row in rows} <= set(federation_session.provider_ids)
    assert all(row["provider"] for row in rows)
    assert all(isinstance(row["sample_id"], int) for row in rows)
