"""Provider scoping for metadata and sample methods."""

import json
from unittest.mock import patch
from urllib.parse import parse_qs, urlparse

import msgpack
import numpy as np
import responses

from malva_client.client import MalvaClient

BASE = "https://test.malva.bio"


def _client(**kwargs):
    with patch.object(MalvaClient, "_test_connection"):
        client = MalvaClient(base_url=BASE, api_token="tok", **kwargs)
    client._endpoint_style = "legacy"
    return client


def _query(request):
    return parse_qs(urlparse(request.url).query)


def _assert_provider(request, expected):
    assert _query(request).get("provider") == [expected]
    return expected


@responses.activate
def test_client_default_provider_is_sent_as_query_param():
    client = _client(provider_id="default-provider")

    def samples(request):
        _assert_provider(request, "default-provider")
        return 200, {}, json.dumps({"samples": [{"sample_id": 1}]})

    def filters(request):
        _assert_provider(request, "default-provider")
        return 200, {}, json.dumps({"values": ["brain"]})

    responses.add_callback(responses.GET, f"{BASE}/samples/api/search", callback=samples)
    responses.add_callback(
        responses.GET,
        f"{BASE}/samples/api/filters/organ/values",
        callback=filters,
    )

    assert client.get_samples()["samples"] == [{"sample_id": 1}]
    assert client.get_filter_values("organ") == ["brain"]


@responses.activate
def test_explicit_provider_wins_over_client_default():
    client = _client(provider_id="default-provider")

    def samples(request):
        _assert_provider(request, "explicit-provider")
        return 200, {}, json.dumps({"samples": []})

    responses.add_callback(responses.GET, f"{BASE}/samples/api/search", callback=samples)

    assert client.get_samples(provider="explicit-provider")["samples"] == []


@responses.activate
def test_set_workspace_updates_provider_default():
    client = _client()
    client.set_workspace(provider_id="workspace-provider")

    def filters(request):
        _assert_provider(request, "workspace-provider")
        return 200, {}, json.dumps({"values": []})

    responses.add_callback(
        responses.GET,
        f"{BASE}/samples/api/filters/organ/values",
        callback=filters,
    )

    assert client.get_filter_values("organ") == []


@responses.activate
def test_dataset_and_study_methods_send_provider():
    client = _client(provider_id="default-provider")

    def hierarchy(request):
        _assert_provider(request, "default-provider")
        return 200, {}, json.dumps({"datasets": []})

    def studies(request):
        _assert_provider(request, "explicit-provider")
        return 200, {}, json.dumps({"studies": [{"study_name": "S1"}]})

    def samples(request):
        _assert_provider(request, "default-provider")
        return 200, {}, json.dumps({"samples": [{"sample_id": 10}]})

    def overview(request):
        _assert_provider(request, "default-provider")
        return 200, {}, json.dumps({"total_studies": 1})

    def details(request):
        _assert_provider(request, "explicit-provider")
        return 200, {}, json.dumps({"sample_uuid": "s1"})

    responses.add_callback(
        responses.GET, f"{BASE}/samples/api/datasets/hierarchy", callback=hierarchy
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/samples/api/datasets/d1/studies",
        callback=studies,
    )
    responses.add_callback(
        responses.GET,
        f"{BASE}/samples/api/datasets/d1/studies/S1/samples",
        callback=samples,
    )
    responses.add_callback(responses.GET, f"{BASE}/samples/api/overview", callback=overview)
    responses.add_callback(responses.GET, f"{BASE}/samples/api/sample/s1", callback=details)

    assert client.get_datasets_hierarchy() == {"datasets": []}
    assert client.get_dataset_studies("d1", provider="explicit-provider")["studies"] == [
        {"study_name": "S1"}
    ]
    assert client.get_study_samples("d1", "S1")["samples"] == [{"sample_id": 10}]
    assert client.get_overview_stats()["total_studies"] == 1
    assert client.get_sample_details("s1", provider="explicit-provider") == {"sample_uuid": "s1"}


@responses.activate
def test_sample_metadata_post_sends_provider():
    client = _client(provider_id="default-provider")

    def metadata(request):
        _assert_provider(request, "explicit-provider")
        return 200, {}, json.dumps({"metadata": {"7": {"organ": "brain"}}})

    responses.add_callback(responses.POST, f"{BASE}/api/sample-metadata", callback=metadata)

    assert client.get_sample_metadata([7], provider="explicit-provider")[7] == {"organ": "brain"}


@responses.activate
def test_cells_by_metadata_post_sends_provider():
    client = _client(provider_id="default-provider")
    payload = msgpack.packb(
        {
            "c": np.array([101], dtype=np.int64).tobytes(),
            "s": np.array([7], dtype=np.int64).tobytes(),
            "n": 1,
        },
        use_bin_type=True,
    )

    def cells(request):
        _assert_provider(request, "explicit-provider")
        return 200, {"Content-Type": "application/x-msgpack"}, payload

    responses.add_callback(responses.POST, f"{BASE}/api/cells/by-filter", callback=cells)

    frame = client.get_cells_by_metadata(
        sample_ids=[7], include_cell_metadata=False, provider="explicit-provider"
    )
    assert frame.to_dict("records") == [{"sample_id": 7, "cell_id": 101}]


def test_provider_argument_can_disable_default():
    client = _client(provider_id="default-provider")
    assert client._provider_params("") == {}
    assert client._provider_params(None) == {"provider": "default-provider"}
    assert client._provider_params("explicit") == {"provider": "explicit"}
