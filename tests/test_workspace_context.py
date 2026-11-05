from unittest.mock import patch

from malva_client.client import MalvaClient
from malva_client.config import Config


def _client(**kwargs):
    with patch.object(MalvaClient, "_test_connection"):
        return MalvaClient(base_url="https://test.malva.bio", api_token="tok", **kwargs)


def test_workspace_context_is_sent_as_headers():
    client = _client(workspace="atlas:human", dataset_id="human", provider_id="p1")

    assert client.session.headers["X-Malva-Workspace"] == "atlas:human"
    assert client.session.headers["X-Malva-Dataset"] == "human"
    assert client.session.headers["X-Malva-Provider"] == "p1"

    client.set_workspace(workspace="atlas:mouse", dataset_id="mouse", provider_id="")
    assert client.session.headers["X-Malva-Workspace"] == "atlas:mouse"
    assert client.session.headers["X-Malva-Dataset"] == "mouse"
    assert "X-Malva-Provider" not in client.session.headers


def test_search_payload_uses_the_active_dataset():
    client = _client(dataset_id="human")
    payload = client._build_search_payload("TP53", False, 10, 1, True, None, 0, 100)

    assert payload["dataset_id"] == "human"


def test_config_round_trips_workspace_context(tmp_path):
    config = Config(str(tmp_path / "config.json"))
    config.server_url = "https://test.malva.bio"
    config.workspace = "atlas:human"
    config.dataset_id = "human"
    config.provider_id = "p1"
    config.save()

    loaded = Config.load(str(tmp_path / "config.json"))
    assert loaded.workspace == "atlas:human"
    assert loaded.dataset_id == "human"
    assert loaded.provider_id == "p1"


def test_get_datasets_lists_search_indexes():
    client = _client(dataset_id="human")

    with patch.object(
        client,
        "_request",
        return_value=[{"dataset_id": "human", "name": "Human"}],
    ) as request:
        datasets = client.get_datasets()

    assert datasets == [{"dataset_id": "human", "name": "Human"}]
    request.assert_called_once_with("GET", "/datasets")
