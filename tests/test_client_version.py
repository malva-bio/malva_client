"""Tests for the client version header and the server minimum-version gate."""

import pytest
import responses

from malva_client import __version__
from malva_client.client import MalvaClient
from malva_client.exceptions import ClientVersionError


@pytest.fixture
def client():
    with responses.RequestsMock() as rsps:
        rsps.add(
            responses.GET,
            "https://malva.test/api/health",
            json={"status": "ok"},
            status=200,
        )
        rsps.add(
            responses.GET,
            "https://malva.test/api/account/quota",
            json={"account_type": "enterprise"},
            status=200,
        )
        yield MalvaClient(base_url="https://malva.test", api_token="token")


def test_client_sends_version_header(client):
    assert client.session.headers["X-Malva-Client-Version"] == __version__
    assert client.session.headers["User-Agent"] == f"malva-client/{__version__}"


@responses.activate
def test_minimum_version_error_is_raised(client):
    responses.add(
        responses.GET,
        "https://malva.test/api/account/quota",
        json={
            "error": "Client version 0.3.4 is no longer supported. Minimum required: 0.4.0",
            "min_client_version": "0.4.0",
            "client_version": "0.3.4",
        },
        status=426,
    )

    with pytest.raises(ClientVersionError) as excinfo:
        client._request("GET", "/api/account/quota")

    assert excinfo.value.min_version == "0.4.0"
    assert excinfo.value.client_version == "0.3.4"
    assert "no longer supported" in str(excinfo.value)


@responses.activate
def test_minimum_version_error_without_body(client):
    responses.add(
        responses.GET,
        "https://malva.test/health",
        body="upgrade required",
        status=426,
    )

    with pytest.raises(ClientVersionError) as excinfo:
        client._request("GET", "/health")

    assert excinfo.value.min_version is None
    assert __version__ in str(excinfo.value)
