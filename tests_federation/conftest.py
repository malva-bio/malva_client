"""Shared fixtures for the opt-in federation suite.

No Explorer imports live here: the suite exercises the Explorer over HTTP
through the legacy ``malva_client`` only.
"""

from unittest.mock import patch

import pytest

from malva_client.client import MalvaClient

FEDERATION_BASE_URL = "https://federation.test"


@pytest.fixture
def client():
    """A legacy client whose connection probe is patched out.

    The suite talks to mocked Explorer endpoints through ``client._request``;
    ``?provider=<id>`` is the routing contract for the metadata proxies.
    """
    with patch.object(MalvaClient, "_test_connection"):
        client = MalvaClient(
            base_url=FEDERATION_BASE_URL,
            api_token="federation-test-token",
        )
    # This suite pins the legacy Explorer routes directly.
    client._endpoint_style = "legacy"
    return client


@pytest.fixture
def base_url():
    return FEDERATION_BASE_URL
