"""Sample-details errors must expose the server body and provider context."""

import json
from urllib.parse import parse_qs, urlparse

import pytest
import responses

from malva_client.exceptions import MalvaAPIError

BASE = "https://test.malva.bio"


@responses.activate
def test_get_sample_details_surfaces_server_error_body(mock_client):
    def fail(request):
        assert parse_qs(urlparse(request.url).query)["provider"] == ["malva-public-index"]
        return 500, {}, json.dumps({"detail": "Failed to retrieve sample metadata"})

    responses.add_callback(
        responses.GET,
        f"{BASE}/samples/api/sample/bad-uuid",
        callback=fail,
    )

    with pytest.raises(
        MalvaAPIError,
        match="Failed to retrieve metadata for sample 'bad-uuid'.*Failed to retrieve sample metadata",
    ):
        mock_client.get_sample_details("bad-uuid", provider="malva-public-index")
