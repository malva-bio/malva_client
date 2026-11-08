"""Values-export timeouts must fall back to the other cell fetchers."""

import logging

import requests
import responses

BASE = "https://test.malva.bio"


def _register_full_search():
    responses.add(
        responses.GET,
        f"{BASE}/search/j1",
        json={
            "status": "completed",
            "job_id": "j1",
            "results": {
                "BRCA1": {
                    "cell": [101],
                    "sample": [10],
                    "expression": [3.0],
                }
            },
        },
    )


@responses.activate
def test_values_export_timeout_falls_back_to_full_search(mock_client, caplog):
    caplog.set_level(logging.WARNING)
    responses.add(
        responses.POST,
        f"{BASE}/exports/cell-expression",
        json={"export_id": "e1", "status": "pending"},
    )
    _register_full_search()

    result = mock_client.retrieve_cells(
        "j1",
        features="BRCA1",
        include_sample_metadata=False,
        poll_interval=0,
        max_wait=0,
    )

    assert result.source == "search_full"
    assert result.to_dataframe()["value"].tolist() == [3.0]
    assert "values export timed out" in caplog.text


@responses.activate
def test_values_export_transport_timeout_falls_back_to_full_search(mock_client, caplog):
    caplog.set_level(logging.WARNING)
    responses.add(
        responses.POST,
        f"{BASE}/exports/cell-expression",
        json={"export_id": "e1", "status": "pending"},
    )

    def timeout_status(request):
        raise requests.exceptions.Timeout("status request timed out")

    responses.add_callback(
        responses.GET,
        f"{BASE}/exports/cell-expression/e1",
        callback=timeout_status,
    )
    _register_full_search()

    result = mock_client.retrieve_cells(
        "j1",
        features="BRCA1",
        include_sample_metadata=False,
        poll_interval=0,
        max_wait=10,
    )

    assert result.source == "search_full"
    assert result.to_dataframe()["value"].tolist() == [3.0]
    assert "timed out" in caplog.text.lower()
