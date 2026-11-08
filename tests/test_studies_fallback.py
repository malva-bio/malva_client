"""Hierarchy fallbacks for the public provider's study routes."""

from unittest.mock import patch

import pandas as pd
import pytest
import responses

from malva_client.client import MalvaClient
from malva_client.exceptions import MalvaAPIError

BASE = "https://test.malva.bio"

HIERARCHY = {
    "datasets": [
        {
            "dataset_id": "malva-public-index",
            "name": "Public index",
            "studies": [
                {
                    "study_name": "Study A",
                    "study_title": "Study A title",
                    "laboratory": "Lab A",
                    "sample_count": 2,
                    "total_cells": 20,
                    "organisms": ["Homo sapiens"],
                    "organs": ["brain"],
                    "diseases": ["normal"],
                },
                {
                    "study_name": "Study B",
                    "study_title": "Study B title",
                    "laboratory": "Lab B",
                    "sample_count": 1,
                    "total_cells": 10,
                    "organisms": ["Homo sapiens"],
                    "organs": ["liver"],
                    "diseases": ["normal"],
                },
            ],
        }
    ]
}


def _client(**kwargs):
    with patch.object(MalvaClient, "_test_connection"):
        client = MalvaClient(base_url=BASE, api_token="tok", **kwargs)
    client._endpoint_style = "legacy"
    return client


@responses.activate
def test_get_studies_falls_back_to_hierarchy():
    client = _client()
    responses.add(
        responses.GET,
        f"{BASE}/samples/api/studies",
        json={"studies": [], "total_count": 0},
    )
    responses.add(
        responses.GET,
        f"{BASE}/samples/api/datasets/hierarchy",
        json=HIERARCHY,
    )

    frame = client.get_studies()

    assert isinstance(frame, pd.DataFrame)
    assert len(frame) == 2
    assert set(["study_name", "study_title", "sample_count", "dataset_id"]).issubset(frame.columns)
    assert frame["dataset_id"].tolist() == ["malva-public-index", "malva-public-index"]
    assert frame["study_name"].tolist() == ["Study A", "Study B"]


@responses.activate
def test_get_dataset_studies_falls_back_to_hierarchy():
    client = _client()
    responses.add(
        responses.GET,
        f"{BASE}/samples/api/datasets/malva-public-index/studies",
        json={
            "dataset_id": "malva-public-index",
            "studies": [],
            "total_count": 0,
            "page": 1,
            "page_size": 1,
            "has_more": False,
        },
    )
    responses.add(
        responses.GET,
        f"{BASE}/samples/api/datasets/hierarchy",
        json=HIERARCHY,
    )

    response = client.get_dataset_studies("malva-public-index", page=1, page_size=1)

    assert response["total_count"] == 2
    assert response["page"] == 1
    assert response["page_size"] == 1
    assert response["has_more"] is True
    assert response["source"] == "dataset_hierarchy"
    assert [study["study_name"] for study in response["studies"]] == ["Study A"]


@responses.activate
def test_get_study_samples_reports_empty_route_for_known_study():
    client = _client()
    responses.add(
        responses.GET,
        f"{BASE}/samples/api/datasets/malva-public-index/studies/Study%20A/samples",
        json={
            "dataset_id": "malva-public-index",
            "study_name": "Study A",
            "samples": [],
            "total_count": 0,
            "page": 1,
            "page_size": 50,
            "has_more": False,
        },
    )
    responses.add(
        responses.GET,
        f"{BASE}/samples/api/datasets/hierarchy",
        json=HIERARCHY,
    )

    with pytest.raises(MalvaAPIError, match="requires a server-side fix"):
        client.get_study_samples("malva-public-index", "Study A")


@responses.activate
def test_get_study_samples_keeps_empty_result_for_unknown_study():
    client = _client()
    responses.add(
        responses.GET,
        f"{BASE}/samples/api/datasets/malva-public-index/studies/Missing/samples",
        json={
            "dataset_id": "malva-public-index",
            "study_name": "Missing",
            "samples": [],
            "total_count": 0,
            "page": 1,
            "page_size": 50,
            "has_more": False,
        },
    )
    responses.add(
        responses.GET,
        f"{BASE}/samples/api/datasets/hierarchy",
        json=HIERARCHY,
    )

    response = client.get_study_samples("malva-public-index", "Missing")

    assert response["samples"] == []
    assert response["total_count"] == 0
