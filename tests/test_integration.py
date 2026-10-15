"""Integration tests against the real Malva platform.

These exercise the configured endpoints end to end. They are skipped unless
``MALVA_API_TOKEN`` is set, and marked ``integration`` so the default test run
(the offline suite) stays fast and network-free.

Local use::

    export MALVA_API_URL="https://malva.mdc-berlin.de"
    export MALVA_API_TOKEN="..."
    pytest -m integration

In CI the token comes from the ``MALVA_API_TOKEN`` repository secret and the
URL from the ``MALVA_API_URL`` repository variable.
"""

import os

import pytest

from malva_client.client import MalvaClient
from malva_client.models import CoverageResult, SearchResult

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(
        not os.environ.get("MALVA_API_TOKEN"),
        reason="MALVA_API_TOKEN is not set",
    ),
]

URL = os.environ.get("MALVA_API_URL", "https://malva.mdc-berlin.de")

# A short real sequence, and a second one to exercise the list/batch path.
SEQUENCE = (
    "AACTTATAAGTGTTGTGAACCCCCACCCAGGCAGTATCCATGAAAGCACAAGTGACTAGTCCTAT"
    "GATGTACAAAGCCTGTATCTCTGTGATGATTTCTGTGC"
)
SEQUENCE_2 = "GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA"
GENE = "FOXP3"
REGION = ("chr12", 52487695, 52491695)


@pytest.fixture(scope="module")
def client():
    return MalvaClient(
        base_url=URL,
        api_token=os.environ["MALVA_API_TOKEN"],
        timeout=180,
    )


@pytest.fixture(scope="module")
def dataset_id(client):
    """A catalog dataset id that also has a coexpression index (and UMAP).

    The dataset ids come from the metadata hierarchy; the coexpression service
    keys its indices by the same id, but not every dataset has an atlas loaded,
    so probe the UMAP endpoint and use the first dataset that answers.
    """
    try:
        hierarchy = client.get_datasets_hierarchy()
    except Exception:
        hierarchy = {}
    ids = [
        entry.get("dataset_id")
        for entry in (hierarchy.get("datasets") or [])
        if entry.get("dataset_id")
    ]
    for candidate in ids:
        try:
            client.get_umap_coordinates(candidate)
            return candidate
        except Exception:
            continue
    pytest.skip("no dataset with a coexpression index / UMAP is available")


@pytest.fixture(scope="module")
def gene_search(client):
    """One gene search reused by the downstream tests (to spare quota)."""
    return client.search(GENE, poll_interval=2, max_wait=240)


class TestConnection:
    def test_is_authenticated(self, client):
        assert client.is_authenticated() is True

    def test_quota_status(self, client):
        assert isinstance(client.get_quota_status(), dict)


class TestSearch:
    def test_search_gene(self, gene_search):
        assert isinstance(gene_search, SearchResult)
        assert gene_search.job_id
        assert gene_search.df is not None

    def test_search_genes_list(self, client):
        result = client.search_genes([GENE], poll_interval=2, max_wait=240)
        assert isinstance(result, SearchResult)

    def test_search_single_sequence(self, client):
        result = client.search_sequences(SEQUENCE, poll_interval=2, max_wait=240)
        assert isinstance(result, SearchResult)

    def test_search_sequence_list(self, client):
        result = client.search_sequences([SEQUENCE, SEQUENCE_2], poll_interval=2, max_wait=240)
        assert isinstance(result, SearchResult)


class TestCoverage:
    def test_region_coverage(self, client):
        chrom, start, end = REGION
        result = client.get_coverage(chrom, start, end, poll_interval=2, max_wait=300)
        assert isinstance(result, CoverageResult)

    def test_sequence_coverage(self, client):
        result = client.get_sequence_coverage(SEQUENCE, poll_interval=2, max_wait=300)
        assert isinstance(result, CoverageResult)


class TestMetadata:
    def test_datasets_hierarchy(self, client):
        assert isinstance(client.get_datasets_hierarchy(), dict)

    def test_overview_stats(self, client):
        assert isinstance(client.get_overview_stats(), dict)

    def test_database_stats(self, client):
        assert isinstance(client.get_database_stats(), dict)

    def test_available_filters(self, client):
        assert client.get_available_filters() is not None

    def test_samples(self, client):
        assert client.get_samples() is not None

    def test_search_samples(self, client):
        assert client.search_samples(query="brain") is not None

    def test_studies(self, client):
        assert client.get_studies() is not None


class TestCellLevel:
    def test_retrieve_cells(self, client, gene_search):
        result = client.retrieve_cells(gene_search, features=GENE)
        assert result is not None

    def test_umap_coordinates(self, client, dataset_id):
        assert client.get_umap_coordinates(dataset_id) is not None


class TestCoexpression:
    def test_coexpression_from_job(self, client, gene_search, dataset_id):
        result = client.get_coexpression(
            gene_search.job_id, dataset_id, poll_interval=2, max_wait=300
        )
        assert result is not None
