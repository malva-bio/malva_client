"""Tests for per-cell metadata on stream-derived AnnData objects."""

from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from malva_client.models import CellExpressionStream

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _dense(adata: ad.AnnData) -> np.ndarray:
    matrix: Any = adata.X
    return np.asarray(matrix.toarray())


def _frame(rows):
    frame = pd.DataFrame(rows, columns=["sample_id", "cell_id", "gene", "value"])
    frame.insert(2, "store_index", 0)
    return frame


def _frames():
    """Three batches; obs order differs from the sorted order of ``to_anndata``."""
    return [
        _frame([(2, 20, "G2", 0.5), (2, 20, "G1", 1.0)]),
        _frame([(1, 10, "G1", 2.0), (1, 11, "G2", 3.0)]),
        _frame([(1, 11, "G1", 4.0)]),
    ]


def _metadata():
    return pd.DataFrame(
        {
            "sample_id": [1, 1, 2, 9],
            "cell_id": [10, 11, 20, 99],
            "cell_type": ["B", "T", "NK", "ignored"],
            "total_counts": [100, 200, 300, 400],
            "cl_id": [7, 8, 9, 10],
        }
    )


class _FakeClient:
    """Records metadata-service calls and returns a fixed cell table."""

    def __init__(self, frame):
        self.frame = frame
        self.calls = []

    def get_cells_by_metadata(self, **kwargs):
        self.calls.append(kwargs)
        return self.frame.copy()


# ---------------------------------------------------------------------------
# Joins
# ---------------------------------------------------------------------------


class TestCellMetadataJoin:
    def test_all_cells_get_metadata_columns(self):
        adata = CellExpressionStream(_frames(), cell_metadata=_metadata()).to_anndata()
        obs: Any = adata.obs

        assert adata.n_obs == 3
        assert {"cell_type", "total_counts", "cl_id"}.issubset(obs.columns)
        assert obs.columns.tolist().count("sample_id") == 1
        assert obs.columns.tolist().count("cell_id") == 1
        assert obs.loc["1:10", "cell_type"] == "B"
        assert obs.loc["1:11", "cell_type"] == "T"
        assert obs.loc["2:20", "cell_type"] == "NK"
        assert obs.loc["2:20", "total_counts"] == 300
        assert obs.loc["1:10", "cl_id"] == 7
        # The metadata row for a cell outside the stream is ignored.
        assert "9:99" not in obs.index

    def test_partial_coverage_fills_nan(self):
        metadata = _metadata().query("cell_id != 11")
        adata = CellExpressionStream(_frames(), cell_metadata=metadata).to_anndata()
        obs: Any = adata.obs

        assert pd.isna(obs.loc["1:11", "cell_type"])
        assert pd.isna(obs.loc["1:11", "total_counts"])
        assert obs.loc["1:10", "cell_type"] == "B"

    def test_obs_order_is_unchanged(self):
        plain = CellExpressionStream(_frames()).to_anndata()
        enriched = CellExpressionStream(_frames(), cell_metadata=_metadata()).to_anndata()

        assert enriched.obs_names.tolist() == plain.obs_names.tolist()
        np.testing.assert_allclose(_dense(enriched), _dense(plain))

    @pytest.mark.parametrize("as_index", [False, True])
    def test_accepts_frame_or_indexed_frame(self, as_index):
        metadata = _metadata()
        if as_index:
            metadata = metadata.set_index(["sample_id", "cell_id"])
        adata = CellExpressionStream(_frames(), cell_metadata=metadata).to_anndata()
        obs: Any = adata.obs

        assert obs.loc["1:10", "cell_type"] == "B"
        assert obs.loc["2:20", "total_counts"] == 300

    def test_accepts_dict_keyed_by_cell(self):
        metadata = {
            (1, 10): {"cell_type": "B", "total_counts": 100},
            (2, 20): {"cell_type": "NK", "total_counts": 300},
        }
        adata = CellExpressionStream(_frames(), cell_metadata=metadata).to_anndata()
        obs: Any = adata.obs

        assert obs.loc["1:10", "cell_type"] == "B"
        assert obs.loc["2:20", "cell_type"] == "NK"
        assert pd.isna(obs.loc["1:11", "cell_type"])

    def test_enrich_accepts_later_metadata(self):
        stream = CellExpressionStream(_frames()).enrich(cell_metadata=_metadata())
        obs: Any = stream.to_anndata().obs

        assert obs.loc["1:10", "cell_type"] == "B"
        assert stream.cell_metadata.shape[0] == 4

    def test_opt_out_skips_the_join(self):
        stream = CellExpressionStream(
            _frames(),
            cell_metadata=_metadata(),
            include_cell_metadata=False,
        )
        obs: Any = stream.to_anndata().obs

        assert not {"cell_type", "total_counts", "cl_id"} & set(obs.columns)


# ---------------------------------------------------------------------------
# Fetching
# ---------------------------------------------------------------------------


class TestCellMetadataFetch:
    def test_enrich_fetches_once_and_restricts_to_requested_cells(self):
        client = _FakeClient(_metadata())
        stream = CellExpressionStream(
            _frames(),
            client=client,
            cells=[(1, 10), (1, 11), (2, 20)],
        ).enrich()

        assert len(client.calls) == 1
        assert client.calls[0]["sample_ids"] == [1, 2]
        assert client.calls[0]["include_cell_metadata"] is True
        assert stream.cell_metadata["cell_id"].tolist() == [10, 11, 20]
        obs: Any = stream.to_anndata().obs
        assert obs.loc["1:11", "cell_type"] == "T"

    def test_enrich_does_not_fetch_when_metadata_was_supplied(self):
        client = _FakeClient(_metadata())
        stream = CellExpressionStream(
            _frames(),
            client=client,
            cells=[(1, 10)],
            cell_metadata=_metadata(),
        ).enrich()

        assert client.calls == []

    def test_opt_out_skips_the_fetch(self):
        client = _FakeClient(_metadata())
        stream = CellExpressionStream(
            _frames(),
            client=client,
            cells=[(1, 10), (2, 20)],
            include_cell_metadata=False,
        ).enrich()

        assert client.calls == []
        assert stream.cell_metadata.empty
        obs: Any = stream.to_anndata().obs
        assert "cell_type" not in obs.columns

    def test_client_stream_forwards_cells_and_metadata(self, mock_client, monkeypatch):
        monkeypatch.setattr(
            mock_client,
            "_iter_cell_expression_batches",
            lambda **kwargs: iter(_frames()),
        )
        client = _FakeClient(_metadata())
        monkeypatch.setattr(mock_client, "get_cells_by_metadata", client.get_cells_by_metadata)

        stream = mock_client.stream_cell_expression(
            cells=[{"sample_id": 1, "cell_id": 10}],
        ).enrich()

        assert client.calls[0]["sample_ids"] == [1]
        assert stream.cell_metadata["cell_id"].tolist() == [10]
        obs: Any = stream.to_anndata().obs
        assert obs.loc["1:10", "cell_type"] == "B"


# ---------------------------------------------------------------------------
# Every conversion path
# ---------------------------------------------------------------------------


class TestCellMetadataEverywhere:
    def test_chunked_iter_anndata_carries_the_same_columns(self):
        chunks = list(
            CellExpressionStream(_frames(), cell_metadata=_metadata()).iter_anndata(chunk_obs=2)
        )

        assert len(chunks) == 2
        expected = ["sample_id", "cell_id", "cell_type", "total_counts", "cl_id"]
        for chunk in chunks:
            assert chunk.obs.columns.tolist() == expected
        first: Any = chunks[0].obs
        second: Any = chunks[1].obs
        assert first.loc["1:10", "cell_type"] == "B"
        assert second.loc["1:11", "cell_type"] == "T"

    def test_per_batch_iter_anndata_carries_metadata(self):
        parts = list(CellExpressionStream(_frames(), cell_metadata=_metadata()).iter_anndata())

        assert all("cell_type" in part.obs.columns for part in parts)
        obs: Any = parts[0].obs
        assert obs.loc["2:20", "cell_type"] == "NK"

    def test_write_h5ad_carries_metadata(self, tmp_path):
        path = tmp_path / "cells.h5ad"
        CellExpressionStream(_frames(), cell_metadata=_metadata()).write_h5ad(path)
        adata = ad.read_h5ad(path)
        obs: Any = adata.obs

        assert obs.loc["1:10", "cell_type"] == "B"
        assert obs.loc["2:20", "total_counts"] == 300

    def test_write_h5ad_opt_out(self, tmp_path):
        path = tmp_path / "cells.h5ad"
        CellExpressionStream(
            _frames(),
            cell_metadata=_metadata(),
            include_cell_metadata=False,
        ).write_h5ad(path)
        adata = ad.read_h5ad(path)

        assert "cell_type" not in adata.obs.columns
