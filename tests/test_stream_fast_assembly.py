"""Fast compact assembly must match the long-format DataFrame path exactly."""

from __future__ import annotations

import struct
import zlib
from typing import Any

import anndata as ad
import msgpack
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from malva_client import models, stream_codec
from malva_client.models import CellExpressionStream, _arrays_to_anndata

# ---------------------------------------------------------------------------
# Compact stream fixtures
# ---------------------------------------------------------------------------


def _leb128_encode(values) -> bytes:
    out = bytearray()
    for value in values:
        while True:
            byte = value & 0x7F
            value >>= 7
            out.append(byte | (0x80 if value else 0))
            if not value:
                break
    return bytes(out)


def _data_frame(store_index, requested, present, rows):
    mask = np.zeros(len(requested), dtype=bool)
    mask[list(present)] = True
    counts = [len(row) for row in rows]
    positions = [position for row in rows for position, _value in row]
    values = [value for row in rows for _position, value in row]
    return {
        "type": "data",
        "store_index": store_index,
        "n": len(positions),
        "cells": np.packbits(mask).tobytes(),
        "nnz": _leb128_encode(counts),
        "feature_width": 2,
        "feature_mode": "absolute",
        "features": np.asarray(positions, dtype="<u2").tobytes(),
        "values": np.asarray(values, dtype="<f4").tobytes(),
    }


def _frames():
    requested = [(7, 10), (7, 11), (7, 12), (8, 20), (8, 21)]
    return requested, [
        {
            "type": "start",
            "n_cells": len(requested),
            "layout": "cells_compact_v1",
            "compression": "zlib",
        },
        {"type": "features", "store_index": 0, "features": ["G2", "G1", "G0"]},
        # Positions are deliberately unsorted and repeat G1 inside the first cell.
        _data_frame(
            0,
            requested,
            [0, 2, 3],
            [
                [(2, 1.5), (0, 2.5), (1, 0.5), (1, 0.25)],
                [(1, 4.0), (2, 5.0)],
                [(0, 6.0)],
            ],
        ),
        {"type": "features", "store_index": 1, "features": ["G1", "G3"]},
        _data_frame(
            1,
            requested,
            [1, 4],
            [
                [(0, 7.0)],
                [(0, 10.0)],
            ],
        ),
        {"type": "done", "n_frames": 2, "n_values": 11},
    ]


def _stream_payload(frames, encoding="zlib"):
    out = bytearray()
    for index, frame in enumerate(frames):
        body = msgpack.packb(frame, use_bin_type=True)
        if index and encoding:
            body = zlib.compress(body)
        out += struct.pack(">Q", len(body)) + body
    return bytes(out)


class _FakeResponse:
    def __init__(self, payload: bytes):
        self._payload = payload
        self.status_code = 200
        self.ok = True
        self.text = ""

    def iter_content(self, chunk_size: int):
        for start in range(0, len(self._payload), chunk_size):
            yield self._payload[start : start + chunk_size]

    def close(self):
        return None


@pytest.fixture
def compact_client(mock_client, monkeypatch):
    payload = _stream_payload(_frames()[1])
    monkeypatch.setattr(
        mock_client.session,
        "post",
        lambda *args, **kwargs: _FakeResponse(payload),
    )
    return mock_client


def _reference_frames(compact_client, requested):
    return list(
        compact_client._iter_cell_expression_batches(
            cells=requested,
            dataset_id="ds1",
            layout="compact",
            compress="zlib",
            fast=False,
        )
    )


def _reference_anndata(frames, uns):
    return _arrays_to_anndata(
        np.concatenate([frame["sample_id"].to_numpy(copy=True) for frame in frames]),
        np.concatenate([frame["cell_id"].to_numpy(copy=True) for frame in frames]),
        np.concatenate([frame["gene"].astype(str).to_numpy(copy=True) for frame in frames]),
        np.concatenate(
            [
                pd.to_numeric(frame["value"], errors="coerce")
                .fillna(0.0)
                .to_numpy(dtype=np.float32)
                for frame in frames
            ]
        ),
        uns=uns,
    )


def _assert_anndata_equal(actual: ad.AnnData, expected: ad.AnnData) -> None:
    actual_matrix: Any = actual.X
    expected_matrix: Any = expected.X
    actual_obs: Any = actual.obs
    expected_obs: Any = expected.obs
    actual_var: Any = actual.var
    expected_var: Any = expected.var
    assert actual.shape == expected.shape
    assert sparse.isspmatrix_csr(actual_matrix)
    assert actual_matrix.dtype == expected_matrix.dtype == np.float32
    assert actual.obs_names.tolist() == expected.obs_names.tolist()
    assert actual.var_names.tolist() == expected.var_names.tolist()
    pd.testing.assert_frame_equal(actual_obs, expected_obs)
    pd.testing.assert_frame_equal(actual_var, expected_var)
    np.testing.assert_allclose(
        np.asarray(actual_matrix.todense()), np.asarray(expected_matrix.todense())
    )
    assert actual.uns["malva"] == expected.uns["malva"]


# ---------------------------------------------------------------------------
# Public stream and to_anndata
# ---------------------------------------------------------------------------


class TestFastCompactStream:
    def test_public_iteration_still_yields_the_long_format_frames(self, compact_client):
        requested = _frames()[0]
        expected = _reference_frames(compact_client, requested)
        actual = list(
            compact_client.stream_cell_expression(
                cells=requested, dataset_id="ds1", layout="compact", compress="zlib"
            )
        )

        assert len(actual) == len(expected) == 2
        for got, want in zip(actual, expected):
            pd.testing.assert_frame_equal(got, want)

    def test_to_anndata_matches_the_previous_implementation(self, compact_client):
        requested = _frames()[0]
        frames = _reference_frames(compact_client, requested)
        uns = {
            "source": "cell_expression_stream",
            "dataset_id": "ds1",
            "sample_ids": [7, 8],
            "gene_ids": None,
            "layout": "cells_by_genes",
        }
        expected = _reference_anndata(frames, uns)
        actual = compact_client.stream_cell_expression(
            cells=requested,
            dataset_id="ds1",
            layout="compact",
            compress="zlib",
            include_cell_metadata=False,
        ).to_anndata()

        _assert_anndata_equal(actual, expected)
        # The duplicate G1 entry inside a cell is summed, as before.
        actual_matrix: Any = actual.X
        expected_matrix: Any = expected.X
        assert actual_matrix.nnz == expected_matrix.nnz

    def test_raw_batches_carry_codes_not_names(self, compact_client):
        requested = _frames()[0]
        raw = list(
            compact_client.stream_cell_expression(
                cells=requested, dataset_id="ds1", layout="compact", compress="zlib"
            )._iter_raw()
        )

        assert all(isinstance(batch, stream_codec.ExpressionBatch) for batch in raw)
        assert all(batch.gene_codes.dtype == np.int64 for batch in raw)
        assert all(batch.values.dtype == np.float32 for batch in raw)
        # Two stores share G1; the second store maps it to the same gene name.
        names = [str(name) for batch in raw for name in batch.gene_names]
        assert "G1" in names and "G3" in names

    def test_default_iter_anndata_matches_the_dataframe_path(self, compact_client):
        requested = _frames()[0]
        frames = _reference_frames(compact_client, requested)
        expected = list(
            CellExpressionStream(frames, dataset_id="ds1", sample_ids=[7, 8]).iter_anndata()
        )
        actual = list(
            compact_client.stream_cell_expression(
                cells=requested,
                dataset_id="ds1",
                layout="compact",
                compress="zlib",
                include_cell_metadata=False,
            ).iter_anndata()
        )

        assert len(actual) == len(expected) == 2
        for got, want in zip(actual, expected):
            _assert_anndata_equal(got, want)
            assert got.uns["malva"]["batch"] is True

    def test_chunked_iter_anndata_matches_the_dataframe_path(self, compact_client):
        requested = _frames()[0]
        frames = _reference_frames(compact_client, requested)
        expected = list(
            CellExpressionStream(frames, dataset_id="ds1", sample_ids=[7, 8]).iter_anndata(
                chunk_obs=3
            )
        )
        actual = list(
            compact_client.stream_cell_expression(
                cells=requested,
                dataset_id="ds1",
                layout="compact",
                compress="zlib",
                include_cell_metadata=False,
            ).iter_anndata(chunk_obs=3)
        )

        assert [chunk.n_obs for chunk in actual] == [chunk.n_obs for chunk in expected]
        for got, want in zip(actual, expected):
            _assert_anndata_equal(got, want)
            assert got.uns["malva"]["chunk"] is True

    def test_chunked_output_concatenates_to_the_full_matrix(self, compact_client):
        requested = _frames()[0]
        full = compact_client.stream_cell_expression(
            cells=requested,
            dataset_id="ds1",
            layout="compact",
            compress="zlib",
            include_cell_metadata=False,
        ).to_anndata()
        chunks = list(
            compact_client.stream_cell_expression(
                cells=requested,
                dataset_id="ds1",
                layout="compact",
                compress="zlib",
                include_cell_metadata=False,
            ).iter_anndata(chunk_obs=3)
        )
        merged = ad.concat(chunks, join="outer", fill_value=0)
        merged = merged[full.obs_names, full.var_names].copy()

        merged_matrix: Any = merged.X
        full_matrix: Any = full.X
        assert merged.shape == full.shape
        assert merged.obs_names.tolist() == full.obs_names.tolist()
        assert merged.var_names.tolist() == full.var_names.tolist()
        np.testing.assert_allclose(
            np.asarray(merged_matrix.todense()), np.asarray(full_matrix.todense())
        )


# ---------------------------------------------------------------------------
# Observation axis fallbacks
# ---------------------------------------------------------------------------


class TestObservationAxisFallbacks:
    def _frame(self, samples, cells, genes, values):
        frame = pd.DataFrame(
            {
                "sample_id": samples,
                "cell_id": cells,
                "store_index": np.zeros(len(samples), dtype=np.int32),
                "gene": genes,
                "value": values,
            }
        )
        return frame

    @pytest.mark.parametrize(
        "samples,cells",
        [
            (["s2", "s2", "s1"], ["c20", "c20", "c10"]),
            ([2.0, 2.0, 1.0], [20.0, 20.0, 10.0]),
            ([2, 2, 1], [1 << 33, 1 << 33, 10]),
        ],
    )
    def test_non_integer_keys_match_the_previous_implementation(self, samples, cells):
        frame = self._frame(samples, cells, ["G2", "G1", "G1"], [0.5, 1.0, 2.0])
        actual = CellExpressionStream([frame]).to_anndata()
        expected = _reference_anndata([frame], actual.uns["malva"])

        _assert_anndata_equal(actual, expected)

    def test_integer_to_object_transition_keeps_rows_and_order(self):
        first = self._frame([2, 2], [20, 20], ["G2", "G1"], [0.5, 1.0])
        second = self._frame(["s1"], ["c10"], ["G2"], [0.7])
        actual = CellExpressionStream([first, second]).to_anndata()
        expected = _reference_anndata([first, second], actual.uns["malva"])

        _assert_anndata_equal(actual, expected)
        assert actual.obs["sample_id"].tolist() == ["2", "s1"]


# ---------------------------------------------------------------------------
# Accumulator behaviour
# ---------------------------------------------------------------------------


def test_accumulator_add_has_no_python_row_loop():
    import ast
    import inspect
    import textwrap

    tree = ast.parse(textwrap.dedent(inspect.getsource(models._AnnDataAccumulator.add)))
    loops = (
        ast.For,
        ast.AsyncFor,
        ast.While,
        ast.ListComp,
        ast.SetComp,
        ast.DictComp,
        ast.GeneratorExp,
    )
    assert not [node for node in ast.walk(tree) if isinstance(node, loops)]


def test_expression_batch_is_accepted_directly():
    requested = [(1, 10), (1, 11)]
    frame = _data_frame(0, requested, [0, 1], [[(0, 1.5), (1, 2.5)], [(1, 3.5)]])
    stream = CellExpressionStream(
        [
            pd.DataFrame(
                {
                    "sample_id": [1, 1, 1],
                    "cell_id": [10, 10, 11],
                    "store_index": [0, 0, 0],
                    "gene": ["G0", "G1", "G1"],
                    "value": [1.5, 2.5, 3.5],
                }
            )
        ]
    )
    expected = stream.to_anndata()

    batch = stream_codec.ExpressionBatch(
        sample_ids=np.array([1, 1, 1], dtype=np.int64),
        cell_ids=np.array([10, 10, 11], dtype=np.int64),
        gene_codes=np.array([0, 1, 1], dtype=np.int64),
        values=np.array([1.5, 2.5, 3.5], dtype=np.float32),
        store_index=0,
        gene_names=["G0", "G1"],
        new_codes=np.array([0, 1], dtype=np.int64),
    )
    accumulator = models._AnnDataAccumulator(uns=expected.uns["malva"])
    accumulator.add(batch)
    actual = accumulator.finish(sort=True)

    _assert_anndata_equal(actual, expected)


def test_accumulator_rejects_unknown_feature_codes():
    batch = stream_codec.ExpressionBatch(
        sample_ids=np.array([1], dtype=np.int64),
        cell_ids=np.array([10], dtype=np.int64),
        gene_codes=np.array([5], dtype=np.int64),
        values=np.array([1.0], dtype=np.float32),
        store_index=0,
        gene_names=["G0"],
        new_codes=np.empty(0, dtype=np.int64),
    )
    with pytest.raises(ValueError, match="unknown feature code"):
        models._AnnDataAccumulator().add(batch)
