"""Tests for the progressive cell-expression accumulator."""

import ast
import inspect
import textwrap
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from malva_client import models
from malva_client.models import (
    CellExpressionStream,
    _AnnDataAccumulator,
    _arrays_to_anndata,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _dense(adata: ad.AnnData) -> np.ndarray:
    matrix: Any = adata.X
    return np.asarray(matrix.toarray())


def _frame(rows):
    """Build a server-shaped long-format frame from (sample, cell, gene, value)."""
    frame = pd.DataFrame(rows, columns=["sample_id", "cell_id", "gene", "value"])
    frame.insert(2, "store_index", 0)
    return frame


def _legacy_to_anndata(frames, *, sample_metadata=None, uns=None):
    """Reproduce the pre-0.8 implementation: concatenate, then build once."""
    sample_parts = []
    cell_parts = []
    gene_parts = []
    value_parts = []
    for batch in frames:
        if batch is None or batch.empty:
            continue
        sample_parts.append(batch["sample_id"].to_numpy(copy=True))
        cell_parts.append(batch["cell_id"].to_numpy(copy=True))
        gene_parts.append(batch["gene"].astype(str).to_numpy(copy=True))
        value_parts.append(
            pd.to_numeric(batch["value"], errors="coerce").fillna(0.0).to_numpy(dtype=np.float32)
        )
    if not sample_parts:
        return models._empty_anndata()
    return _arrays_to_anndata(
        np.concatenate(sample_parts),
        np.concatenate(cell_parts),
        np.concatenate(gene_parts),
        np.concatenate(value_parts),
        sample_metadata=sample_metadata,
        uns=uns,
    )


def _tricky_frames():
    """Out-of-order cells, duplicate entries, and a cell split across batches."""
    return [
        _frame(
            [
                (2, 20, "G2", 0.5),
                (2, 20, "G1", 1.0),
                (2, 20, "G1", 0.25),
                (1, 11, "G3", 2.0),
                (1, 10, "G1", 3.0),
            ]
        ),
        _frame(
            [
                (1, 10, "G2", 4.0),  # same cell continues in the next batch
                (1, 10, "G3", np.nan),
                (3, 30, "G4", 5.0),
            ]
        ),
        _frame([(3, 30, "G1", 6.0), (3, 31, "G2", 7.0)]),
    ]


def _synthetic_frames(n_batches=4, cells_per_batch=100, n_genes=20):
    frames = []
    for batch_index in range(n_batches):
        cells = batch_index * cells_per_batch + np.arange(cells_per_batch)
        rows = np.repeat(cells, n_genes)
        genes = [f"G{gene:03d}" for gene in range(n_genes)] * cells_per_batch
        values = np.arange(rows.size, dtype=np.float32) % 7 + 0.5
        frames.append(
            pd.DataFrame(
                {
                    "sample_id": np.ones(rows.size, dtype=np.int64),
                    # Zero-padded ids keep first-seen and lexicographic order equal,
                    # so chunked output can be compared with sorted per-batch parts.
                    "cell_id": [f"c{row:05d}" for row in rows],
                    "store_index": np.zeros(rows.size, dtype=np.int32),
                    "gene": genes,
                    "value": values,
                }
            )
        )
    return frames


# ---------------------------------------------------------------------------
# to_anndata parity with the previous implementation
# ---------------------------------------------------------------------------


class TestToAnndataParity:
    def test_matches_legacy_arrays_and_order(self):
        frames = _tricky_frames()
        metadata = pd.DataFrame(
            {"organ": ["brain", "liver", "kidney"]},
            index=pd.Index([1, 2, 3], name="sample_id"),
        )
        stream = CellExpressionStream(
            frames,
            dataset_id="ds1",
            sample_ids=[1, 2, 3],
            gene_ids=["G1", "G2", "G3", "G4"],
            sample_metadata=metadata,
        )
        legacy_uns = stream._uns()

        adata = stream.to_anndata()
        expected = _legacy_to_anndata(
            frames,
            sample_metadata=metadata,
            uns=legacy_uns,
        )

        assert isinstance(adata.X, sparse.csr_matrix)
        assert adata.X.dtype == np.float32
        assert adata.shape == expected.shape
        assert adata.obs_names.tolist() == expected.obs_names.tolist()
        assert adata.var_names.tolist() == expected.var_names.tolist()
        obs: Any = adata.obs
        expected_obs: Any = expected.obs
        pd.testing.assert_frame_equal(obs, expected_obs)
        np.testing.assert_allclose(_dense(adata), _dense(expected))
        assert adata.uns["malva"] == expected.uns["malva"]

    def test_keeps_sorted_axes(self):
        adata = CellExpressionStream(_tricky_frames()).to_anndata()

        assert adata.obs_names.tolist() == sorted(adata.obs_names)
        assert adata.var_names.tolist() == sorted(adata.var_names)

    def test_single_batch(self):
        frame = _frame([(1, 10, "G1", 1.0), (1, 10, "G2", 2.0)])
        adata = CellExpressionStream([frame]).to_anndata()

        assert adata.shape == (1, 2)
        assert adata.obs_names.tolist() == ["1:10"]
        np.testing.assert_allclose(_dense(adata), [[1.0, 2.0]])

    def test_empty_stream_and_empty_batches(self):
        assert CellExpressionStream([]).to_anndata().shape == (0, 0)
        assert CellExpressionStream([None, pd.DataFrame()]).to_anndata().shape == (0, 0)

    def test_split_cell_is_one_observation(self):
        frames = [
            _frame([(1, 10, "G1", 1.0)]),
            _frame([(1, 10, "G2", 2.0)]),
        ]
        adata = CellExpressionStream(frames).to_anndata()

        assert adata.n_obs == 1
        assert adata.var_names.tolist() == ["G1", "G2"]
        np.testing.assert_allclose(_dense(adata), [[1.0, 2.0]])

    def test_missing_columns_are_rejected(self):
        frame = pd.DataFrame({"sample_id": [1], "cell_id": [10], "gene": ["G1"]})
        with pytest.raises(ValueError, match="must contain sample_id"):
            CellExpressionStream([frame]).to_anndata()


# ---------------------------------------------------------------------------
# Accumulator state
# ---------------------------------------------------------------------------


class TestAccumulator:
    def test_progress_properties_follow_ingestion(self):
        accumulator = _AnnDataAccumulator()
        assert (accumulator.n_obs, accumulator.n_vars, accumulator.bytes_estimate) == (0, 0, 0)

        accumulator.add(_frame([(1, 10, "G1", 1.0), (1, 11, "G2", 2.0)]))
        assert accumulator.n_obs == 2
        assert accumulator.n_vars == 2
        assert accumulator.bytes_estimate > 0

        accumulator.add(_frame([(1, 11, "G1", 3.0), (2, 20, "G3", 4.0)]))
        assert accumulator.n_obs == 3
        assert accumulator.n_vars == 3

    def test_flush_keeps_the_variable_axis(self):
        accumulator = _AnnDataAccumulator()
        accumulator.add(_frame([(1, 10, "G2", 1.0)]))
        first = accumulator.flush()

        assert accumulator.n_obs == 0
        assert accumulator.n_vars == 1
        accumulator.add(_frame([(2, 20, "G1", 2.0)]))
        second = accumulator.flush()

        assert first.var_names.tolist() == ["G2"]
        assert second.var_names.tolist() == ["G2", "G1"]
        np.testing.assert_allclose(_dense(second), [[0.0, 2.0]])

    def test_add_has_no_python_row_loop(self):
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


# ---------------------------------------------------------------------------
# Chunked iter_anndata
# ---------------------------------------------------------------------------


class TestChunkedIterAnndata:
    def test_chunks_concatenate_to_the_per_batch_result(self):
        frames = _synthetic_frames()
        expected = ad.concat(
            CellExpressionStream(frames).iter_anndata(), join="outer", fill_value=0
        )
        chunks = list(CellExpressionStream(frames).iter_anndata(chunk_obs=300))

        assert len(chunks) == 2
        assert [chunk.n_obs for chunk in chunks] == [300, 100]
        merged = ad.concat(chunks, join="outer", fill_value=0)
        assert merged.obs_names.tolist() == expected.obs_names.tolist()
        assert merged.var_names.tolist() == expected.var_names.tolist()
        matrix: Any = merged.X
        assert matrix.dtype == np.float32
        np.testing.assert_allclose(_dense(merged), _dense(expected))

    def test_chunks_are_nested_on_the_variable_axis(self):
        frames = [
            _frame([(1, 10, "G2", 1.0), (1, 11, "G2", 2.0)]),
            _frame([(1, 12, "G1", 3.0), (2, 20, "G1", 4.0)]),
            _frame([(2, 21, "G3", 5.0)]),
        ]
        chunks = list(CellExpressionStream(frames).iter_anndata(chunk_obs=2))

        assert [chunk.var_names.tolist() for chunk in chunks] == [
            ["G2"],
            ["G2", "G1"],
            ["G2", "G1", "G3"],
        ]
        for previous, current in zip(chunks, chunks[1:]):
            assert current.var_names[: len(previous.var_names)].tolist() == (
                previous.var_names.tolist()
            )
        np.testing.assert_allclose(_dense(chunks[1]), [[0.0, 3.0], [0.0, 4.0]])

    def test_last_chunk_may_be_short(self):
        frames = _synthetic_frames(n_batches=3, cells_per_batch=10, n_genes=2)
        chunks = list(CellExpressionStream(frames).iter_anndata(chunk_obs=12))

        assert [chunk.n_obs for chunk in chunks] == [20, 10]
        assert all(chunk.uns["malva"]["chunk"] is True for chunk in chunks)
        assert all(chunk.uns["malva"]["source"] == "cell_expression_stream" for chunk in chunks)

    def test_chunk_equals_concat_for_an_empty_stream(self):
        assert list(CellExpressionStream([]).iter_anndata(chunk_obs=5)) == []
        assert list(CellExpressionStream([None, pd.DataFrame()]).iter_anndata(5)) == []

    def test_chunked_mode_sums_a_cell_split_inside_a_chunk(self):
        frames = [_frame([(1, 10, "G1", 1.0)]), _frame([(1, 10, "G2", 2.0)])]
        chunks = list(CellExpressionStream(frames).iter_anndata(chunk_obs=10))

        assert len(chunks) == 1
        assert chunks[0].obs_names.tolist() == ["1:10"]
        np.testing.assert_allclose(_dense(chunks[0]), [[1.0, 2.0]])

    @pytest.mark.parametrize("chunk_obs", [0, -3, True, 2.5, "10"])
    def test_invalid_chunk_obs_is_rejected(self, chunk_obs):
        with pytest.raises(ValueError, match="positive integer"):
            list(CellExpressionStream(_tricky_frames()).iter_anndata(chunk_obs=chunk_obs))

    def test_default_mode_is_still_one_part_per_batch(self):
        frames = _tricky_frames()
        parts = list(CellExpressionStream(frames).iter_anndata())

        assert len(parts) == len(frames)
        assert all(part.uns["malva"]["batch"] is True for part in parts)
        assert all(part.var_names.tolist() == sorted(part.var_names) for part in parts)


# ---------------------------------------------------------------------------
# Memory shape
# ---------------------------------------------------------------------------


def test_to_anndata_folds_batches_without_a_legacy_concatenation(monkeypatch):
    frames = _synthetic_frames()
    stacks = []
    real_vstack = models.sparse.vstack

    def spy_vstack(blocks, *args, **kwargs):
        stacks.append(len(list(blocks)))
        return real_vstack(list(blocks), *args, **kwargs)

    def forbidden(*args, **kwargs):
        raise AssertionError("to_anndata still uses the concatenating legacy path")

    monkeypatch.setattr(models.sparse, "vstack", spy_vstack)
    monkeypatch.setattr(models, "_arrays_to_anndata", forbidden)
    adata = CellExpressionStream(frames).to_anndata()

    assert adata.shape == (400, 20)
    # Every batch is kept as its own block and folded exactly once.
    assert stacks == [len(frames)]
    source = inspect.getsource(models.CellExpressionStream.to_anndata)
    assert "np.concatenate" not in source
    assert "_collect_parts" not in source
