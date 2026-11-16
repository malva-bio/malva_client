"""Opt-in stream encodings: varints, compact frames, compression, name mapping."""

from __future__ import annotations

import importlib
import struct
import zlib

import msgpack
import numpy as np
import pandas as pd
import pytest

from malva_client import stream_codec
from malva_client.client import _map_stream_gene_names
from malva_client.exceptions import MalvaAPIError


def _leb128_encode(values) -> bytes:
    """Reference encoder, used only by the tests."""
    out = bytearray()
    for value in values:
        while True:
            byte = value & 0x7F
            value >>= 7
            out.append(byte | (0x80 if value else 0))
            if not value:
                break
    return bytes(out)


def _delta_encode(positions: np.ndarray, counts: np.ndarray) -> np.ndarray:
    deltas = positions.astype(np.int64, copy=True)
    deltas[1:] -= positions[:-1]
    starts = np.concatenate(([0], np.cumsum(counts)[:-1])).astype(np.int64)
    deltas[starts] = positions[starts]
    return deltas


class TestVarints:
    def test_round_trip(self):
        values = [0, 1, 127, 128, 300, 16383, 16384, 2**21, 2**31 - 1, 0]
        encoded = _leb128_encode(values)

        decoded = stream_codec.unpack_leb128(encoded, len(values))

        assert decoded.tolist() == values

    def test_multi_byte_boundaries(self):
        payload = _leb128_encode([208, 8090, 1113984])
        decoded = stream_codec.unpack_leb128(payload, 3)
        assert decoded.tolist() == [208, 8090, 1113984]
        assert decoded.dtype == np.int64

    def test_empty_buffer(self):
        assert stream_codec.unpack_leb128(b"", 0).size == 0

    def test_count_mismatch_is_rejected(self):
        with pytest.raises(ValueError):
            stream_codec.unpack_leb128(_leb128_encode([1, 2, 3]), 2)


class TestRowDeltas:
    def test_round_trip_with_uneven_rows(self):
        counts = np.array([3, 1, 4], dtype=np.int64)
        positions = np.array([2, 5, 9, 7, 0, 1, 2, 40], dtype=np.int32)

        deltas = _delta_encode(positions, counts)
        restored = stream_codec.delta_decode_rows(deltas, counts)

        assert restored.tolist() == positions.tolist()

    def test_single_row_keeps_absolute_first_value(self):
        counts = np.array([4], dtype=np.int64)
        positions = np.array([100, 101, 130, 131], dtype=np.int32)
        assert _delta_encode(positions, counts).tolist() == [100, 1, 29, 1]


class TestDecompression:
    def test_plain_body_is_returned_unchanged(self):
        assert stream_codec.decompress_body(b"abc", None) == b"abc"

    def test_zlib_round_trip(self):
        body = b"x" * 1000
        assert stream_codec.decompress_body(zlib.compress(body), "zlib") == body

    def test_unknown_encoding_is_rejected(self):
        with pytest.raises(ValueError):
            stream_codec.decompress_body(b"abc", "brotli")

    def test_zstd_round_trip_when_a_codec_is_installed(self):
        if not stream_codec.has_zstd():  # pragma: no cover - depends on the install
            pytest.skip("no zstd codec installed")
        body = b"y" * 500
        assert stream_codec.decompress_body(_compress(body, "zstd"), "zstd") == body


class TestCompactFrame:
    def _frame(self, requested, rows, values_dtype="<f4", mode="delta", width=2):
        counts = np.array([len(row) for row in rows], dtype=np.int64)
        positions = np.array([p for row in rows for p, _v in row], dtype=np.int64)
        values = np.array([v for row in rows for _p, v in row], dtype=np.float32)
        # Present cells are a subsequence of the request order (a cell with no
        # non-zeros emits nothing), so take indices that always fit.
        stride = max(1, len(requested) // max(1, len(rows)))
        present = np.arange(len(rows)) * stride
        mask = np.zeros(len(requested), dtype=bool)
        mask[present] = True
        encoded = _delta_encode(positions, counts) if mode == "delta" else positions
        return {
            "type": "data",
            "store_index": 0,
            "n": int(values.size),
            "cells": np.packbits(mask).tobytes(),
            "nnz": _leb128_encode(counts.tolist()),
            "feature_width": width,
            "feature_mode": mode,
            "features": encoded.astype(f"<u{width}").tobytes(),
            "values": values.astype(values_dtype).tobytes(),
        }

    def test_expands_to_the_default_framing(self):
        requested = [(7, 10), (7, 11), (7, 12), (7, 13)]
        rows = [[(0, 1.5), (3, 2.5)], [(4, 0.5)]]
        frame = self._frame(requested, rows)

        samples, cells, genes, values = stream_codec.decode_compact_frame(
            frame, requested, ["G0", "G1", "G2", "G3", "G4"]
        )

        assert samples.tolist() == [7, 7, 7]
        assert cells.tolist() == [10, 10, 12]
        assert list(genes) == ["G0", "G3", "G4"]
        assert values.tolist() == [1.5, 2.5, 0.5]

    def test_absolute_feature_mode(self):
        requested = [(1, 2), (1, 3)]
        frame = self._frame(requested, [[(5, 1.0)], [(2, 2.0)]], mode="absolute")
        _s, _c, genes, values = stream_codec.decode_compact_frame(
            frame, requested, ["G0", "G1", "G2", "G3", "G4", "G5"]
        )
        assert list(genes) == ["G5", "G2"]
        assert values.tolist() == [1.0, 2.0]

    def test_float16_values_are_upcast(self):
        requested = [(1, 2), (1, 3)]
        frame = self._frame(requested, [[(1, 0.25), (2, 0.5)]], values_dtype="<f2")
        _s, _c, _g, values = stream_codec.decode_compact_frame(
            frame, requested, ["G0", "G1", "G2"], "f16"
        )
        assert values.dtype == np.float32
        assert values.tolist() == [0.25, 0.5]

    def test_empty_frame(self):
        requested = [(1, 2), (1, 3)]
        empty = np.empty(0, dtype=np.int64)
        frame = {
            "type": "data",
            "store_index": 0,
            "n": 0,
            "cells": np.packbits(np.zeros(len(requested), dtype=bool)).tobytes(),
            "nnz": _leb128_encode([]),
            "feature_width": 2,
            "feature_mode": "delta",
            "features": b"",
            "values": b"",
        }
        samples, cells, _g, values = stream_codec.decode_compact_frame(frame, requested, ["G0"])
        assert samples.size == 0 and cells.size == 0 and values.size == 0
        assert empty.size == 0


class TestGeneNameMapping:
    @pytest.mark.parametrize(
        "positions,features,requested",
        [
            ([0, 2, 1], ["G0", "G1", "G2"], []),
            ([5, 9], [], ["G5", "G9"]),
            ([], ["G0"], []),
        ],
    )
    def test_matches_the_reference_mapping(self, positions, features, requested):
        expected = _map_stream_gene_names(positions, features, requested)
        got = stream_codec.map_gene_names(
            np.asarray(positions, dtype=np.int32), features, requested
        )
        assert list(got) == expected

    def test_out_of_range_position_falls_back(self):
        got = stream_codec.map_gene_names(np.asarray([7], dtype=np.int32), ["G0"], [])
        assert list(got) == ["feature_7"]


class _FakeResponse:
    def __init__(self, payload: bytes, status_code: int = 200):
        self._payload = payload
        self.status_code = status_code
        self.ok = status_code < 400
        self.text = ""

    def iter_content(self, chunk_size: int):
        for start in range(0, len(self._payload), chunk_size):
            yield self._payload[start : start + chunk_size]

    def close(self):
        return None


def _compress(body: bytes, encoding: str) -> bytes:
    if encoding == "zlib":
        return zlib.compress(body)
    if encoding == "zstd":
        try:
            zstandard = importlib.import_module("zstandard")
        except ImportError:  # pragma: no cover - depends on the install
            import numcodecs

            return bytes(numcodecs.Zstd().encode(body))
        return zstandard.ZstdCompressor(level=1).compress(body)
    raise ValueError(encoding)


def _stream_payload(frames, encoding=None):
    out = bytearray()
    for index, frame in enumerate(frames):
        body = msgpack.packb(frame, use_bin_type=True)
        if index and encoding:
            body = _compress(body, encoding)
        out += struct.pack(">Q", len(body)) + body
    return bytes(out)


class TestClientStreamOptions:
    def _client(self, monkeypatch, mock_client, frames, encoding=None, captured=None):
        client = mock_client

        def fake_post(url, json=None, stream=None, timeout=None):
            if captured is not None:
                captured.update(json or {})
            return _FakeResponse(_stream_payload(frames, encoding))

        monkeypatch.setattr(client.session, "post", fake_post)
        return client

    def _frames(self, layout=None, compression=None, values="f32"):
        start = {"type": "start", "n_cells": 2, "gene_filter": None}
        if layout:
            start["layout"] = layout
            start["values"] = values
        if compression:
            start["compression"] = compression
        if layout:
            data = {
                "type": "data",
                "store_index": 0,
                "n": 3,
                "cells": np.packbits(np.array([True, True])).tobytes(),
                "nnz": _leb128_encode([2, 1]),
                "feature_width": 2,
                "feature_mode": "delta",
                "features": np.asarray([0, 1, 1], dtype="<u2").tobytes(),
                "values": np.asarray([1.0, 2.0, 3.0], dtype="<f4").tobytes(),
            }
        else:
            data = {
                "type": "data",
                "store_index": 0,
                "n": 3,
                "sample_ids": np.asarray([7, 7, 7], dtype="<i8").tobytes(),
                "cell_ids": np.asarray([10, 10, 11], dtype="<i8").tobytes(),
                "feature_index": np.asarray([0, 1, 1], dtype="<i4").tobytes(),
                "values": np.asarray([1.0, 2.0, 3.0], dtype="<f4").tobytes(),
            }
        return [
            start,
            {"type": "features", "store_index": 0, "features": ["G0", "G1"]},
            data,
            {"type": "done", "n_frames": 1, "n_values": 3},
        ]

    def test_compact_frames_decode_to_the_same_dataframe(self, monkeypatch, mock_client):
        captured = {}
        client = self._client(
            monkeypatch, mock_client, self._frames("cells_compact_v1"), captured=captured
        )

        batches = list(
            client.stream_cell_expression(
                cells=[(7, 10), (7, 11)], dataset_id="ds1", layout="compact"
            )
        )

        assert captured["layout"] == "compact"
        frame = pd.concat(batches, ignore_index=True)
        assert frame.columns.tolist() == [
            "sample_id",
            "cell_id",
            "store_index",
            "gene",
            "value",
        ]
        assert frame["cell_id"].tolist() == [10, 10, 11]
        assert frame["gene"].tolist() == ["G0", "G1", "G1"]
        assert frame["value"].tolist() == [1.0, 2.0, 3.0]

    def test_compressed_compact_stream(self, monkeypatch, mock_client):
        captured = {}
        client = self._client(
            monkeypatch,
            mock_client,
            self._frames("cells_compact_v1", "zlib"),
            encoding="zlib",
            captured=captured,
        )

        frame = pd.concat(
            client.stream_cell_expression(
                cells=[(7, 10), (7, 11)],
                dataset_id="ds1",
                layout="compact",
                compress="zlib",
            ),
            ignore_index=True,
        )

        assert captured["compress"] == "zlib"
        assert frame["value"].tolist() == [1.0, 2.0, 3.0]

    def test_default_request_sends_no_options(self, monkeypatch, mock_client):
        captured = {}
        client = self._client(monkeypatch, mock_client, self._frames(), captured=captured)

        list(client.stream_cell_expression(cells=[(7, 10), (7, 11)], dataset_id="ds1"))

        assert "layout" not in captured
        assert "compress" not in captured
        assert "values" not in captured

    def test_options_need_a_cell_list(self, mock_client):
        with pytest.raises(ValueError):
            mock_client.stream_cell_expression(sample_ids=[1], compress="zstd")

    def test_half_precision_needs_the_compact_layout(self, mock_client):
        with pytest.raises(ValueError):
            mock_client.stream_cell_expression(cells=[(7, 10)], values="f16")

    @pytest.mark.parametrize("encoding", ["zstd", "zlib"])
    def test_compressed_modes_decode_to_the_same_frame(self, monkeypatch, mock_client, encoding):
        captured = {}
        client = self._client(
            monkeypatch,
            mock_client,
            self._frames("cells_compact_v1", encoding),
            encoding=encoding,
            captured=captured,
        )

        frame = pd.concat(
            client.stream_cell_expression(
                cells=[(7, 10), (7, 11)],
                dataset_id="ds1",
                layout="compact",
                compress=encoding,
            ),
            ignore_index=True,
        )

        assert captured["compress"] == encoding
        assert frame["value"].tolist() == [1.0, 2.0, 3.0]
        assert frame["gene"].tolist() == ["G0", "G1", "G1"]

    def test_zlib_does_not_need_a_zstd_decoder(self, monkeypatch, mock_client):
        captured = {}
        client = self._client(
            monkeypatch,
            mock_client,
            self._frames("cells_compact_v1", "zlib"),
            encoding="zlib",
            captured=captured,
        )
        monkeypatch.setattr(stream_codec, "has_zstd", lambda: False)

        frame = pd.concat(
            client.stream_cell_expression(
                cells=[(7, 10), (7, 11)],
                dataset_id="ds1",
                layout="compact",
                compress="zlib",
            ),
            ignore_index=True,
        )

        assert captured["compress"] == "zlib"
        assert frame["value"].tolist() == [1.0, 2.0, 3.0]

    def test_missing_zstd_decoder_falls_back_to_legacy(self, monkeypatch, mock_client):
        captured = {}
        client = self._client(monkeypatch, mock_client, self._frames(), captured=captured)
        monkeypatch.setattr(stream_codec, "has_zstd", lambda: False)

        frame = pd.concat(
            client.stream_cell_expression(
                cells=[(7, 10), (7, 11)],
                dataset_id="ds1",
                layout="compact",
                compress="zstd",
            ),
            ignore_index=True,
        )

        assert "compress" not in captured
        assert "layout" not in captured
        assert frame["value"].tolist() == [1.0, 2.0, 3.0]

    def test_old_server_that_ignores_compact_is_decoded_as_legacy(
        self, monkeypatch, mock_client, caplog
    ):
        captured = {}
        client = self._client(monkeypatch, mock_client, self._frames(), captured=captured)

        with caplog.at_level("DEBUG", logger="malva_client"):
            frame = pd.concat(
                client.stream_cell_expression(
                    cells=[(7, 10), (7, 11)], dataset_id="ds1", layout="compact"
                ),
                ignore_index=True,
            )

        assert captured["layout"] == "compact"
        assert frame["gene"].tolist() == ["G0", "G1", "G1"]
        assert "did not honour layout" in caplog.text

    def test_corrupt_compact_frame_raises_a_clear_error(self, monkeypatch, mock_client):
        frames = self._frames("cells_compact_v1")
        frames[2] = dict(frames[2], feature_width=3)
        client = self._client(monkeypatch, mock_client, frames)

        with pytest.raises(MalvaAPIError) as excinfo:
            list(client.stream_cell_expression(cells=[(7, 10), (7, 11)], layout="compact"))

        message = str(excinfo.value)
        assert "frame 2" in message
        assert "feature width" in message

    def test_declared_count_mismatch_raises(self, monkeypatch, mock_client):
        frames = self._frames()
        frames[2] = dict(frames[2], n=99)
        client = self._client(monkeypatch, mock_client, frames)

        with pytest.raises(MalvaAPIError) as excinfo:
            list(client.stream_cell_expression(cells=[(7, 10), (7, 11)]))

        assert "declares 99" in str(excinfo.value)

    def test_compact_and_legacy_streams_build_the_same_anndata(self, monkeypatch, mock_client):
        legacy_client = self._client(monkeypatch, mock_client, self._frames())
        legacy = legacy_client.stream_cell_expression(
            cells=[(7, 10), (7, 11)], dataset_id="ds1", include_cell_metadata=False
        ).to_anndata()

        compact_client = self._client(monkeypatch, mock_client, self._frames("cells_compact_v1"))
        compact = compact_client.stream_cell_expression(
            cells=[(7, 10), (7, 11)],
            dataset_id="ds1",
            layout="compact",
            include_cell_metadata=False,
        ).to_anndata()

        assert legacy.obs_names.tolist() == compact.obs_names.tolist()
        assert legacy.var_names.tolist() == compact.var_names.tolist()
        np.testing.assert_array_equal(legacy.X.toarray(), compact.X.toarray())
        assert compact.X.dtype == legacy.X.dtype


class TestDecodeValidation:
    def test_values_count_mismatch_is_rejected(self):
        requested = [(1, 2)]
        frame = {
            "type": "data",
            "n": 2,
            "cells": np.packbits(np.array([True])).tobytes(),
            "nnz": _leb128_encode([2]),
            "feature_width": 2,
            "feature_mode": "absolute",
            "features": np.asarray([0, 1], dtype="<u2").tobytes(),
            "values": np.asarray([1.0], dtype="<f4").tobytes(),
        }
        with pytest.raises(ValueError):
            stream_codec.decode_compact_frame(frame, requested, ["G0", "G1"])

    def test_position_beyond_the_panel_is_rejected(self):
        requested = [(1, 2)]
        frame = {
            "type": "data",
            "n": 1,
            "cells": np.packbits(np.array([True])).tobytes(),
            "nnz": _leb128_encode([1]),
            "feature_width": 2,
            "feature_mode": "absolute",
            "features": np.asarray([9], dtype="<u2").tobytes(),
            "values": np.asarray([1.0], dtype="<f4").tobytes(),
        }
        with pytest.raises(ValueError):
            stream_codec.decode_compact_frame(frame, requested, ["G0", "G1"])

    def test_negative_positions_are_not_wrapped_by_the_names(self):
        got = stream_codec.map_gene_names(np.asarray([-1], dtype=np.int64), ["G0", "G1"], [])
        assert list(got) == ["feature_-1"]

    def test_large_feature_indices_decode(self):
        requested = [(1, 2), (3, 4)]
        values = np.asarray([1.0, 2.0], dtype="<f4")
        features = np.asarray([65_536, 65_537], dtype="<u4")
        frame = {
            "type": "data",
            "n": 2,
            "cells": np.packbits(np.array([True, True])).tobytes(),
            "nnz": _leb128_encode([1, 1]),
            "feature_width": 4,
            "feature_mode": "delta",
            "features": features.tobytes(),
            "values": values.tobytes(),
        }
        _s, _c, genes, decoded = stream_codec.decode_compact_frame(
            frame, requested, [f"G{i}" for i in range(65_538)]
        )
        assert list(genes) == ["G65536", "G65537"]
        np.testing.assert_array_equal(decoded, values)
