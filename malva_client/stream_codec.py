"""Decoders for the opt-in cell-expression stream encodings.

The server's default framing is four flat little-endian arrays per data frame
(``sample_ids``/``cell_ids``/``feature_index``/``values``).  Two opt-in modes
shrink it, and both are advertised in the stream's ``start`` frame rather than
negotiated up front:

``layout="compact"``
    Cell identity comes from the request order.  Each data frame carries a
    bitmask over the requested cells, then per present cell its non-zero count
    and its feature positions as unsigned LEB128 varints, delta coded inside the
    row.  Values stay ``float32`` unless ``values="f16"`` asked for half
    precision.

``compress="zstd"|"zlib"``
    Every frame body after the ``start`` frame is compressed as one blob.

Nothing here runs for a default request: the client only decodes these forms
when the server says it produced them, so released servers and clients keep
working unchanged.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

#: Layout name advertised by the server in the ``start`` frame.
COMPACT_LAYOUT = "cells_compact_v1"


def unpack_leb128(buffer: bytes, count: int) -> np.ndarray:
    """Decode ``count`` unsigned LEB128 varints into an ``int64`` array."""
    raw = np.frombuffer(buffer, dtype=np.uint8)
    if count == 0:
        return np.empty(0, dtype=np.int64)
    ends = np.flatnonzero((raw & 0x80) == 0)
    if int(ends.size) != count:
        raise ValueError(f"expected {count} varints in the stream, found {int(ends.size)}")
    starts = np.empty(count, dtype=np.int64)
    starts[0] = 0
    starts[1:] = ends[:-1] + 1
    lengths = ends - starts + 1
    offsets = np.arange(int(raw.size), dtype=np.int64) - np.repeat(starts, lengths)
    shifted = (raw & 0x7F).astype(np.uint64) << (7 * offsets).astype(np.uint64)
    return np.add.reduceat(shifted, starts).astype(np.int64, copy=False)


def delta_decode_rows(deltas: np.ndarray, counts: np.ndarray) -> np.ndarray:
    """Undo the per-row delta coding of a concatenated position array."""
    if deltas.size == 0:
        return deltas.astype(np.int64, copy=False)
    cumulative = np.cumsum(deltas, dtype=np.int64)
    starts = np.empty(int(counts.size), dtype=np.int64)
    starts[0] = 0
    np.cumsum(counts[:-1], out=starts[1:])
    bases = np.empty(int(counts.size), dtype=np.int64)
    bases[0] = 0
    bases[1:] = cumulative[starts[1:] - 1]
    return cumulative - np.repeat(bases, counts)


def decompress_body(body: bytes, encoding: Optional[str]) -> bytes:
    """Decode one compressed frame body; ``None`` means it was not compressed."""
    if not encoding:
        return body
    name = str(encoding).strip().lower()
    if name == "zlib":
        import zlib

        return zlib.decompress(body)
    if name in ("zstd", "zstandard"):
        try:
            import zstandard

            return zstandard.ZstdDecompressor().decompressobj().decompress(body)
        except ImportError:
            pass
        try:
            import numcodecs

            return bytes(numcodecs.Zstd().decode(body))
        except ImportError as exc:  # pragma: no cover - depends on the install
            raise RuntimeError(
                "the server compressed the stream with zstd but neither "
                "'zstandard' nor 'numcodecs' is installed; install one of them "
                "(pip install zstandard) or request compress=None"
            ) from exc
    raise ValueError(f"unknown stream encoding {encoding!r}")


def decode_compact_frame(
    frame: Dict[str, Any],
    requested: Sequence[Tuple[int, int]],
    features: Sequence[str],
    values_dtype: str = "f32",
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Expand one compact data frame into the default framing's four arrays.

    Returns ``(sample_ids, cell_ids, feature_names, values)``, each one entry
    per non-zero, in the order the default framing would have emitted them
    (request order, ascending features within a cell).
    """
    n_requested = len(requested)
    mask = np.unpackbits(np.frombuffer(frame["cells"], dtype=np.uint8), count=n_requested)
    present = np.flatnonzero(mask[:n_requested])
    nnz = unpack_leb128(frame["nnz"], int(present.size))
    total = int(nnz.sum())
    if total == 0:
        empty = np.empty(0, dtype=np.int64)
        return empty, empty, np.empty(0, dtype=object), np.empty(0, dtype=np.float32)
    width = int(frame.get("feature_width", 4))
    if width not in (2, 4):
        raise ValueError(f"unsupported compact feature width {width}")
    positions = np.frombuffer(frame["features"], dtype=f"<u{width}").astype(np.int64)
    if positions.size != total:
        raise ValueError(
            f"compact frame declares {total} positions but carries {positions.size}"
        )
    if frame.get("feature_mode") == "delta":
        positions = delta_decode_rows(positions, nnz)
    dtype = "<f2" if values_dtype == "f16" else "<f4"
    values = np.frombuffer(frame["values"], dtype=dtype)
    if values.size != total:
        raise ValueError(f"compact frame declares {total} values but carries {values.size}")
    if values_dtype == "f16":
        values = values.astype(np.float32)
    samples = np.fromiter((requested[i][0] for i in present.tolist()), np.int64, present.size)
    cells = np.fromiter((requested[i][1] for i in present.tolist()), np.int64, present.size)
    if features is not None and len(features):
        index = positions
        highest = int(index.max(initial=-1))
        if highest >= len(features):
            raise ValueError(
                f"compact frame references feature {highest} "
                f"but the store has {len(features)}"
            )
        names = np.asarray(features, dtype=object)[index]
    else:
        names = np.asarray([f"feature_{p}" for p in positions.tolist()], dtype=object)
    return np.repeat(samples, nnz), np.repeat(cells, nnz), names, values


def map_gene_names(
    positions: np.ndarray, features: Sequence[str], requested_genes: List[str]
) -> np.ndarray:
    """Vectorised equivalent of the per-entry gene-name mapping.

    Released servers send positions into the frame's ``features`` list; very old
    ones send store ordinals, resolved against the requested panel.  Both are
    served by one lookup array instead of a Python loop over every entry.
    """
    index = np.asarray(positions, dtype=np.int64)
    count = int(index.size)
    if count == 0:
        return np.empty(0, dtype=object)
    if features is not None and len(features) and int(index.max(initial=-1)) < len(features):
        return np.asarray(features, dtype=object)[index]
    unique = np.unique(index)
    if requested_genes and int(unique.size) == len(requested_genes):
        table = np.asarray(requested_genes, dtype=object)
        slots = np.searchsorted(unique, index)
        return table[slots]
    return np.asarray([f"feature_{p}" for p in index.tolist()], dtype=object)
