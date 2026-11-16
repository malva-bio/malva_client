"""Decoders for the cell-expression stream encodings.

The server's released framing is four flat little-endian arrays per data frame
(``sample_ids``/``cell_ids``/``feature_index``/``values``).  Two encodings shrink
it, and both are advertised in the stream's ``start`` frame rather than
negotiated up front:

``layout="compact"``
    Cell identity comes from the request order.  Each data frame carries a
    bitmask over the requested cells, then per present cell its non-zero count
    and its feature positions as fixed-width unsigned integers, delta coded
    inside the row.  Values stay ``float32`` unless ``values="f16"`` asked for
    half precision.

``compress="zstd"|"zlib"``
    Every frame body after the ``start`` frame is compressed as one blob.

``stream_cell_expression(cells=...)`` asks for compact plus zstd by default;
``layout="legacy"`` and ``compress=None`` restore the released framing.  The
decoder is driven by the ``start`` frame, never by the request, so a server that
ignores the options still decodes correctly.  Decoding stays vectorised: every
compact field is one ``np.frombuffer`` over the raw column plus
``np.repeat``/``np.cumsum``, with no Python loop over values.
"""

from __future__ import annotations

import threading
import zlib
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

#: Layout name advertised by the server in the ``start`` frame.
COMPACT_LAYOUT = "cells_compact_v1"

#: Decoders cache one codec object per thread.  ``zstandard`` codec objects
#: keep scratch buffers and are not safe to share across threads; the zlib and
#: numcodecs paths benefit from the same cache because a stream decodes one
#: frame after another.
_DECODER_LOCAL = threading.local()


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


def _zstd_decoder() -> Tuple[Optional[str], Any]:
    """Return ``(name, decoder)`` for the first importable zstd implementation."""
    cached = getattr(_DECODER_LOCAL, "zstd", None)
    if cached is not None:
        return cached
    decoder: Tuple[Optional[str], Any] = (None, None)
    try:
        import zstandard

        decoder = ("zstandard", zstandard.ZstdDecompressor())
    except ImportError:
        try:
            import numcodecs

            decoder = ("numcodecs", numcodecs.Zstd())
        except ImportError:
            decoder = (None, None)
    _DECODER_LOCAL.zstd = decoder
    return decoder


def has_zstd() -> bool:
    """Whether this install can decode zstd frame bodies.

    ``zstandard`` is preferred; ``numcodecs`` (already present with the zarr
    stack most users have) ships the same codec.  When neither is importable a
    defaulted request keeps the compact framing and sends no compression, which
    is still smaller and faster than the released framing.
    """
    return _zstd_decoder()[0] is not None


def decompress_body(body: bytes, encoding: Optional[str]) -> bytes:
    """Decode one compressed frame body; ``None`` means it was not compressed."""
    if not encoding:
        return body
    name = str(encoding).strip().lower()
    if name == "zlib":
        return zlib.decompress(body)
    if name in ("zstd", "zstandard"):
        codec_name, decoder = _zstd_decoder()
        if decoder is None:
            raise RuntimeError(
                "the server compressed the stream with zstd but neither "
                "'zstandard' nor 'numcodecs' is installed; install one of them "
                "(pip install zstandard) or request compress=None"
            )
        if codec_name == "zstandard":
            return decoder.decompress(body)
        return bytes(decoder.decode(body))
    raise ValueError(f"unknown stream encoding {encoding!r}")


def _feature_names(positions: np.ndarray, features: Sequence[str]) -> np.ndarray:
    """Map panel positions to names, refusing positions the panel cannot hold."""
    highest = int(positions.max(initial=-1))
    lowest = int(positions.min(initial=0))
    if features is not None and len(features):
        if lowest < 0 or highest >= len(features):
            raise ValueError(
                f"compact frame references feature {highest if highest >= 0 else lowest} "
                f"but the store has {len(features)}"
            )
        return np.asarray(features, dtype=object)[positions]
    return np.asarray([f"feature_{p}" for p in positions.tolist()], dtype=object)


def decode_compact_columns(
    frame: Dict[str, Any],
    requested: Any,
    values_dtype: str = "f32",
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Expand a compact frame into aligned columns, without naming features.

    Returns ``(present, nnz, sample_ids, cell_ids, positions, values)``.
    ``present`` indexes the request order and ``nnz`` is its per-cell non-zero
    count. ``requested`` may be a sequence of pairs or a pre-built ``(n, 2)``
    integer array; the latter avoids re-converting it for every frame.
    """
    requested_array = (
        requested
        if isinstance(requested, np.ndarray) and requested.ndim == 2
        else np.asarray(requested, dtype=np.int64)
    )
    n_requested = int(requested_array.shape[0])
    raw_mask = np.frombuffer(frame["cells"], dtype=np.uint8)
    if int(raw_mask.size) * 8 < n_requested:
        raise ValueError(
            f"compact frame carries a {int(raw_mask.size) * 8}-bit mask "
            f"for {n_requested} requested cells"
        )
    mask = np.unpackbits(raw_mask, count=n_requested)
    present = np.flatnonzero(mask[:n_requested])
    nnz = unpack_leb128(frame["nnz"], int(present.size))
    if int(nnz.min(initial=0)) < 0:
        raise ValueError("compact frame carries a negative non-zero count")
    total = int(nnz.sum())
    declared = int(frame.get("n") or 0)
    if declared and declared != total:
        raise ValueError(f"compact frame declares n={declared} but its counts sum to {total}")
    empty = np.empty(0, dtype=np.int64)
    if total == 0:
        return (
            present,
            nnz,
            empty,
            empty,
            empty,
            np.empty(0, dtype=np.float32),
        )
    width = int(frame.get("feature_width", 4))
    if width not in (2, 4):
        raise ValueError(f"unsupported compact feature width {width}")
    positions = np.frombuffer(frame["features"], dtype=f"<u{width}").astype(np.int64)
    if positions.size != total:
        raise ValueError(f"compact frame declares {total} positions but carries {positions.size}")
    if frame.get("feature_mode") == "delta":
        positions = delta_decode_rows(positions, nnz)
    dtype = "<f2" if values_dtype == "f16" else "<f4"
    values = np.frombuffer(frame["values"], dtype=dtype)
    if values.size != total:
        raise ValueError(f"compact frame declares {total} values but carries {values.size}")
    if values_dtype == "f16":
        values = values.astype(np.float32)
    samples = np.repeat(requested_array[present, 0], nnz)
    cells = np.repeat(requested_array[present, 1], nnz)
    return present, nnz, samples, cells, positions, values


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
    _present, _nnz, samples, cells, positions, values = decode_compact_columns(
        frame, requested, values_dtype
    )
    return samples, cells, _feature_names(positions, features), values


class ExpressionBatch:
    """One decoded server frame with feature codes instead of repeated names.

    Every non-zero entry carries its cell in ``sample_ids``/``cell_ids`` and its
    gene as an integer code into ``gene_names``. That keeps the per-entry cost
    to three numeric arrays; names are resolved once per batch when a caller
    asks for the public DataFrame or when the AnnData accumulator needs a new
    gene. ``new_codes`` lists the codes that have not appeared in an earlier
    batch of the same stream, in first-appearance order.
    """

    __slots__ = (
        "sample_ids",
        "cell_ids",
        "gene_codes",
        "values",
        "store_index",
        "gene_names",
        "new_codes",
    )

    def __init__(
        self,
        sample_ids: np.ndarray,
        cell_ids: np.ndarray,
        gene_codes: np.ndarray,
        values: np.ndarray,
        store_index: int,
        gene_names: Sequence[str],
        new_codes: np.ndarray,
    ):
        self.sample_ids = sample_ids
        self.cell_ids = cell_ids
        self.gene_codes = gene_codes
        self.values = values
        self.store_index = int(store_index)
        self.gene_names = gene_names
        self.new_codes = new_codes

    def __len__(self) -> int:
        return int(self.values.size)

    def to_frame(self) -> Any:
        """Return the public long-format DataFrame for this frame."""
        import pandas as pd

        names = np.asarray(self.gene_names, dtype=object)[self.gene_codes]
        count = int(self.values.size)
        return pd.DataFrame(
            {
                "sample_id": self.sample_ids,
                "cell_id": self.cell_ids,
                "store_index": np.full(count, self.store_index, dtype=np.int32),
                "gene": names,
                "value": self.values,
            }
        )


def map_gene_names(
    positions: np.ndarray, features: Sequence[str], requested_genes: List[str]
) -> np.ndarray:
    """Vectorised equivalent of the per-entry gene-name mapping.

    Released servers send positions into the frame's ``features`` list; very old
    ones send store ordinals, resolved against the requested panel.  The common
    path is one fancy index into the feature list instead of a Python loop over
    every entry; the rare ordinal path keeps the released client's exact
    first-appearance semantics.
    """
    index = np.asarray(positions, dtype=np.int64)
    count = int(index.size)
    if count == 0:
        return np.empty(0, dtype=object)
    if features is not None and len(features):
        highest = int(index.max(initial=-1))
        lowest = int(index.min(initial=0))
        if lowest >= 0 and highest < len(features):
            return np.asarray(features, dtype=object)[index]
    unique = list(dict.fromkeys(index.tolist()))
    names: Dict[int, str] = {}
    if requested_genes and len(unique) == len(requested_genes):
        names = {position: requested_genes[k] for k, position in enumerate(unique)}
    return np.asarray([names.get(position, f"feature_{position}") for position in index.tolist()])
