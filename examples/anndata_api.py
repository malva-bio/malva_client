#!/usr/bin/env python3
"""Convert Malva results to standard AnnData objects.

This example shows the simplified, AnnData-friendly API:

* ``SearchResult.enrich().to_anndata()``
* ``CellExpressionStream.to_anndata()``
* ``CellExpressionStream.write_h5ad(path)``
* ``CellExpressionStream.iter_anndata()`` + ``anndata.concat``

Example:
    export MALVA_API_TOKEN="..."
    python examples/anndata_api.py \
        --gene SPP1 --gene CD3D --sample-id 30064771072 \
        --dataset Human_Cell_Atlas --output cells.h5ad
"""

from __future__ import annotations

import argparse
import os

import anndata as ad

from malva_client.client import MalvaClient


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene", action="append", default=[], help="Gene symbol (repeatable).")
    parser.add_argument(
        "--sample-id",
        action="append",
        type=int,
        default=[],
        help="Sample id to stream (repeatable).",
    )
    parser.add_argument("--dataset", default=None, help="Malva dataset id.")
    parser.add_argument("--output", default="cells.h5ad", help="Output .h5ad path.")
    parser.add_argument(
        "--url",
        default=os.environ.get("MALVA_URL", "http://127.0.0.1:18000"),
        help="Malva API base URL.",
    )
    args = parser.parse_args()

    client = MalvaClient(
        base_url=args.url,
        api_token=os.environ.get("MALVA_API_TOKEN"),
        timeout=300,
    )

    # 1. Aggregate search result -> sample × cell-type AnnData.
    if args.gene:
        query = ", ".join(args.gene)
        search_adata = client.search(query).enrich().to_anndata()
        print("Search AnnData:")
        print(search_adata)

    if not args.sample_id:
        return

    # 2. Cell-expression stream -> cells × genes AnnData.
    stream = client.stream_cell_expression(
        sample_ids=args.sample_id,
        gene_ids=args.gene or None,
        dataset_id=args.dataset,
    )
    cells_adata = stream.to_anndata()
    print("Streaming AnnData:")
    print(cells_adata)

    # 3. Bounded-memory write and lazy read.
    client.stream_cell_expression(
        sample_ids=args.sample_id,
        gene_ids=args.gene or None,
        dataset_id=args.dataset,
    ).write_h5ad(args.output)
    backed = ad.read_h5ad(args.output, backed="r")
    print("Backed file:")
    print(backed)
    backed.file.close()

    # 4. Build the matrix yourself from per-batch parts.
    parts = client.stream_cell_expression(
        sample_ids=args.sample_id,
        gene_ids=args.gene or None,
        dataset_id=args.dataset,
    ).iter_anndata()
    combined = ad.concat(parts, join="outer", fill_value=0)
    print("Concatenated batches:")
    print(combined)


if __name__ == "__main__":
    main()
