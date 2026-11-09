#!/usr/bin/env python3
"""Stream sparse per-cell expression for selected samples.

Example:
    export MALVA_API_TOKEN="$(security find-generic-password -s malva_client_test -a MALVA_API_TOKEN -w)"
    python examples/cell_expression_stream.py \
        --sample-id 30064771072 \
        --gene CD3D --gene CD3E --gene CD247
"""

from __future__ import annotations

import argparse
import os

from malva_client.client import MalvaClient


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sample-id",
        action="append",
        type=int,
        required=True,
        help="Global sample id to stream (repeat for multiple samples).",
    )
    parser.add_argument(
        "--gene",
        action="append",
        default=[],
        help="Gene symbol to include (repeat for a panel). Omit to stream all genes.",
    )
    parser.add_argument(
        "--dataset",
        default=None,
        help="Dataset id, e.g. Human_Cell_Atlas. Required only when multiple datasets exist.",
    )
    parser.add_argument(
        "--url",
        default=os.environ.get("MALVA_URL", "http://127.0.0.1:18000"),
        help="Malva API base URL.",
    )
    parser.add_argument("--batch-bytes", type=int, default=4 * 1024 * 1024)
    args = parser.parse_args()

    client = MalvaClient(
        base_url=args.url,
        api_token=os.environ.get("MALVA_API_TOKEN"),
        timeout=300,
    )

    total = 0
    for batch in client.stream_cell_expression(
        sample_ids=args.sample_id,
        gene_ids=args.gene or None,
        dataset_id=args.dataset,
        max_batch_bytes=args.batch_bytes,
    ):
        total += len(batch)
        print(batch.head().to_string(index=False))

    print(f"Streamed {total:,} non-zero values.")


if __name__ == "__main__":
    main()
