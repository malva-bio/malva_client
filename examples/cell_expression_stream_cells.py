#!/usr/bin/env python3
"""Stream full gene expression for an ordered list of cells.

Cell ids are ``(sample_id, cell_id)`` pairs from a search result.  The input
order is preserved, so you can sort cells by expression and fetch them in
high-to-low batches.

The dataset must have a CSR (cell-major) expression store.  If it does not,
the API returns a clear "not available yet" error.
"""

from __future__ import annotations

import argparse
import os

from malva_client.client import MalvaClient


def _parse_cell(value: str) -> tuple[int, int]:
    try:
        sample, cell = value.split(":", 1)
        return int(sample), int(cell)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError("cells must look like SAMPLE_ID:CELL_ID") from exc


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cell",
        action="append",
        required=True,
        type=_parse_cell,
        help="Cell as SAMPLE_ID:CELL_ID; repeat in the order you want",
    )
    parser.add_argument("--gene", action="append", help="Gene symbol; repeat to select a panel")
    parser.add_argument("--dataset", default="Human_Cell_Atlas")
    parser.add_argument("--url", default=os.environ.get("MALVA_URL", "http://127.0.0.1:18000"))
    parser.add_argument("--batch-bytes", type=int, default=4 * 1024 * 1024)
    args = parser.parse_args()

    client = MalvaClient(
        base_url=args.url,
        api_token=os.environ.get("MALVA_API_TOKEN"),
        timeout=300,
    )

    total = 0
    for batch in client.stream_cell_expression(
        cells=args.cell,
        gene_ids=args.gene or None,
        dataset_id=args.dataset,
        max_batch_bytes=args.batch_bytes,
    ):
        total += len(batch)
        print(batch.head())

    print(f"Streamed {total:,} non-zero values.")


if __name__ == "__main__":
    main()
