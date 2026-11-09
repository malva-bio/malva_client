#!/usr/bin/env python3
"""Find the Malva metacells most similar to each cluster in an AnnData file.

The AnnData object stays local: the client averages cells per cluster and sends
only the resulting cluster pseudobulk profiles to Malva.

Example:
    export MALVA_API_TOKEN="..."
    python examples/compare_clusters_to_metacells.py \
        --h5ad my_data.h5ad --groupby leiden --dataset Human_Cell_Atlas
"""

from __future__ import annotations

import argparse
import os

import anndata as ad

from malva_client.client import MalvaClient


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True, help="AnnData .h5ad file to analyse.")
    parser.add_argument(
        "--groupby",
        default="leiden",
        help="Categorical obs column with the cluster labels (default: leiden).",
    )
    parser.add_argument(
        "--dataset",
        required=True,
        help="Malva dataset / coexpression index id, e.g. Human_Cell_Atlas.",
    )
    parser.add_argument("--top-n", type=int, default=50, help="Neighbours per cluster.")
    parser.add_argument(
        "--layer",
        default=None,
        help="AnnData layer to use instead of adata.X.",
    )
    parser.add_argument(
        "--url",
        default=os.environ.get("MALVA_URL", "http://127.0.0.1:18000"),
        help="Malva API base URL.",
    )
    args = parser.parse_args()

    adata = ad.read_h5ad(args.h5ad)
    client = MalvaClient(
        base_url=args.url,
        api_token=os.environ.get("MALVA_API_TOKEN"),
        timeout=300,
    )

    matches = client.find_metacells(
        adata,
        groupby=args.groupby,
        dataset_id=args.dataset,
        top_n=args.top_n,
        layer=args.layer,
    )
    print(matches.to_dataframe().head(30).to_string(index=False))

    labels = matches.assign_labels(
        adata,
        column="malva_cell_type",
        confidence_column="malva_confidence",
    )
    print(labels.to_string(index=False))
    print(adata.obs[[args.groupby, "malva_cell_type", "malva_confidence"]].head())


if __name__ == "__main__":
    main()
