# Malva Client

[![PyPI](https://img.shields.io/pypi/v/malva-client)](https://pypi.org/project/malva-client/)
[![Documentation](https://readthedocs.org/projects/malva-client/badge/?version=latest)](https://malva-client.readthedocs.io)
[![GitHub](https://img.shields.io/badge/GitHub-malva--bio%2Fmalva__client-blue)](https://github.com/malva-bio/malva_client)

Python client for the [Malva](https://malva.bio) genomic search platform. Search genes, sequences, and natural language queries across >7,000 single-cell and spatial transcriptomics samples.

For full documentation, visit [malva-client.readthedocs.io](https://malva-client.readthedocs.io).

## Installation

```bash
pip install malva-client
```

For single-cell analysis workflows:
```bash
pip install malva-client scanpy
```

## Authentication

Generate an API token at [malva.bio](https://malva.bio) (login with ORCID, then go to Profile > Generate API Token), then configure the client:

```bash
malva_client config --server https://malva.mdc-berlin.de --token YOUR_API_TOKEN
```

## Quick Start (CLI)

```bash
malva_client search "CD3D" --output results.csv
malva_client search "ATCGATCGATCGATCGATCGATCG" --output sequence_results.json --format json
malva_client search "CD4 T cells in brain tissue" --output nl_results.csv --format csv
```

## Quick Start (Python)

```python
from malva_client import MalvaClient

client = MalvaClient()

# Search for a gene and inspect aggregate expression rows
results = client.search("CD3D")
print(results.df.head())

# Search a DNA sequence of at least 24 nt
sequence_results = client.search_sequences("ATCGATCGATCGATCGATCGATCG")
print(sequence_results.df.head())
```

### Per-cell matrices

```python
# First run a normal expression search
result = client.search("SPP1")

# Retrieve positive cells from that search job
cells = client.retrieve_cells(result, include_sample_metadata=False)
sample_id = int(cells.cells.iloc[0]["sample_id"])

# Inspect cells from one real encoded sample ID
sample_cells = client.retrieve_cells(result, sample_ids=[sample_id], include_sample_metadata=True)
cell_ids = sample_cells.get_cell_ids()
cell_df = sample_cells.to_dataframe(include_sample_metadata=True)
print(cell_ids.head())
print(cell_df.head())

# Fetch all cells for that sample once from metadata for denominator analyses
all_cells = client.get_cells_by_metadata(sample_ids=[sample_id])
print(all_cells.head())

# Full-database retrieval is also available; it queues a server-side export.
# all_cells = client.get_cells_by_metadata(include_all_database_cells=True)
```

### Working with Results

```python
# Enrich results with metadata
results.enrich_with_metadata()

# Filter and aggregate
filtered = results.filter_by(disease='normal', organ='brain')
by_cell_type = results.aggregate_by('cell_type', agg_func='mean')
print(filtered.df.head())
print(by_cell_type.head())
```

See the [tutorials](https://malva-client.readthedocs.io/en/latest/tutorials.html) for coverage analysis, dataset discovery, cell-level searches, and more.

### Streaming expression to AnnData

```python
import anndata as ad
from malva_client import MalvaClient

client = MalvaClient()
sample_id = 30064771072  # replace with a sample id from client.search_samples()

# Terminal conversion: obs names are "sample:cell", var names are genes,
# X is CSR float32. Each server batch is folded into a sparse block as it
# arrives, so peak memory stays close to one copy of the expression data.
adata = client.stream_cell_expression(
    sample_ids=[sample_id],
    gene_ids=["SPP1", "CD3D"],
).to_anndata()
print(adata)

# Chunked streaming: yield a chunk once the accumulated batches hold at
# least 5,000 cells. A chunk carries every gene discovered so far, and a
# gene keeps its column once it appears.
chunks = client.stream_cell_expression(
    sample_ids=[sample_id],
).iter_anndata(chunk_obs=5_000)
for chunk in chunks:
    print(chunk.shape)  # e.g. (5000, 812), (5000, 830), (1200, 830)

# A gene first seen after a chunk was emitted is missing from that chunk.
# Combine chunks with an outer join when you need one rectangular matrix.
combined = ad.concat(
    client.stream_cell_expression(sample_ids=[sample_id]).iter_anndata(chunk_obs=5_000),
    join="outer",
    fill_value=0,
)
print(combined)
```

Per-cell metadata is joined into `obs` for every conversion method. With an
explicit cell list, `enrich()` fetches it once from the metadata service:

```python
cells = [(sample_id, cell_id)]  # replace with your (sample_id, cell_id) pairs

adata = client.stream_cell_expression(
    cells=cells,
    gene_ids=["SPP1", "CD3D"],
).enrich().to_anndata()
print(adata.obs[["cell_type", "total_counts"]].head())
```

Pass `cell_metadata=frame` to supply your own table keyed by `sample_id` and
`cell_id`, or `include_cell_metadata=False` to skip both the fetch and the join.


## Indexing Your Own Data

For local indexing and quantification, see [Malva Tools](https://github.com/malva-bio/malva) (`malva` CLI).

## Citation

If you use Malva in your research, please cite our [paper in Nature](https://www.nature.com/articles/s41586-026-10975-w):

```latex
@article{LenPerin2026,
  title = {Ultrafast and reference-free sequence discovery in single-cell data},
  ISSN = {1476-4687},
  url = {http://dx.doi.org/10.1038/s41586-026-10975-w},
  DOI = {10.1038/s41586-026-10975-w},
  journal = {Nature},
  publisher = {Springer Science and Business Media LLC},
  author = {León-Periñán,  Daniel and Karaiskos,  Nikos and Rajewsky,  Nikolaus},
  year = {2026},
  month = Aug 
}
```

## License

The Clear BSD License - Copyright (c) 2025-2026 Malva
