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

## Indexing Your Own Data

For local indexing and quantification, see [Malva Tools](https://github.com/malva-bio/malva) (`malva` CLI).

## Citation

If you use Malva in your research, please cite:

> [TBA]

## License

The Clear BSD License - Copyright (c) 2025-2026 Malva
