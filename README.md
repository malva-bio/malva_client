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
malva_client config --server https://malva.bio --token YOUR_API_TOKEN
```

## Quick Start (CLI)

```bash
malva_client search "CD3D" --output results.csv
malva_client search "ATCGATCGATCGATCGATCGATCG" --format json
malva_client search "CD4 T cells in brain tissue"
```

## Quick Start (Python)

```python
from malva_client import MalvaClient

client = MalvaClient("https://malva.bio", "YOUR_API_TOKEN")

# Search for genes, sequences, or natural language queries
results = client.search("CD3D")

# Enrich and visualize
results.enrich_with_metadata()
fig = results.plot_expression_summary("cell_type")

# Filter and aggregate
filtered = results.filter_by(disease='normal', organ='brain')
```

See the [tutorials](https://malva-client.readthedocs.io/en/latest/tutorials.html) for coverage analysis, dataset discovery, cell-level searches, and more.

## Indexing Your Own Data

For local indexing and quantification, see [Malva Tools](https://github.com/malva-bio/malva) (`malva` CLI).

## License

The Clear BSD License - Copyright (c) 2025-2026 Malva
