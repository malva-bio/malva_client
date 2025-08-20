# Malva Client

Python client for the Malva genomic search platform, enabling programmatic access to harmonized single-cell and spatial transcriptomics data.

## Installation

```bash
pip install malva_client
```

## Quick Start

### Authentication

Get your API token from [malva.mdc-berlin.de](https://malva.mdc-berlin.de) → Profile → Generate API Token

```bash
# Configure client
malva_client config --server https://malva.mdc-berlin.net --token YOUR_API_TOKEN
```

### Python API

```python
from malva_client import MalvaClient

# Initialize client
client = MalvaClient("https://malva.mdc-berlin.net", "YOUR_API_TOKEN")

# Search for genes or sequences
results = client.search("BRCA1")
results = client.search("find cells with hallmarks of neurodegeneration")
results = client.search("ATCGATCGATCG")  # DNA sequence

# Get coverage data
coverage = client.get_coverage("chr1", 1000000, 2000000)

# Download samples
sample = client.download_sample("sample-uuid")
```

### Command Line Interface

```bash
# Search operations
malva_client search "BRCA1 expression"
malva_client search "CD4 T cells" --output results.csv --format csv

# Async operations
malva_client search "FOXP3" --no-wait  # Returns job ID
malva_client status <job_id>
malva_client results <job_id> --output results.xlsx

# Account management
malva_client quota
malva_client history
malva_client config # shows the current configuration
```

## Query Types

- **Gene symbols**: `"BRCA1"`, `"TP53"`
- **Natural language**: `"find cells with hallmarks of neurodegeneration"`
- **DNA sequences**: `"ATCGATCGATCG"`

## Output Formats

- CSV: `--format csv`
- JSON: `--format json`
- Excel: `--format excel`
- Table: `--format table` (default)

## Configuration

Configuration via file (`~/.malva/config.json`), environment variables, or command line:

```bash
# Environment variables
export MALVA_API_URL="https://malva.mdc-berlin.net"
export MALVA_API_TOKEN="your_token"

# Command line
malva_client config --server URL --token TOKEN
```

## Requirements

- Python 3.8+
- pandas, requests, click, rich

## Documentation

Full documentation and examples: [Link to documentation]

## License

Copyright 2025 Malva