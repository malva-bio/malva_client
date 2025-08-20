# Malva Client

Python client Malva Index search platform.

## Installation

```bash
pip install malva-client
```

## Quick Start

### Python API

```python
from malva_client import MalvaClient

# Initialize client
client = MalvaClient('https://malva.mdc-berlin.net')

# Login (will prompt for credentials)
client.login()

# Search for genes
results = client.search('BRCA1 expression in breast cancer')
print(results.head())

# Search with specific query type
results = client.search('ATCGATCGATCG', query_type='dna')

# Async search (don't wait for completion)
job_id = client.submit_search('CD4 T cells')
# ... do other work ...
results = client.get_results(job_id)
```

### Command Line Interface

```bash
# Login
malva-client login

# Search
malva-client search "BRCA1 expression"

# Search with output to file
malva-client search "CD4 T cells" --output results.csv

# Check search status
malva-client status <job_id>

# Get results from previous search
malva-client results <job_id> --output results.xlsx

# View search history
malva-client history

# Check quota
malva-client quota
```

## Authentication

The client supports multiple authentication methods:

1. **Interactive login**: Use `client.login()` or `malva-client login` to authenticate interactively
2. **API token**: Pass token directly to client: `MalvaClient(url, api_token='your_token')`
3. **Environment variable**: Set `MALVA_API_TOKEN` environment variable

Credentials are securely stored using the system keyring.

## Configuration

Configuration can be set via:

1. **Configuration file**: `~/.malva_config.json`
2. **Environment variables**: `MALVA_SERVER_URL`, `MALVA_API_TOKEN`, `MALVA_VERIFY_SSL`
3. **Command line arguments**: `--server`, `--token`, `--no-ssl-verify`

Example config file:
```json
{
  "server_url": "https://malva.mdc-berlin.net",
  "verify_ssl": true,
  "timeout": 30
}
```

## Query Types

The client supports three types of queries:

- **Natural language**: `"BRCA1 expression in breast cancer"`
- **Gene names**: `"BRCA1"` 
- **DNA sequences**: `"ATCGATCGATCG"`

Query types are auto-detected but can be specified explicitly.

## Output Formats

Results can be exported in multiple formats:

- **CSV**: `--format csv` or `--output file.csv`
- **Excel**: `--format excel` or `--output file.xlsx`  
- **JSON**: `--format json` or `--output file.json`
- **Table**: `--format table` (default for CLI)

## Python API Examples

### Basic Usage

```python
from malva_client import MalvaClient
import pandas as pd

# Initialize and login
client = MalvaClient('https://malva.mdc-berlin.net')
client.login()

# Simple search
results = client.search('BRCA1')
print(f"Found {len(results)} cells")
print(results.groupby('cell_type').size())
```

### Advanced Usage

```python
# Submit multiple searches asynchronously
job_ids = []
queries = ['BRCA1', 'TP53', 'MYC']

for query in queries:
    job_id = client.submit_search(query)
    job_ids.append(job_id)

# Wait for all to complete and collect results
all_results = []
for job_id in job_ids:
    # Wait for completion
    while True:
        status = client.get_status(job_id)
        if status.status == 'completed':
            break
        elif status.status == 'error':
            print(f"Job {job_id} failed: {status.error}")
            continue
        time.sleep(5)
    
    # Get results
    results = client.get_results(job_id)
    all_results.append(results)

# Combine results
combined = pd.concat(all_results, ignore_index=True)
```

### Context Manager

```python
# Automatically cleanup connections
with MalvaClient('https://malva.mdc-berlin.net') as client:
    client.login()
    results = client.search('BRCA1')
    # Connection automatically closed
```

## CLI Examples

### Basic Searches

```bash
# Simple search
malva search "BRCA1 expression"

# DNA sequence search
malva search "ATCGATCGATCG" --type dna

# Save results to file
malva search "CD4 T cells" --output cd4_results.csv
```

### Asynchronous Searches

```bash
# Submit search without waiting
malva search "BRCA1" --no-wait
# Job submitted with ID: ABC123...

# Check status later
malva status ABC123

# Get results when ready
malva results ABC123 --output brca1_results.xlsx
```

### Managing History

```bash
# View recent searches
malva history --limit 20

# Check quota status
malva quota

# Show configuration
malva config-show
```

## Error Handling

```python
from malva_client import MalvaClient, SearchError, AuthenticationError

try:
    client = MalvaClient('https://malva.mdc-berlin.net')
    results = client.search('BRCA1')
except AuthenticationError:
    print("Please login first")
    client.login()
except SearchError as e:
    print(f"Search failed: {e}")
except Exception as e:
    print(f"Unexpected error: {e}")
```

## Requirements

- Python 3.8+
- pandas
- requests
- click (for CLI)
- rich (for CLI formatting)
- keyring (for secure credential storage)

## License

Copyright 2025 Malva