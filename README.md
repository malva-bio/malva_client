# Malva Client

[![Documentation](https://readthedocs.org/projects/malva-client/badge/?version=latest)](https://malva-client.readthedocs.io)
[![GitHub](https://img.shields.io/badge/GitHub-malva--bio%2Fmalva__client-blue)](https://github.com/malva-bio/malva_client)

Remote client for the [Malva](https://malva.bio) search platform.

## System Requirements

### Software Dependencies
See `requirements.txt` or `pyproject.toml`

### Operating Systems
- Linux (Ubuntu 18.04+, CentOS 7+)
- macOS 10.15+
- Windows 10+

### Hardware Requirements
- Internet connection required
- No special hardware requirements
- Standard desktop computer sufficient for API queries

## Installation

Installation time: <1 minute on a standard desktop computer

```bash
cd <this folder>
pip install .
```

For full functionality including single-cell analysis:
```bash
pip install scanpy
```

## Authentication

The Malva platform provides two complementary interfaces:
1. **Web Interface**: Browser-based access at https://malva.bio for interactive exploration
2. **API/Command Line**: Programmatic access for automated analyses and integration into computational workflows

### Generating Your API Token

To access Malva programmatically, you have two options:

**Option A**: Use a token provided by the Malva team

**Option B**: Generate your own token via ORCID authentication:
1. Navigate to [https://malva.bio](https://malva.bio)
2. Login with your ORCID account
3. Click the menu button in the upper left corner
4. Select "Profile" from the dropdown
5. Click "Generate API Token"
6. Copy the generated token

### Configure the Client

```bash
# Configure client with your API token
malva_client config --server https://malva.bio --token YOUR_API_TOKEN
```

**Note**: The platform is currently in beta phase. Users receive 20 queries per day (extensible upon request).

## Quick Start (CLI)

### Basic Usage Examples

```bash
# Gene search (runtime: 1-5 seconds)
malva_client search "CD3D" --output results.csv

# Sequence search (runtime: 1-5 seconds)
malva_client search "ATCGATCGATCGATCGATCGATCG" --format json

# Circular RNA search (runtime: 1-5 seconds)
# the results from this search query are provided as an example, see cdr1as_results.csv
malva_client search "CDR1as" -o cdr1as_results.csv --format csv

# Natural language query (runtime: 5-10 seconds)
malva_client search "CD4 T cells in brain tissue"
```

### Output Formats

- CSV: `--format csv` (tabular data)
- JSON: `--format json` (structured data)
- Excel: `--format excel` (spreadsheet format)
- Table: `--format table` (default, console display)

### Expected Output

- **Gene searches**: Cell counts and expression levels across cell types and tissues
- **Sequence searches**: Matching cells with pseudocount quantification
- **Results include**: Cell type annotations, tissue origin, disease status, technology metadata

## Quick Start (Python)

### Basic Usage

```python
from malva_client import MalvaClient

# Initialize client
API_TOKEN = "YOUR_API_TOKEN"  # Place your token here
client = MalvaClient("https://malva.bio", API_TOKEN)

# Search for genes
results = client.search("CD3D")
print(results)

# Search for sequences
results = client.search("ATCGATCGATCGCCACATGGACTTGAC")

# Natural language queries
results = client.search("cells expressing markers of neurodegeneration")
```

### Working with Results

```python
# Enrich results with metadata
results.enrich_with_metadata()

# Visualize expression patterns
fig = results.plot_expression_summary("cell_type")

# Filter and aggregate
filtered = results.filter_by(disease='normal', organ='brain')
aggregated = filtered.aggregate_by('cell_type', agg_func='mean')
```

### Coverage Analysis

```python
# Get genomic region coverage
coverage = client.get_coverage("chr1", 1000000, 2000000)
coverage.plot()

# Get sequence coverage
seq_cov = client.get_sequence_coverage("ATCGATCG")
df = seq_cov.to_dataframe()
```

### Dataset Discovery

```python
# Browse datasets and studies
hierarchy = client.get_datasets_hierarchy()
details = client.get_sample_details("sample-uuid")
stats = client.get_overview_stats()
```

### Advanced Features

```python
import dnaio
from malva_client.tools import mask_sequence

# Load sequence from file
with dnaio.open("sequence.fna") as f_in:
    for s in f_in:
        seq = s.sequence
        break

# Optional: mask low-complexity regions (requires BLAST)
# seq = mask_sequence(seq)

# Search with cell-level resolution
results = client.search_cells(seq)
df_cells = results.enrich_with_metadata()

# Download complete sample for downstream analysis
sample_uuid = df_cells['uuid'].unique()[0]
sample = client.download_sample(sample_uuid)
```

## Indexing Your Own Data

For local indexing and quantification, use **Malva Tools** (`malva` CLI), not this client:

```bash
malva index --fasta transcriptome.fa --output my_index
malva quant --index my_index --reads sample_R2.fastq.gz --output counts.h5ad
```

See the [Malva Tools documentation](https://malva.readthedocs.io) and [source code](https://github.com/malva-bio/malva).

## Use Cases

The Malva client enables sequence searches across a harmonized index of >7,000 single-cell and spatial transcriptomics samples:

- **Cross-study comparisons**: Identify expression patterns across experimental conditions
- **Rare event detection**: Find low-frequency sequences missed in individual studies
- **Viral/bacterial detection**: Quantify pathogen transcripts in tissue samples
- **Circular RNA analysis**: Detect back-splicing events
- **Novel junction discovery**: Identify rare splicing events or fusion transcripts

## Expected Runtime

- Individual gene/sequence queries: 1-20 seconds
- Metadata enrichment: 1-5 seconds
- Sample download: 10-60 seconds (depending on sample size)

## License

The Clear BSD License - Copyright (c) 2025-2026 Malva
