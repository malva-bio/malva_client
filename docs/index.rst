Malva Client Documentation
==========================

Python client for the Malva genomic search platform, enabling programmatic access to harmonized single-cell and spatial transcriptomics data.

Installation
------------

.. code-block:: bash

   pip install malva_client

Quick Start
-----------

**Authentication**

Get your API token from malva.mdc-berlin.de → Profile → Generate API Token

.. code-block:: bash

   # Configure client
   malva_client config --server https://malva.mdc-berlin.de --token YOUR_TOKEN

   # Or login interactively
   malva_client login

**Python API**

.. code-block:: python

   from malva_client import MalvaClient

   # Initialize client
   client = MalvaClient("https://malva.mdc-berlin.de", "YOUR_API_TOKEN")

   # Search for genes or sequences
   results = client.search("BRCA1")
   results = client.search("find cells with hallmarks of neurodegeneration")
   results = client.search("ATCGATCGATCG")  # DNA sequence

   # Get coverage data
   coverage = client.get_coverage("chr1", 1000000, 2000000)

   # Download samples
   sample = client.download_sample("sample-uuid")

**Command Line Interface**

.. code-block:: bash

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
   malva_client config --show

Features
--------

- **Gene symbols**: Search for specific genes like "BRCA1", "TP53"
- **Natural language**: Use plain English queries like "find cells with neurodegeneration markers"
- **DNA sequences**: Search for specific sequences like "ATCGATCGATCG"
- **Coverage analysis**: Get expression coverage across genomic regions
- **Sample downloads**: Access and download sample data

Output Formats
--------------

- CSV: ``--format csv``
- JSON: ``--format json``
- Excel: ``--format excel``
- Table: ``--format table`` (default)

Configuration
-------------

Configuration via file (``~/.malva/config.json``), environment variables, or command line:

.. code-block:: bash

   # Environment variables
   export MALVA_API_URL="https://malva.mdc-berlin.de"
   export MALVA_API_TOKEN="your_token"

   # Command line
   malva_client config --server URL --token TOKEN

.. toctree::
   :maxdepth: 2
   :caption: Documentation

   about
   installation
   quickstart
   api/modules
   examples/index

.. toctree::
   :maxdepth: 1
   :caption: Links

   Malva Platform <https://malva.mdc-berlin.de>
   GitHub Repository <https://github.com/malva_bio/malva_client>