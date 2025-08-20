Quick Start
===========

Authentication
--------------

Get your API token from malva.mdc-berlin.de → Profile → Generate API Token

.. code-block:: bash

   # Configure client
   malva_client config --server https://malva.mdc-berlin.net --token YOUR_API_TOKEN

   # Or login interactively
   malva_client login

Python API
----------

.. code-block:: python

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

Command Line Interface
----------------------

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