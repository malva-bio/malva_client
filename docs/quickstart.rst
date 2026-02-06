Quick Start
===========

Authentication
--------------

Get your API token from `malva.bio <https://malva.bio>`_ → Profile → Generate API Token.

.. code-block:: bash

   # Configure via CLI
   malva_client config --server https://malva.bio --token YOUR_API_TOKEN

   # Or set an environment variable
   export MALVA_API_TOKEN=YOUR_API_TOKEN

   # Or login interactively (opens browser for ORCID authentication)
   malva_client login

Connecting
----------

.. code-block:: python

   from malva_client import MalvaClient

   # Reads token from CLI config / environment automatically
   client = MalvaClient()

   # Or pass credentials explicitly
   client = MalvaClient("https://malva.bio", "YOUR_API_TOKEN")

Search Query Types
------------------

.. grid:: 2

   .. grid-item-card:: Gene Symbol
      :text-align: center

      Search by gene name — ``"BRCA1"``, ``"TP53"``

   .. grid-item-card:: Natural Language
      :text-align: center

      Free-text queries — ``"CD4 T cells in brain tissue"``

   .. grid-item-card:: DNA Sequence
      :text-align: center

      Raw nucleotide sequences up to 500 kb

   .. grid-item-card:: Batch Search
      :text-align: center

      Multiple sequences in a single call via ``search_sequences()``

Gene Search
^^^^^^^^^^^

.. code-block:: python

   results = client.search("BRCA1")
   results.enrich_with_metadata()
   results.plot_expression_summary("cell_type")

Natural Language Search
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   results = client.search("find cells with hallmarks of neurodegeneration")

Sequence Search
^^^^^^^^^^^^^^^

.. code-block:: python

   results = client.search_sequence("ATCGATCGATCGATCGATCGATCG")

Multi-gene Search
^^^^^^^^^^^^^^^^^

.. code-block:: python

   results = client.search_genes(["BRCA1", "TP53"])

Batch Sequence Search
^^^^^^^^^^^^^^^^^^^^^

Each sequence counts as one API query.

.. code-block:: python

   seqs = [
       "ATCGATCGATCGATCGATCGATCG",
       "GCTAGCTAGCTAGCTAGCTAGCTA",
   ]
   batch = client.search_sequences(seqs)

   for seq, result in batch.items():
       print(seq[:20], "→", result)

Tuning Search Parameters
-------------------------

Malva's index uses fixed-length k-mers (k = 24).  Two parameters let you
control the sensitivity/specificity trade-off:

* **window_size** — number of consecutive k-mers evaluated per sliding window
* **threshold** — fraction of k-mers in a window that must match (0.0–1.0)

See :doc:`query_parameters` for recommended values per use case.

.. code-block:: python

   # Transcript expression (default-like, permissive)
   results = client.search("BRCA1", window_size=96, threshold=0.55)

   # Splice junction / exact match
   results = client.search_sequence(junction_seq, window_size=24, threshold=1.0)

   # Strand-specific search
   results = client.search_sequence(probe, stranded=True)

Working with Results
--------------------

.. code-block:: python

   # Enrich with harmonized metadata
   results.enrich_with_metadata()

   # Filter by metadata fields
   filtered = results.filter_by(disease="normal", organ="brain")

   # Aggregate across a category
   aggregated = filtered.aggregate_by("cell_type", agg_func="mean")

   # Visualize
   fig = results.plot_expression_summary("cell_type")

   # Export
   df = results.to_pandas()
   df.to_csv("results.csv")

Single-Cell Resolution
----------------------

Use ``search_cells()`` to retrieve individual cell-level hits instead of
aggregated expression.

.. code-block:: python

   cells = client.search_cells("CDR1as")
   cells.enrich_with_metadata()
   df = cells.to_pandas()
   print(df[["cell_type", "organ", "pseudocount"]].head())

Coverage Analysis
-----------------

Visualize k-mer coverage across a genomic region or arbitrary sequence.

.. code-block:: python

   # Genomic region
   coverage = client.get_coverage("chr1", 1000000, 2000000)
   coverage.plot()

   # Arbitrary sequence
   seq_cov = client.get_sequence_coverage("ATCGATCG" * 10, sequence_name="my_probe")
   df = seq_cov.to_dataframe()

   # Download as WIG for genome browser
   coverage.download_wig("coverage.wig")

Dataset Discovery
-----------------

Browse the dataset → study → sample hierarchy.

.. code-block:: python

   hierarchy = client.get_datasets_hierarchy()

   studies = client.get_dataset_studies("human_cell_atlas")
   samples = client.get_study_samples("human_cell_atlas", "Smith2023")
   details = client.get_sample_details("sample-uuid")

   stats = client.get_overview_stats()

Async Operations
----------------

Submit a search without blocking, then retrieve results later.

.. code-block:: python

   job_id = client.submit_search("FOXP3")

   # … do other work …

   status = client.get_job_status(job_id)
   if status["status"] == "completed":
       results = client.get_job_results(job_id)

   # Or wait for completion
   results = client.wait_for_job(job_id, max_wait=600)

Sample Downloads
----------------

Download a complete single-cell sample for downstream analysis.

.. code-block:: python

   # Save to disk
   path = client.download_sample("sample-uuid", output_path="sample.h5ad")

   # Or load directly into memory as AnnData
   adata = client.download_sample("sample-uuid")

Command-Line Interface
----------------------

.. code-block:: bash

   # Basic search
   malva_client search "BRCA1" --output results.csv --format csv

   # With tuning parameters
   malva_client search "ATCGATCG..." --window-size 24 --threshold 1.0 --stranded

   # Async workflow
   malva_client search "FOXP3" --no-wait
   malva_client status <job_id>
   malva_client results <job_id> --output results.xlsx --format excel

   # Account management
   malva_client quota
   malva_client history
   malva_client config --show
