Quick Start
===========

This guide walks you through connecting to the Malva API, running your first
searches, and exploring the results.  For installation instructions see
:doc:`installation`.

Authentication
--------------

Before you can query Malva you need an API token.  Get one from
`malva.bio <https://malva.bio>`_ **> Profile > Generate API Token**, then
store it in one of the following ways:

.. tab-set::

   .. tab-item:: CLI Config (recommended)

      .. code-block:: bash

         malva_client config --server https://malva.mdc-berlin.de --token YOUR_API_TOKEN

   .. tab-item:: Environment Variable

      .. code-block:: bash

         export MALVA_API_TOKEN=YOUR_API_TOKEN

   .. tab-item:: Interactive Login

      .. code-block:: bash

         malva_client login

      Opens a browser window for ORCID authentication.

Connecting
----------

Once your token is in place, create a client instance.  The constructor
reads credentials from your CLI config or environment automatically:

.. code-block:: python

   from malva_client import MalvaClient

   # Picks up token from CLI config / environment
   client = MalvaClient()

   # Or pass credentials explicitly
   client = MalvaClient("https://malva.mdc-berlin.de", "YOUR_API_TOKEN")


Running a Search
----------------

Malva accepts four kinds of queries.  Every search returns a
:class:`~malva_client.models.SearchResult` that wraps a pandas DataFrame,
so you can filter, aggregate, and plot straight away.

.. tab-set::

   .. tab-item:: Gene symbol

      Search by gene name.  Malva resolves the symbol to its transcriptomic
      sequence and returns expression across all indexed samples.

      .. code-block:: python

         result = client.search("BRCA1")
         result.enrich_with_metadata()
         result.plot_expression_summary("cell_type")

   .. tab-item:: Natural language

      Free-text queries interpreted by Malva's query engine.  Describe the
      biology you are looking for.

      .. code-block:: python

         result = client.search(
             "find cells with hallmarks of neurodegeneration"
         )

   .. tab-item:: DNA sequence

      Raw nucleotide sequences up to 500 kb.  Useful for probes, splice
      junctions, viral sequences, or any arbitrary DNA.

      .. code-block:: python

         result = client.search_sequence(
             "ATCGATCGATCGATCGATCGATCG"
         )

   .. tab-item:: Batch

      Query multiple sequences or genes in a single API call.  See
      :ref:`batch-searches` below for details on working with the results.

      .. code-block:: python

         result = client.search_sequences([
             "ATCGATCGATCGATCGATCGATCG",
             "GCTAGCTAGCTAGCTAGCTAGCTA",
         ])

         result = client.search_genes(["BRCA1", "TP53"])


.. _batch-searches:

Batch Searches
--------------

``search_sequences()`` and ``search_genes()`` send all items in **one
request** and return a single :class:`~malva_client.models.SearchResult`.
The underlying DataFrame contains a ``gene_sequence`` column that
identifies which query each row belongs to, so you can inspect results
per sequence or per gene.

Batch sequences
^^^^^^^^^^^^^^^

.. code-block:: python

   seqs = [
       "ATCGATCGATCG" * 5,   # seq_1
       "GCTAGCTAGCTA" * 5,   # seq_2
   ]
   result = client.search_sequences(seqs)
   result.enrich_with_metadata()

   # See which sequences are in the result
   print(result.df["gene_sequence"].unique())

   # Filter to a single sequence
   seq1_only = result.filter_by(gene_sequence="seq_1")
   seq1_only.plot_expression_by("cell_type", limit=10)

   # Compare sequences side by side
   comparison = result.aggregate_by(["gene_sequence", "cell_type"])
   print(comparison.head(10))

.. tip::

   The total nucleotide count across all sequences must be 100 kb or
   less.  If you need to search longer sequences, use
   ``search_sequence()`` individually.

Batch genes
^^^^^^^^^^^

.. code-block:: python

   result = client.search_genes(["BRCA1", "TP53", "MYC"])
   result.enrich_with_metadata()

   # Filter to a single gene
   brca1 = result.filter_by(gene_sequence="BRCA1")
   brca1.plot_expression_summary("organ")

   # Aggregate across all genes
   by_gene = result.aggregate_by("gene_sequence")
   print(by_gene)

.. tip::

   Up to 10 genes can be queried per call.

Tuning Search Parameters
-------------------------

Malva's index uses fixed-length k-mers (k = 24).  Two parameters let you
control the sensitivity/specificity trade-off:

* **window_size** -- number of consecutive k-mers evaluated per sliding window
* **threshold** -- fraction of k-mers in a window that must match (0.0--1.0)

See :doc:`query_parameters` for recommended values per use case.

.. tab-set::

   .. tab-item:: Expression (permissive)

      .. code-block:: python

         result = client.search("BRCA1", window_size=96, threshold=0.55)

   .. tab-item:: Exact match (junctions / SNVs)

      .. code-block:: python

         result = client.search_sequence(
             junction_seq, window_size=24, threshold=1.0
         )

   .. tab-item:: Strand-specific

      .. code-block:: python

         result = client.search_sequence(probe, stranded=True)

Working with Results
--------------------

Every search returns a :class:`~malva_client.models.SearchResult` that
behaves like a pandas DataFrame.  Call ``enrich_with_metadata()`` to pull
in harmonized sample annotations (organ, disease, species, study, etc.),
then filter, aggregate, and plot:

.. code-block:: python

   # Enrich with harmonized metadata
   result.enrich_with_metadata()

   # Filter by metadata fields
   brain = result.filter_by(disease="normal", organ="brain")

   # Aggregate across a category
   by_cell_type = brain.aggregate_by("cell_type", agg_func="mean")

   # Visualize
   fig = result.plot_expression_summary("cell_type")

   # Export
   df = result.to_pandas()
   df.to_csv("results.csv")

Single-Cell Resolution
----------------------

By default, Malva aggregates expression per sample and cell type.  Use
``search_cells()`` to retrieve individual cell-level hits instead:

.. code-block:: python

   cells = client.search_cells("CDR1as")
   cells.enrich_with_metadata()
   df = cells.to_pandas()
   print(df[["cell_type", "organ", "pseudocount"]].head())

Coverage Analysis
-----------------

The coverage API returns k-mer match density across a genomic region or
an arbitrary sequence, broken down by cell type.  The result is a
**position x cell-type matrix** that you can plot, convert to a DataFrame,
or export as a WIG file.

Basic usage
^^^^^^^^^^^

.. code-block:: python

   # Genomic region
   coverage = client.get_coverage("chr1", 1000000, 2000000)
   coverage.plot()

   # Arbitrary sequence
   seq_cov = client.get_sequence_coverage(
       "ATCGATCG" * 10, sequence_name="my_probe"
   )
   df = seq_cov.to_dataframe()
   print(df.head())  # columns = cell types, index = positions

Filtering by metadata
^^^^^^^^^^^^^^^^^^^^^

By default, coverage is aggregated across **all** samples in the
database.  To restrict the analysis to a subset -- for example, only
heart samples or a specific study -- pass ``metadata_filters`` at search
time:

.. code-block:: python

   # Only heart samples
   heart_cov = client.get_coverage(
       "chr1", 1000000, 2000000,
       metadata_filters={"organs": ["Heart"]}
   )
   heart_cov.plot()

   # Only samples from a specific study
   study_cov = client.get_coverage(
       "chr1", 1000000, 2000000,
       metadata_filters={"studies": ["Roussos-Human-10x3pv3"]}
   )

   # Combine filters
   filtered_cov = client.get_coverage(
       "chr1", 1000000, 2000000,
       metadata_filters={
           "organs": ["Brain"],
           "disease": ["Normal"],
       }
   )

.. note::

   Metadata filters are applied **before** cell-type aggregation on the
   server, so only matching samples contribute to the coverage matrix.
   To discover which filter values are available for an existing job, use
   ``coverage.get_filter_options()``.

Plotting specific cell types
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once you have a :class:`~malva_client.models.CoverageResult`, you can
select which cell types to display:

.. code-block:: python

   coverage = client.get_coverage("chr1", 1000000, 2000000)

   # Plot only selected cell types
   coverage.plot(cell_types=["Neuron", "Astrocyte"])

   # Convert to DataFrame for custom analysis
   df = coverage.to_dataframe()
   df[["Neuron", "Astrocyte"]].plot(title="Coverage comparison")

Exporting as WIG
^^^^^^^^^^^^^^^^

Download coverage tracks for viewing in a genome browser:

.. code-block:: python

   coverage.download_wig("coverage.wig")

Dataset Discovery
-----------------

Browse the dataset, study, and sample hierarchy:

.. code-block:: python

   hierarchy = client.get_datasets_hierarchy()

   studies = client.get_dataset_studies("human_cell_atlas")
   samples = client.get_study_samples(
       "human_cell_atlas", "Roussos-Human-10x3pv3"
   )
   details = client.get_sample_details(
       "34f13021-4ea8-4fae-b990-33b4d6442621"
   )

   stats = client.get_overview_stats()

Async Operations
----------------

Searches over large sequences or many genes can take a while.  Instead of
blocking, you can submit a job and retrieve results later.

Single query
^^^^^^^^^^^^

.. code-block:: python

   job_id = client.submit_search("FOXP3")

   # ... do other work ...

   status = client.get_job_status(job_id)
   if status["status"] == "completed":
       result = client.get_job_results(job_id)

   # Or block until the job finishes
   result = client.wait_for_job(job_id, max_wait=600)

Async batch
^^^^^^^^^^^

``search_sequences()`` and ``search_genes()`` accept
``wait_for_completion=False``.  The returned
:class:`~malva_client.models.SearchResult` will have ``status='pending'``
and a ``job_id`` you can poll:

.. code-block:: python

   # Submit a batch of sequences without waiting
   pending = client.search_sequences(
       ["ATCGATCG" * 10, "GCTAGCTA" * 10],
       wait_for_completion=False,
   )
   print(pending.job_id)   # UUID of the pending job
   print(pending.status)   # "pending"

   # Check status any time
   status = client.get_job_status(pending.job_id)
   print(status["status"])  # "pending", "running", or "completed"

   # Block until done and get the full result
   result = client.wait_for_job(pending.job_id)
   result.enrich_with_metadata()
   result.plot_expression_by("cell_type")

The same works for genes:

.. code-block:: python

   pending = client.search_genes(
       ["BRCA1", "TP53", "MYC"],
       wait_for_completion=False,
   )
   result = client.wait_for_job(pending.job_id)

Sample Downloads
----------------

Download a complete single-cell sample for downstream analysis:

.. code-block:: python

   # Save to disk
   path = client.download_sample("sample-uuid", output_path="sample.h5ad")

   # Or load directly into memory as AnnData
   adata = client.download_sample("sample-uuid")

Command-Line Interface
----------------------

All core functionality is also available from the terminal:

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
