Quick Start
===========

This guide shows a minimal Malva client workflow that can be run from top to
bottom once your API credentials are configured. For installation instructions
see :doc:`installation`.

Authentication
--------------

Before querying Malva, create an API token at `malva.bio <https://malva.bio>`_
(**Profile > Generate API Token**) and store it with the CLI:

.. code-block:: bash

   malva_client config --server https://malva.mdc-berlin.de --token YOUR_API_TOKEN

Alternatively, set ``MALVA_API_TOKEN`` in your shell before running Python.

Connect
-------

.. code-block:: python

   from malva_client import MalvaClient

   client = MalvaClient()
   print(client.is_authenticated())

Run Searches
------------

Gene search
^^^^^^^^^^^

Search by gene symbol and inspect the aggregate expression table:

.. code-block:: python

   result = client.search("BRCA1")
   df = result.df
   print(df.head())
   print(df.columns.tolist())

The default search result is aggregated by sample and cell type. Common columns
include ``sample_id``, ``cell_type``, ``gene_sequence``, ``rel``, ``exp``,
``pct``, ``raw_kmers``, and ``cell_count``. These expression columns match
the Expression Explorer display modes.

Sequence search
^^^^^^^^^^^^^^^

Search a DNA sequence. The example sequence is at least 24 nt, matching the
index k-mer size.

.. code-block:: python

   sequence = "ATCGATCGATCGATCGATCGATCG"
   sequence_result = client.search_sequences(sequence)
   print(sequence_result.df.head())

Batch search
^^^^^^^^^^^^

Submit multiple genes in one request and compare the returned rows by query:

.. code-block:: python

   batch = client.search_genes(["BRCA1", "TP53"])
   batch_df = batch.df
   print(batch_df[["gene_sequence", "sample_id", "cell_type", "rel"]].head())

Tune K-mer Filters
------------------

Use ``min_kmer_presence``, ``max_kmer_presence``, and ``stranded`` to control
which k-mers participate in sequence matching:

.. code-block:: python

   filtered = client.search_sequences(
       sequence,
       max_kmer_presence=10000,
       stranded=False,
   )
   print(filtered.df.head())

See :doc:`query_parameters` for guidance on choosing these values.

Work with Results
-----------------

Enrich aggregate search results with sample metadata, then filter and aggregate
using pandas-like methods:

.. code-block:: python

   result = client.search("SPP1")
   result.enrich_with_metadata()

   print(result.available_filter_fields())

   # This may be empty if no matching rows exist in the current index.
   brain = result.filter_by(organ="brain")
   print(brain.to_pandas().head())

   by_cell_type = result.aggregate_by("cell_type", agg_func="mean")
   print(by_cell_type.head())


Expression Columns
------------------

Aggregate search results use one row per ``sample_id`` × ``cell_type`` × query.
The main expression columns are:

``rel``
   Relative normalized expression used by the default Explorer view.

``exp``
   Raw aggregate expression value from the search result payload.

``pct``
   Percent of cells positive for the query in that sample × cell-type group.
   Divide by 100 to obtain fraction expressing.

``raw_kmers``
   Mean raw k-mer hit count per expressing cell, without normalization.

``cell_count``
   Number of positive cells in that sample × cell-type group.

Per-cell retrieval is different: ``retrieve_cells()`` returns positive cells,
and its ``value`` column contains raw per-cell expression/k-mer counts when the
server has stored per-cell values for that search. Missing cell × feature
entries are zero.

Retrieve Cells
--------------

Use ``retrieve_cells()`` when you need cell IDs for downstream analysis. Start
with a normal aggregate search, then choose an encoded sample ID from the cells
that were returned.

.. code-block:: python

   result = client.search("SPP1")
   cells = client.retrieve_cells(result, include_sample_metadata=False)

   sample_id = int(cells.cells.iloc[0]["sample_id"])
   sample_cells = client.retrieve_cells(
       result,
       sample_ids=[sample_id],
       include_sample_metadata=True,
   )

   cell_ids = sample_cells.get_cell_ids()
   cell_df = sample_cells.to_dataframe(include_sample_metadata=True)
   print(cell_ids.head())
   print(cell_df.head())

For aggregate searches, the ``value`` column is a positive-cell indicator:
``1`` means the cell was positive for the query. To get the denominator cells
for downstream fraction-expressing or mean-including-zero calculations, fetch
the metadata-defined cell universe independently and reuse it across searches:

.. code-block:: python

   all_cells = client.get_cells_by_metadata(sample_ids=[sample_id])
   print(all_cells[["sample_id", "cell_id", "cell_type", "total_counts"]].head())

Fetching every cell in the database is also supported through
``get_cells_by_metadata(include_all_database_cells=True)``, but it is a large
operation and is intentionally not part of the quick-start workflow.

Coverage Results
----------------

Coverage queries return aggregate tracks, not single-cell coverage. Use
``to_dataframe()`` for the cell-type × position display matrix and
``to_long_dataframe()`` when you need sample-aware downstream filtering:

.. code-block:: python

   coverage = client.get_coverage("chr11", 67435510, 67439682)

   matrix = coverage.to_dataframe()
   per_sample = coverage.to_long_dataframe()

   sample_rows = per_sample[per_sample["sample_id"] == per_sample.iloc[0]["sample_id"]]
   print(matrix.head())
   print(sample_rows.head())

The long table has one row per genomic position × sample × cell type and
includes ``raw_signal``, ``cell_count``, and ``mean_signal``.

Asynchronous Jobs
-----------------

Submit without waiting, poll status, and then fetch results:

.. code-block:: python

   pending = client.search("FOXP3", wait_for_completion=False)
   job_id = pending.job_id

   status = client.get_job_status(job_id)
   print(status["status"])

   completed = client.wait_for_job(job_id, max_wait=600)
   completed_df = completed.df
   print(completed_df.head())

Command-Line Interface
----------------------

The same basic search workflow is available from the terminal:

.. code-block:: bash

   malva_client search "BRCA1" --output results.csv --format csv
   malva_client search "ATCGATCGATCGATCGATCGATCG" --output sequence_results.json --format json
   malva_client quota

Next Steps
----------

Coverage, coexpression, dataset discovery, and sample download workflows depend
on the datasets and samples you want to analyze. See the dedicated tutorials and
API reference for those workflows.
