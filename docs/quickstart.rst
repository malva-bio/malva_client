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

Select a workspace
------------------

A workspace identifies one project index (dataset). Pass it explicitly or set
the environment variables once:

.. code-block:: bash

   export MALVA_WORKSPACE="atlas:human"
   export MALVA_DATASET_ID="human"
   export MALVA_PROVIDER_ID="malva-public-index"

.. code-block:: python

   client = MalvaClient(workspace="atlas:human", dataset_id="human")
   # or later:
   client.set_workspace(workspace="atlas:human", dataset_id="human")

The client sends these as ``X-Malva-Workspace`` / ``X-Malva-Dataset`` /
``X-Malva-Provider`` headers and includes ``dataset_id`` in searches, so the
Explorer resolves the same project/index the Console has selected. Per-call
``dataset_id`` arguments still win over the configured default.

Provider scoping
----------------

Metadata and sample methods accept a ``provider`` argument:

.. code-block:: python

   samples = client.get_samples(provider="malva-public-index")
   studies = client.get_studies(provider="malva-public-index")
   values = client.get_filter_values("organ", provider="malva-public-index")

The value is sent as the ``provider`` query parameter. An explicit argument
wins over the client-level ``provider_id`` configured in the constructor, via
``set_workspace(provider_id=...)`` or ``MALVA_PROVIDER_ID``.

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

Federated search
^^^^^^^^^^^^^^^^

When the Explorer has several providers configured, ``search_federated``
submits the query to every shard and returns one table with the provider
origin as the first columns:

.. code-block:: python

   df = client.search_federated("FOXP3", dataset_id="main")
   print(df[["provider", "provider_id", "sample_id", "cell_type", "rel"]])

``provider`` is the human-readable provider name; ``provider_id`` is the
stable id used to address that shard. Pass ``provider_ids=[...]`` to restrict
the search; requested ids are validated against ``get_workspace()`` before the
job is submitted.

Every result, including an empty successful search, exposes provenance in
``df.attrs``: ``job_id``, ``status``, ``shards``, ``providers``, ``errors``
and ``shard_errors``. Failed shards emit a warning and are listed in
``attrs["errors"]`` instead of being silently dropped. Pass ``strict=True`` to
raise ``MalvaAPIError`` when any shard fails.

Study discovery
^^^^^^^^^^^^^^^

``get_studies`` and ``get_dataset_studies`` fall back to the dataset hierarchy
when a provider's study endpoint returns no rows. In that case ``get_studies``
adds a ``dataset_id`` column so the study can be passed directly to
``get_study_samples(dataset_id, study_name)``. If a hierarchy-derived study
cannot be read through the provider's sample route, the client raises a clear
``MalvaAPIError`` rather than returning a misleading empty list.

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

``enrich()`` attaches sample metadata and returns the same result object, so it
chains naturally. Only sample IDs that are not already available locally are
requested from Malva:

.. code-block:: python

   result = client.search("SPP1")
   result.enrich()

   print(result.available_filter_fields())

   # This may be empty if no matching rows exist in the current index.
   brain = result.filter_by(organ="brain")
   print(brain.to_dataframe().head())

   by_cell_type = result.aggregate_by("cell_type", agg_func="mean")
   print(by_cell_type.head())

The older ``enrich_with_metadata()`` name still works and is an alias.

Convert search results to AnnData
---------------------------------

Use ``to_anndata()`` to get a standard in-memory AnnData object. The layout is
fixed and easy to remember: observations are unique ``sample_id`` ×
``cell_type`` pairs, variables are genes, and ``X`` is a sparse float32 matrix
of relative expression. Sample and cell-type metadata are in ``obs`` and
provenance is in ``uns["malva"]``:

.. code-block:: python

   adata = client.search("SPP1, CD3D").enrich().to_anndata()
   print(adata)
   print(adata.obs.head())
   print(adata.var_names[:10])


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

Stream per-cell expression
--------------------------

For selected samples, ``stream_cell_expression()`` streams sparse non-zero
per-cell values in batches instead of building one large matrix in memory.
Pass a gene panel whenever possible; without one the server has to scan the
full expression store.

.. code-block:: python

   for batch in client.stream_cell_expression(
       sample_ids=[sample_id],
       gene_ids=["SPP1", "CD3D"],
   ):
       print(batch.head())

Each batch is a pandas ``DataFrame`` with ``sample_id``, ``cell_id``,
``store_index``, ``gene``, and ``value`` columns. Only non-zero entries are
returned. Iterate over the batches rather than wrapping large streams in
``list()``.

The returned object is a small ``CellExpressionStream`` wrapper. You can still
loop over it exactly as above, or convert it to AnnData. ``to_anndata()`` reads
the whole stream and builds a cells × genes sparse float32 matrix; observation
names are ``sample_id:cell_id`` and variables are genes:

.. code-block:: python

   stream = client.stream_cell_expression(
       sample_ids=[sample_id],
       gene_ids=["SPP1", "CD3D"],
   )
   adata = stream.to_anndata()
   print(adata)

For large selections, use ``write_h5ad()``. It writes the stream to disk with
bounded memory and can be read lazily afterwards:

.. code-block:: python

   client.stream_cell_expression(
       sample_ids=[sample_id],
       gene_ids=["SPP1", "CD3D"],
   ).write_h5ad("cells.h5ad")

   import anndata as ad
   backed = ad.read_h5ad("cells.h5ad", backed="r")
   print(backed)

If you prefer to build the matrix yourself, ``iter_anndata()`` yields one
AnnData per server batch. Combine the parts with ``anndata.concat``:

.. code-block:: python

   import anndata as ad

   parts = client.stream_cell_expression(
       sample_ids=[sample_id],
       gene_ids=["SPP1", "CD3D"],
   ).iter_anndata()
   combined = ad.concat(parts, join="outer", fill_value=0)
   print(combined)


Get expression for specific cells
--------------------------------

If you already have a list of cells, pass ``(sample_id, cell_id)`` pairs.
Cells come back in the same order:

.. code-block:: python

    cells = sorted_cells  # list of (sample_id, cell_id) pairs
    for batch in client.stream_cell_expression(
        cells=cells,
        gene_ids=["CD3D", "CD3E"],
    ):
        print(batch.head())

This is useful when you have sorted a cell table by expression and want to
fetch the cells in high-to-low order. Leave out ``gene_ids`` to get all
features for each cell.

The server needs a cell-major CSR index in the expression store. If the store
does not have one, the call fails with a clear "not available yet" error
instead of scanning the whole store.

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

Coverage-like Cell Retrieval
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The coverage API currently returns aggregate coverage tracks and does not yet
provide a dedicated endpoint for retrieving every per-cell hit behind a coverage
region. When you need a cell-level approximation for a sequence or genomic
segment, one practical workaround is to submit the sequence as overlapping
search probes in a normal search job, then call ``retrieve_cells()`` on that
job. This returns cells positive for one or more sliding-window probes.

The following is example code and may need to be adjusted for your sequence,
window size, step size, strandedness, and downstream filtering:

.. code-block:: python

   def sliding_windows(sequence, window=48, step=24):
       sequence = sequence.upper().replace(" ", "").replace("\n", "")
       return [
           sequence[i:i + window]
           for i in range(0, max(len(sequence) - window + 1, 0), step)
           if len(sequence[i:i + window]) == window
       ]

   # Example sequence; replace with your region/sequence of interest.
   sequence = (
       "ACAGGCTGAAGGAGTTTTATTTCAAATGGTTGCCTGATGCCTGTGGTT"
       "AACCACAGGCATCAGGCAACCATTTGAAATAAAACTCCTTCAGCCTGT"
   )

   probes = sliding_windows(sequence, window=48, step=24)
   result = client.search_sequences(probes, stranded=False)

   positive_cells = client.retrieve_cells(
       result,
       include_sample_metadata=True,
   ).to_dataframe(include_sample_metadata=True)

   print(positive_cells.head())

If you need denominator cells for fraction-expressing or mean-including-zero
calculations, fetch the metadata-defined cell universe for the samples observed
in ``positive_cells``:

.. code-block:: python

   sample_ids = positive_cells["sample_id"].dropna().astype(int).unique().tolist()
   denominator_cells = client.get_cells_by_metadata(sample_ids=sample_ids)
   print(denominator_cells.head())

Compare Your Clusters to Malva Metacells
-----------------------------------------

If you have your own AnnData object, you can ask which Malva cells look most
like each of your clusters. The data stays local: the client averages cells
inside each cluster and sends only the small cluster × gene matrix, not every
cell.

.. code-block:: python

   import scanpy as sc

   # Normalise and log-transform first (the client warns if X still looks
   # like raw integer counts). Your object is never modified by the client.
   sc.pp.normalize_total(adata, target_sum=1e4)
   sc.pp.log1p(adata)

   matches = client.find_metacells(
       adata,
       groupby="leiden",          # any categorical obs column
       dataset_id="Human_Cell_Atlas",
       top_n=50,
   )
   print(matches.to_dataframe().head())

The result is a small ``MetacellMatchResult`` object with one row per cluster
and ranked metacell. Use ``to_dataframe()`` whenever you want a plain pandas
table. Columns are ``cluster``, ``rank``, ``metacell_id``, ``sample_id``,
``cell_type``, ``n_cells``, ``similarity``, ``distance``
(``1 - similarity``), ``enrichment_score`` (a per-cluster z-score of
similarity across all metacells), and ``percentile`` (the fraction of atlas
metacells with a similarity at or below that neighbour).

``similarity`` is a cosine in CLR (centred-log-ratio) space over a fixed panel
of variable, non-technical genes, not a cosine of raw expression. The server
divides each profile by its panel total, log-transforms it, and centres it, so
sequencing depth, metacell size and most global technical shifts cancel out.
Raw counts are accepted too; the server converts them to proportions before
scoring.

Sample metadata is merged locally when you already have it, so no extra
metadata request is made:

.. code-block:: python

   sample_meta = my_sample_table.set_index("sample_id")  # your own table
   matches = client.find_metacells(
       adata,
       groupby="leiden",
       dataset_id="Human_Cell_Atlas",
       sample_metadata=sample_meta,
   ).enrich()

Genes that are not present in the Malva dataset are ignored, and the client
warns which ones were missing.

Assign one label per cluster
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``assign_labels()`` combines the neighbour hits into one label per cluster.
With CLR-scored results it defaults to reciprocal-rank weights (``1/rank``),
which keep the best single match dominant while still using the whole
neighbourhood; this avoids abundant but merely related cell types outvoting the
top hit. Legacy raw-mean results keep similarity weighting. It returns a small
table with the label, confidence, and number of supporting metacells. Pass
``adata`` to also write the labels (and optionally the confidence) into
``obs``:

.. code-block:: python

   labels = matches.assign_labels(
       adata,
       column="malva_cell_type",
       confidence_column="malva_confidence",
   )
   print(labels)

   print(adata.obs[["leiden", "malva_cell_type", "malva_confidence"]].head())

Use ``weight="similarity"`` to reproduce the previous summed-similarity
behaviour, or ``weight="rank"`` to request reciprocal-rank weighting
explicitly. ``weight="auto"`` is the default.

``compare_clusters_to_metacells()`` remains available as a deprecated alias for
``find_metacells()``.

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
