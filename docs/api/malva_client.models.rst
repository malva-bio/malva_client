Result Objects
==============

Search and coverage methods return rich result objects that wrap pandas
DataFrames with built-in filtering, aggregation, and plotting.

SearchResult
------------

Returned by :meth:`~malva_client.client.MalvaClient.search`,
:meth:`~malva_client.client.MalvaClient.search_sequence`,
:meth:`~malva_client.client.MalvaClient.search_sequences`, and
:meth:`~malva_client.client.MalvaClient.search_genes`.  Inherits all
DataFrame-like methods from :class:`~malva_client.models.MalvaDataFrame`.

.. autoclass:: malva_client.models.SearchResult
   :no-members:

   .. automethod:: enrich_with_metadata

   .. autoproperty:: job_id
   .. autoproperty:: status
   .. autoproperty:: results
   .. autoproperty:: total_cells

MalvaDataFrame
--------------

Base class shared by :class:`~malva_client.models.SearchResult`.  Provides
filtering, aggregation, and plotting on top of a pandas DataFrame.

.. autoclass:: malva_client.models.MalvaDataFrame
   :no-members:

   .. automethod:: __init__

   .. rubric:: Data Access

   .. autoproperty:: df
   .. automethod:: to_pandas
   .. automethod:: available_fields
   .. automethod:: available_filter_fields
   .. automethod:: field_info
   .. automethod:: unique_values
   .. automethod:: value_counts

   .. rubric:: Filtering and Aggregation

   .. automethod:: filter_by
   .. automethod:: aggregate_by

   .. rubric:: Plotting

   .. automethod:: plot_expression_by
   .. automethod:: plot_expression_summary

SingleCellResult
----------------

Returned by :meth:`~malva_client.client.MalvaClient.search_cells`.
Contains individual cell-level hits rather than aggregated expression.

.. autoclass:: malva_client.models.SingleCellResult
   :no-members:

   .. automethod:: __init__

   .. rubric:: Properties

   .. autoproperty:: is_completed
   .. autoproperty:: is_pending
   .. autoproperty:: has_results
   .. autoproperty:: cell_count
   .. autoproperty:: sample_count
   .. autoproperty:: query_gene

   .. rubric:: Data Access

   .. automethod:: to_dataframe
   .. automethod:: get_cell_ids
   .. automethod:: get_expression_values
   .. automethod:: get_sample_ids
   .. automethod:: get_expression_stats

   .. rubric:: Filtering

   .. automethod:: filter_by_expression
   .. automethod:: filter_by_samples
   .. automethod:: get_top_expressing_cells

   .. rubric:: Aggregation and Export

   .. automethod:: aggregate_by_sample
   .. automethod:: enrich_with_metadata
   .. automethod:: save_to_csv

CoverageResult
--------------

Returned by :meth:`~malva_client.client.MalvaClient.get_coverage` and
:meth:`~malva_client.client.MalvaClient.get_sequence_coverage`.  Contains
a position-by-cell-type coverage matrix.

.. autoclass:: malva_client.models.CoverageResult
   :no-members:

   .. automethod:: __init__
   .. automethod:: to_dataframe
   .. automethod:: plot
   .. automethod:: download_wig
   .. automethod:: get_filter_options
