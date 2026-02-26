Models
======

SearchResult
------------

.. autoclass:: malva_client.models.SearchResult
   :no-members:
   :show-inheritance:

   .. automethod:: __init__
   .. autoproperty:: job_id
   .. autoproperty:: status
   .. autoproperty:: results
   .. autoproperty:: total_cells
   .. automethod:: enrich_with_metadata

MalvaDataFrame
--------------

.. autoclass:: malva_client.models.MalvaDataFrame
   :no-members:
   :show-inheritance:

   .. automethod:: __init__
   .. autoproperty:: df
   .. automethod:: filter_by
   .. automethod:: aggregate_by
   .. automethod:: plot_expression_by
   .. automethod:: plot_expression_summary
   .. automethod:: to_pandas
   .. automethod:: available_fields
   .. automethod:: available_filter_fields
   .. automethod:: field_info
   .. automethod:: unique_values
   .. automethod:: value_counts

SingleCellResult
----------------

.. autoclass:: malva_client.models.SingleCellResult
   :no-members:
   :show-inheritance:

   .. automethod:: __init__
   .. autoproperty:: is_completed
   .. autoproperty:: is_pending
   .. autoproperty:: has_results
   .. autoproperty:: cell_count
   .. autoproperty:: sample_count
   .. autoproperty:: query_gene
   .. automethod:: to_dataframe
   .. automethod:: get_cell_ids
   .. automethod:: get_expression_values
   .. automethod:: get_sample_ids
   .. automethod:: filter_by_expression
   .. automethod:: filter_by_samples
   .. automethod:: get_top_expressing_cells
   .. automethod:: aggregate_by_sample
   .. automethod:: get_expression_stats
   .. automethod:: enrich_with_metadata
   .. automethod:: save_to_csv

CoverageResult
--------------

.. autoclass:: malva_client.models.CoverageResult
   :no-members:
   :show-inheritance:

   .. automethod:: __init__
   .. automethod:: to_dataframe
   .. automethod:: get_filter_options
   .. automethod:: download_wig
   .. automethod:: plot

CoexpressionResult
------------------

.. autoclass:: malva_client.models.CoexpressionResult
   :no-members:
   :show-inheritance:

   .. automethod:: __init__
   .. automethod:: genes_to_dataframe
   .. automethod:: scores_to_dataframe
   .. automethod:: umap_to_dataframe
   .. automethod:: go_to_dataframe
   .. automethod:: cell_type_enrichment_to_dataframe
   .. automethod:: tissue_breakdown_to_dataframe
   .. automethod:: get_top_genes
   .. automethod:: plot_umap
   .. automethod:: plot_top_genes
   .. automethod:: plot_go_enrichment

UMAPCoordinates
---------------

.. autoclass:: malva_client.models.UMAPCoordinates
   :no-members:
   :show-inheritance:

   .. automethod:: __init__
   .. automethod:: to_dataframe
   .. automethod:: plot
