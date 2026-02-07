Client
======

.. autoclass:: malva_client.client.MalvaClient
   :no-members:
   :show-inheritance:

   .. automethod:: __init__

   .. rubric:: Search

   .. automethod:: search
   .. automethod:: search_cells
   .. automethod:: search_sequence
   .. automethod:: search_sequences
   .. automethod:: search_genes

   .. rubric:: Job Management

   .. automethod:: submit_search
   .. automethod:: get_job_status
   .. automethod:: get_job_results
   .. automethod:: wait_for_job
   .. automethod:: cancel_job
   .. automethod:: list_jobs

   .. rubric:: Coverage Analysis

   .. automethod:: get_coverage
   .. automethod:: get_sequence_coverage
   .. automethod:: get_coverage_data
   .. automethod:: download_coverage_wig
   .. automethod:: get_coverage_filter_options

   .. rubric:: Samples & Datasets

   .. automethod:: get_samples
   .. automethod:: search_samples
   .. automethod:: get_studies
   .. automethod:: download_sample
   .. automethod:: check_sample_availability
   .. automethod:: get_sample_metadata
   .. automethod:: get_sample_details

   .. rubric:: Dataset Discovery

   .. automethod:: get_datasets_hierarchy
   .. automethod:: get_dataset_studies
   .. automethod:: get_study_samples
   .. automethod:: get_filter_values
   .. automethod:: get_overview_stats

   .. rubric:: Utilities

   .. automethod:: get_quota_status
   .. automethod:: is_authenticated
   .. automethod:: get_database_stats
   .. automethod:: get_available_filters
   .. automethod:: print_dict_summary

Convenience Functions
---------------------

.. autofunction:: malva_client.client.search_gene

.. autofunction:: malva_client.client.search_sequence
