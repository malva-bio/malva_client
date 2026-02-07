Client
======

The :class:`~malva_client.client.MalvaClient` class is the main entry point
for interacting with the Malva API.  It handles authentication, search
requests, result polling, coverage analysis, dataset discovery, and sample
downloads.

MalvaClient
-----------

.. autoclass:: malva_client.client.MalvaClient
   :no-members:

   .. automethod:: __init__

   .. rubric:: Search

   .. automethod:: search
   .. automethod:: search_sequence
   .. automethod:: search_sequences
   .. automethod:: search_genes
   .. automethod:: search_cells

   .. rubric:: Job Management

   .. automethod:: submit_search
   .. automethod:: get_job_status
   .. automethod:: get_job_results
   .. automethod:: wait_for_job
   .. automethod:: list_jobs
   .. automethod:: cancel_job

   .. rubric:: Coverage Analysis

   .. automethod:: get_coverage
   .. automethod:: get_sequence_coverage
   .. automethod:: get_coverage_data
   .. automethod:: download_coverage_wig
   .. automethod:: get_coverage_filter_options

   .. rubric:: Dataset Discovery

   .. automethod:: get_datasets_hierarchy
   .. automethod:: get_dataset_studies
   .. automethod:: get_study_samples
   .. automethod:: get_sample_details
   .. automethod:: get_overview_stats

   .. rubric:: Samples

   .. automethod:: get_samples
   .. automethod:: search_samples
   .. automethod:: get_studies
   .. automethod:: download_sample
   .. automethod:: check_sample_availability
   .. automethod:: get_sample_metadata
   .. automethod:: get_available_filters
   .. automethod:: get_filter_values

   .. rubric:: Account

   .. automethod:: get_quota_status
   .. automethod:: is_authenticated
   .. automethod:: get_database_stats

Convenience Functions
---------------------

These module-level functions create a temporary client and run a single
search.  Useful for quick one-off queries.

.. autofunction:: malva_client.client.search_gene

.. autofunction:: malva_client.client.search_sequence
