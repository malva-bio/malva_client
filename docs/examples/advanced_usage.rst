Advanced Usage Examples
========================

Async Search Operations
------------------------

.. code-block:: python

    # Submit search without waiting
    job_id = client.submit_search("BRCA1")

    # Check status later
    status = client.get_job_status(job_id)

    # Get results when ready
    if status['status'] == 'completed':
        results = client.get_job_results(job_id)

    # Or wait for completion
    results = client.wait_for_job(job_id)

    # List recent jobs
    jobs = client.list_jobs(limit=10)

    # Cancel a running job
    client.cancel_job(job_id)

Coverage Analysis
-----------------

.. code-block:: python

    # Get genomic region coverage
    coverage = client.get_coverage("chr1", 1000000, 2000000, strand='+')

    # Convert to DataFrame
    df = coverage.to_dataframe()

    # Plot coverage across cell types
    coverage.plot()

    # Get sequence coverage
    seq_coverage = client.get_sequence_coverage("ATCGATCG", sequence_name="my_seq")

    # Download as WIG file
    coverage.download_wig("output.wig")

    # Get filter options for a coverage result
    filters = coverage.get_filter_options()

Dataset Discovery
-----------------

.. code-block:: python

    # Get full dataset hierarchy
    hierarchy = client.get_datasets_hierarchy()

    # Browse studies in a dataset
    studies = client.get_dataset_studies("dataset_id")

    # Browse samples in a study
    samples = client.get_study_samples("dataset_id", "study_name")

    # Get detailed sample metadata
    details = client.get_sample_details("sample-uuid")

    # Search for filter values
    organs = client.get_filter_values("organ", search="brain")

    # Get platform-wide statistics
    stats = client.get_overview_stats()

Sample Downloads
----------------

.. code-block:: python

    # Download sample data
    sample = client.download_sample("sample-uuid")

    # Or save to file
    client.download_sample("sample-uuid", "output.h5ad")
