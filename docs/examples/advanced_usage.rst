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

Coverage Analysis
-----------------

.. code-block:: python

    # Get genomic coverage
    coverage = client.get_coverage("chr1", 1000000, 2000000)
    
    # Get sequence coverage
    seq_coverage = client.get_sequence_coverage("ATCGATCG")

Sample Downloads
----------------

.. code-block:: python

    # Download sample data
    sample = client.download_sample("sample-uuid")
    
    # Or save to file
    client.download_sample("sample-uuid", "output.h5ad")