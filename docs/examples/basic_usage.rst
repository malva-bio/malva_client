Basic Usage Examples
====================

Quick Start
-----------

.. code-block:: python

    from malva_client import MalvaClient
    
    client = MalvaClient("https://malva.mdc-berlin.net", "YOUR_TOKEN")
    results = client.search("BRCA1")

Search Examples
---------------

.. code-block:: python

    # Gene search
    results = client.search("TP53")
    
    # Natural language search
    results = client.search("find cells with neurodegeneration markers")
    
    # DNA sequence search
    results = client.search("ATCGATCGATCG")