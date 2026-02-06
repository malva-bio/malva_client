Malva Client
============

Python client for the `Malva <https://malva.bio>`_ genomic search platform, enabling programmatic access to harmonized single-cell and spatial transcriptomics data.

.. grid:: 2

   .. grid-item-card:: Getting Started
      :link: installation
      :link-type: doc

      Install the package and configure authentication.

   .. grid-item-card:: Tutorials
      :link: tutorials/index
      :link-type: doc

      Step-by-step Jupyter notebooks for common workflows.

   .. grid-item-card:: Examples
      :link: examples/index
      :link-type: doc

      Code examples for search, coverage, and data access.

   .. grid-item-card:: API Reference
      :link: api/modules
      :link-type: doc

      Full reference for all classes and methods.

Key Features
------------

.. grid:: 3

   .. grid-item-card:: Search
      :text-align: center

      Gene symbols, natural language, and DNA sequences.

   .. grid-item-card:: Coverage
      :text-align: center

      Genomic region and sequence coverage across cell types.

   .. grid-item-card:: Dataset Discovery
      :text-align: center

      Browse datasets, studies, samples, and metadata.

Quick Example
-------------

.. code-block:: python

   from malva_client import MalvaClient

   client = MalvaClient("https://malva.bio", "YOUR_API_TOKEN")

   # Gene expression search
   results = client.search("BRCA1")
   results.enrich_with_metadata()
   results.plot_expression_by("cell_type")

   # Genomic coverage
   coverage = client.get_coverage("chr1", 1000000, 2000000)
   coverage.plot()

   # Dataset hierarchy
   hierarchy = client.get_datasets_hierarchy()

.. toctree::
   :maxdepth: 2
   :caption: Guides
   :hidden:

   about
   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: Tutorials
   :hidden:

   tutorials/index

.. toctree::
   :maxdepth: 2
   :caption: Examples
   :hidden:

   examples/index

.. toctree::
   :maxdepth: 2
   :caption: API Reference
   :hidden:

   api/modules

.. toctree::
   :maxdepth: 1
   :caption: Links
   :hidden:

   Malva Platform <https://malva.bio>
   GitHub <https://github.com/malva-bio/malva_client>
   Malva Tools <https://malva.readthedocs.io>
   Changelog <https://github.com/malva-bio/malva_client/blob/main/CHANGELOG.md>
