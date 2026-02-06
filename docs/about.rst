About Malva
===========

Malva is a genomic search platform that provides access to harmonized single-cell and spatial transcriptomics data. The platform enables researchers to search across thousands of samples for gene expression patterns, sequences, and biological signatures.

Key Features
------------

.. grid:: 2

   .. grid-item-card:: Large-Scale Data
      :text-align: center

      Access to over 7,000 single-cell samples across diverse tissues and conditions.

   .. grid-item-card:: Real-time Search
      :text-align: center

      Instantaneous results without reprocessing raw data, powered by a k-mer index.

   .. grid-item-card:: Multiple Query Types
      :text-align: center

      Gene symbols, natural language queries, DNA sequences, and batch search.

   .. grid-item-card:: Harmonized Metadata
      :text-align: center

      Standardized metadata across studies for consistent cross-dataset analysis.

Use Cases
---------

.. grid:: 3

   .. grid-item-card:: Expression & Biomarkers
      :text-align: center

      Gene expression analysis, biomarker discovery, and cell-type characterization across studies.

   .. grid-item-card:: Sequence Events
      :text-align: center

      Circular RNA detection, splice-junction discovery, SNV screening, and fusion transcript identification.

   .. grid-item-card:: Pathogen & Rare Events
      :text-align: center

      Viral/bacterial transcript quantification and low-frequency sequence detection missed in individual studies.

Data Sources
------------

Malva integrates data from:

- Human Cell Atlas datasets
- Published single-cell studies
- Spatial transcriptomics experiments
- Disease-focused research cohorts

All data is processed through standardized pipelines with consistent quality control and metadata annotation.

Indexing Your Own Data
----------------------

For indexing and quantifying your own data locally, use **Malva Tools** (the ``malva`` command-line suite), not the Malva Client. Malva Tools provides local indexing, quantification, and analysis capabilities:

.. code-block:: bash

   # Index a reference transcriptome
   malva index --fasta transcriptome.fa --output my_index

   # Quantify single-cell data against the index
   malva quant --index my_index --reads sample_R2.fastq.gz --output counts.h5ad

Malva Tools documentation: `malva.readthedocs.io <https://malva.readthedocs.io>`_

Malva Tools source code: `github.com/malva-bio/malva <https://github.com/malva-bio/malva>`_

Contact
-------

For questions or support, contact the Malva team at hello@malva.bio.

More information available at `malva.bio <https://malva.bio>`_.
