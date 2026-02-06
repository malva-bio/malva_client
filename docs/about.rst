About Malva
===========

Malva is a genomic search platform that provides access to harmonized single-cell and spatial transcriptomics data. The platform enables researchers to search across thousands of samples for gene expression patterns, sequences, and biological signatures.

Key Features
------------

.. grid:: 2
   :gutter: 3

   .. grid-item-card::
      :class-card: sd-border-primary sd-rounded-3
      :class-header: sd-bg-primary sd-text-white

      Large-Scale Data
      ^^^
      Access to over 7,000 single-cell samples across diverse tissues,
      conditions, and species — all queryable through a single API.

   .. grid-item-card::
      :class-card: sd-border-success sd-rounded-3
      :class-header: sd-bg-success sd-text-white

      Real-time Search
      ^^^
      Instantaneous results powered by a k-mer index.  No reprocessing of raw
      data — just query and get answers.

   .. grid-item-card::
      :class-card: sd-border-info sd-rounded-3
      :class-header: sd-bg-info sd-text-white

      Multiple Query Types
      ^^^
      Gene symbols, natural language queries, DNA sequences, and batch search.
      From single genes to genome-scale probes.

   .. grid-item-card::
      :class-card: sd-border-warning sd-rounded-3
      :class-header: sd-bg-warning sd-text-dark

      Harmonized Metadata
      ^^^
      Standardized metadata across studies for consistent cross-dataset
      analysis — cell types, organs, diseases, and more.

Use Cases
---------

.. grid:: 3
   :gutter: 2

   .. grid-item-card:: Expression & Biomarkers
      :text-align: center
      :class-card: sd-rounded-3 sd-shadow-sm

      .. raw:: html

         <div style="font-size: 2rem; margin-bottom: 0.5rem;">&#x1F9EC;</div>

      Gene expression analysis, biomarker discovery, and cell-type
      characterization across studies.

   .. grid-item-card:: Sequence Events
      :text-align: center
      :class-card: sd-rounded-3 sd-shadow-sm

      .. raw:: html

         <div style="font-size: 2rem; margin-bottom: 0.5rem;">&#x1F52C;</div>

      Circular RNA detection, splice-junction discovery, SNV screening,
      and fusion transcript identification.

   .. grid-item-card:: Pathogen & Rare Events
      :text-align: center
      :class-card: sd-rounded-3 sd-shadow-sm

      .. raw:: html

         <div style="font-size: 2rem; margin-bottom: 0.5rem;">&#x1F9A0;</div>

      Viral/bacterial transcript quantification and low-frequency sequence
      detection missed in individual studies.

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
