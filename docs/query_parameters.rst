Query Parameters
================

How Search Works
----------------

Malva indexes every sample with fixed-length **k-mers** (k = 24 nucleotides).
When you submit a query sequence, the engine:

1. Extracts all k-mers from the query.
2. Looks up each k-mer in the index and identifies which cells contain it.
3. Aggregates per-cell hit counts and normalises them into expression scores.

The three user-facing parameters let you filter which k-mers participate in
the search and whether both strands are considered.

Parameters
----------

.. grid:: 3
   :gutter: 3

   .. grid-item-card::
      :class-card: sd-border-primary sd-rounded-3
      :class-header: sd-bg-primary sd-text-white

      min_kmer_presence
      ^^^
      Exclude k-mers that appear in **fewer** than this many cells across
      the entire database.

      * ``0`` (default) — no lower filter; all k-mers used
      * ``10`` – ``100`` — removes k-mers that are extremely rare
        (likely sequencing errors or ultra-low-coverage regions)

   .. grid-item-card::
      :class-card: sd-border-success sd-rounded-3
      :class-header: sd-bg-success sd-text-white

      max_kmer_presence
      ^^^
      Exclude k-mers that appear in **more** than this many cells.

      * ``100000`` (default) — retains nearly all k-mers
      * ``10000`` – ``50000`` — removes highly repetitive / ubiquitous
        k-mers that would inflate scores for unrelated cell types

   .. grid-item-card::
      :class-card: sd-border-info sd-rounded-3
      :class-header: sd-bg-info sd-text-white

      Strandedness
      ^^^
      By default the forward strand only is searched.  Set
      ``stranded=False`` to include reverse-complement k-mers as well.

      Useful for ambiguous queries or when the orientation of your
      sequence is not known.

Recommended Settings
--------------------

.. list-table::
   :header-rows: 1
   :widths: 30 20 20 30
   :class: sd-table-hover

   * - Use Case
     - min_kmer_presence
     - max_kmer_presence
     - Notes
   * - Transcript expression
     - ``0``
     - ``100000``
     - Defaults are fine for most genes
   * - Splice junction detection
     - ``0``
     - ``10000``
     - Junction k-mers should be moderately rare
   * - Viral / pathogen sequences
     - ``0``
     - ``5000``
     - Pathogen k-mers uncommon in human transcriptomes
   * - Repeat / low-complexity regions
     - ``10``
     - ``5000``
     - Discard both error k-mers and repeat elements
   * - SNV / short variant detection
     - ``0``
     - ``50000``
     - Variant k-mers rare; exclude highly repetitive flanking k-mers

Probe Design Guidelines
-----------------------

.. tab-set::

   .. tab-item:: Transcript Expression

      Use the full coding sequence or a representative exon.  The
      default parameters work well.

      .. code-block:: python

         results = client.search("BRCA1")

         # Or with explicit filtering
         results = client.search_sequences(
             "ATCGATCGATCG" * 20,
             max_kmer_presence=50000,
         )

   .. tab-item:: Splice Junction

      Design a probe of ~48 nt centred on the junction (24 nt from each
      exon).  Reducing ``max_kmer_presence`` ensures only junction-specific
      k-mers contribute.

      .. code-block:: python

         junction = "ACGTACGT" * 6  # 48 nt spanning junction
         results = client.search_sequences(
             junction,
             max_kmer_presence=10000,
         )

   .. tab-item:: Circular RNA (BSJ)

      Same principle: design a probe that spans the back-splice junction.

      .. code-block:: python

         bsj_probe = "CTAG" * 12  # 48 nt across BSJ
         results = client.search_sequences(
             bsj_probe,
             max_kmer_presence=10000,
         )

   .. tab-item:: SNV Detection

      Centre a 48 nt probe on the variant position.

      .. code-block:: python

         snv_probe = ref_context[:24] + alt_base + ref_context[25:48]
         results = client.search_sequences(
             snv_probe,
             max_kmer_presence=50000,
         )

.. warning::

   **3'-biased protocols** — Many scRNA-seq protocols (10x Chromium, Drop-seq)
   capture only the 3' end of transcripts.  If your query targets a region far
   from the 3' UTR, you may observe lower signal or no hits even when the gene
   is expressed.  Consider designing probes against the 3' end of the transcript
   when working with 3'-biased datasets.

Stranded Search
---------------

By default Malva searches the **forward strand only**, which is correct for
most RNA-seq queries where the query sequence is in the same orientation as the
transcript.

Pass ``stranded=False`` when you want to include reverse-complement k-mers:

.. code-block:: python

   # Both strands — useful when orientation is ambiguous
   results = client.search_sequences(probe, stranded=False)

   # Forward only (default)
   results = client.search_sequences(probe, stranded=True)
   results = client.search_sequences(probe)  # same as above

Default Behaviour
-----------------

When ``min_kmer_presence``, ``max_kmer_presence``, or ``stranded`` are omitted,
the server applies the defaults shown above (``0``, ``100000``, forward only).
These work well for the majority of expression queries.  Adjust them when you
need to reduce noise from rare k-mers (``min_kmer_presence``) or filter out
repetitive/ubiquitous k-mers (``max_kmer_presence``).
