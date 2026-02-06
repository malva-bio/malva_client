Query Parameters
================

How Search Works
----------------

Malva indexes every sample with fixed-length **k-mers** (k = 24 nucleotides).
When you submit a query sequence, the engine:

1. Extracts all k-mers from the query.
2. Slides a window of **w** consecutive k-mers across the query.
3. For each window, counts how many k-mers are present in a cell's index.
4. Reports a hit when the fraction of matching k-mers meets or exceeds the
   **threshold** (tau).

The two user-facing parameters — ``window_size`` and ``threshold`` — let you
control this sensitivity/specificity trade-off.

Parameters
----------

.. grid:: 3
   :gutter: 3

   .. grid-item-card::
      :class-card: sd-border-primary sd-rounded-3
      :class-header: sd-bg-primary sd-text-white

      Window Size (w)
      ^^^
      Number of consecutive k-mers per sliding window.  A window of *w* k-mers
      covers *w + 23* nucleotides.

      * **Larger** (64--120) — robust to mismatches, best for expression
      * **Smaller** (24) — contiguous region, best for junctions/SNVs

   .. grid-item-card::
      :class-card: sd-border-success sd-rounded-3
      :class-header: sd-bg-success sd-text-white

      Threshold (tau)
      ^^^
      Minimum fraction of k-mers in a window that must match (0.0--1.0).

      * **Lower** (0.5--0.65) — tolerates mismatches, expression queries
      * **Higher** (0.9--1.0) — near-exact, junctions/SNV screening

   .. grid-item-card::
      :class-card: sd-border-info sd-rounded-3
      :class-header: sd-bg-info sd-text-white

      Strandedness
      ^^^
      By default both strands are searched.  Set ``stranded=True`` to restrict
      to the forward strand only.

      Useful for antisense oligos or distinguishing sense vs. antisense
      transcription.

Recommended Parameters
----------------------

The table below summarises recommended settings for common use cases
(based on Supplementary Table 2 of the Malva manuscript).

.. list-table::
   :header-rows: 1
   :widths: 30 20 20 30
   :class: sd-table-hover

   * - Use Case
     - Window Size (w)
     - Threshold (tau)
     - Notes
   * - Transcript expression
     - 64--120
     - 0.50--0.65
     - Default-like; tolerant of mismatches
   * - Splice junction detection
     - 24
     - 1.0
     - Exact match over junction sequence
   * - Circular RNA (back-splice)
     - 24
     - 1.0
     - Query should span the BSJ
   * - SNV / short variant detection
     - 24
     - 1.0
     - Probe covers variant +/- 12 nt
   * - Isoform-level analysis
     - 64--96
     - 0.50--0.65
     - Use isoform-specific exon sequences

Probe Design Guidelines
-----------------------

.. tab-set::

   .. tab-item:: Transcript Expression

      Use the full coding sequence or a representative exon.  Larger windows
      average over local noise.

      .. code-block:: python

         results = client.search("BRCA1", window_size=96, threshold=0.55)

   .. tab-item:: Splice Junction

      Design a probe of ~48 nt centred on the junction (24 nt from each exon).
      Use ``window_size=24`` and ``threshold=1.0`` so only reads spanning the
      exact junction are counted.

      .. code-block:: python

         junction = "ACGTACGT" * 6  # 48 nt spanning junction
         results = client.search_sequence(
             junction, window_size=24, threshold=1.0
         )

   .. tab-item:: Circular RNA (BSJ)

      Same principle: design a probe that spans the back-splice junction and
      use exact matching.

      .. code-block:: python

         bsj_probe = "CTAG" * 12  # 48 nt across BSJ
         results = client.search_sequence(
             bsj_probe, window_size=24, threshold=1.0
         )

   .. tab-item:: SNV Detection

      Centre a 48 nt probe on the variant position.  The variant must fall
      within the window so that exact matching distinguishes reference from
      alternative alleles.

      .. code-block:: python

         snv_probe = ref_context[:24] + alt_base + ref_context[25:48]
         results = client.search_sequence(
             snv_probe, window_size=24, threshold=1.0
         )

.. warning::

   **3'-biased protocols** — Many scRNA-seq protocols (10x Chromium, Drop-seq)
   capture only the 3' end of transcripts.  If your query targets a region far
   from the 3' UTR, you may observe lower signal or no hits even when the gene
   is expressed.  Consider designing probes against the 3' end of the transcript
   when working with 3'-biased datasets.

Stranded Search
---------------

Pass ``stranded=True`` when the orientation of your query matters:

.. code-block:: python

   # Antisense probe — only count forward-strand matches
   results = client.search_sequence(antisense_oligo, stranded=True)

Omitting the parameter (or setting it to ``False``) searches both strands,
which is the correct default for most expression queries.

Default Behaviour
-----------------

When ``window_size``, ``threshold``, or ``stranded`` are omitted, the server
applies its own defaults (currently tuned for general transcript-level
expression queries).  You only need to set these parameters when your use case
requires tighter or looser matching than the default.
