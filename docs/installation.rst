Installation
============

Requires **Python 3.8** or newer.

.. tab-set::

   .. tab-item:: PyPI (recommended)

      Install the latest stable release from PyPI:

      .. code-block:: bash

         pip install malva_client

      This pulls in all required dependencies automatically.

   .. tab-item:: From Source

      Clone the repository and install in editable mode for development:

      .. code-block:: bash

         git clone https://github.com/malva-bio/malva_client
         cd malva_client
         pip install -e .

      Editable mode (``-e``) means changes you make to the source code take
      effect immediately without reinstalling.

   .. tab-item:: Optional Extras

      For full functionality you may want these additional packages:

      .. code-block:: bash

         # Excel export support
         pip install openpyxl

         # Download samples as AnnData objects
         pip install anndata h5py

         # Single-cell analysis workflows
         pip install scanpy

Verify Installation
-------------------

Run this one-liner to confirm everything is working:

.. code-block:: bash

   python -c "import malva_client; print(malva_client.__version__)"

You should see the installed version number (e.g. ``0.2.0``).

Next Steps
----------

.. grid:: 1
   :gutter: 2

   .. grid-item-card::
      :link: quickstart
      :link-type: doc
      :class-card: sd-border-success
      :class-header: sd-bg-success sd-text-white

      Quick Start
      ^^^
      Configure authentication and run your first search in minutes.
