Installation
============

Requirements
------------

* Python 3.8+
* pandas, requests, click, rich

Install from PyPI
-----------------

.. code-block:: bash

   pip install malva_client

Install from source
-------------------

.. code-block:: bash

   git clone https://github.com/malva-bio/malva_client
   cd malva_client
   pip install -e .

Dependencies
------------

Required dependencies:

- requests>=2.25.0
- pandas>=1.3.0
- click>=8.0.0
- rich>=10.0.0
- keyring>=23.0.0

Optional dependencies for full functionality:

- openpyxl>=3.0.0 (Excel export)
- anndata>=0.8.0 (sample downloads)
- h5py>=3.0.0 (H5AD files)