#!/bin/bash

# Clean previous build
rm -rf _build/
rm -rf api/

# Generate API documentation
sphinx-apidoc -o api/ ../malva_client --force --separate --no-toc

# Ensure api/modules.rst exists with correct content
cat > api/modules.rst << 'EOF'
API Reference
=============

.. toctree::
   :maxdepth: 4

   malva_client.client
   malva_client.models
   malva_client.exceptions
   malva_client.auth
   malva_client.config
   malva_client.storage
   malva_client.tools
   malva_client.cli
EOF

# Build HTML documentation
sphinx-build -b html . _build/html

echo "Documentation built in _build/html/"
echo "Open _build/html/index.html in your browser"