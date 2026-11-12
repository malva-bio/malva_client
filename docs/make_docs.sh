#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
rm -rf _build
sphinx-build -W --keep-going -b html . _build/html
echo "Documentation built in _build/html/"
echo "Open _build/html/index.html in your browser"
