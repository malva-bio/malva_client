import os
import sys
sys.path.insert(0, os.path.abspath('..'))

project = 'Malva Client'
copyright = '2025-2026, Malva Team'
author = 'Malva Team'
try:
    from importlib.metadata import version as _get_version
    release = _get_version("malva_client")
except Exception:
    release = '0.2.0'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx_autodoc_typehints',
    'myst_parser',
    'sphinx_click',
    'sphinx_copybutton',
    'sphinx_design',
    'nbsphinx',
]

nbsphinx_execute = 'never'

# FURO THEME
html_theme = 'furo'
html_title = "Malva Client"
html_logo = "_static/malva_logo.svg"

html_theme_options = {
    "sidebar_hide_name": False,
    "navigation_with_keys": True,
    "source_repository": "https://github.com/malva-bio/malva_client/",
    "source_branch": "main",
    "source_directory": "docs/",
}

# Static files
html_static_path = ['_static']

# Settings
suppress_warnings = ['autodoc.duplicate_object']
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'undoc-members': True,
    'exclude-members': '__weakref__'
}

napoleon_google_docstring = True
napoleon_numpy_docstring = True
source_suffix = {'.rst': None, '.md': 'myst_parser'}
master_doc = 'index'
html_search_language = 'en'

autodoc_mock_imports = ['h5py', 'anndata', 'dnaio']

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
}
