"""
Malva Client - Python client for the Malva genomic search platform
"""

try:
    from importlib.metadata import version as _get_version
    __version__ = _get_version("malva_client")
except Exception:
    __version__ = "0.2.0"
__author__ = "Malva Team"
__email__ = "hello@malva.bio"

# Define what should be available when importing, but don't actually import yet
__all__ = [
    'MalvaClient',
    'MalvaAPIError',
    'AuthenticationError',
    'SearchError',
    'QuotaExceededError',
    'ValidationError',
    'ConfigurationError',
    'SearchResult',
    'CoverageResult',
    'SingleCellResult',
    'CoexpressionResult',
    'UMAPCoordinates',
    'Config',
    'search_gene',
    'search_sequence',
    'login',
    'logout',
    'get_stored_token'
]

def __getattr__(name):
    """
    Lazy import mechanism - only import modules when they're actually used
    """
    if name == 'MalvaClient':
        from .client import MalvaClient
        return MalvaClient

    elif name in ('MalvaAPIError', 'AuthenticationError', 'SearchError',
                  'QuotaExceededError', 'ValidationError', 'ConfigurationError'):
        from .exceptions import (
            MalvaAPIError, AuthenticationError, SearchError,
            QuotaExceededError, ValidationError, ConfigurationError
        )
        return locals()[name]

    elif name in ('SearchResult', 'CoverageResult', 'SingleCellResult',
                  'CoexpressionResult', 'UMAPCoordinates'):
        from .models import (SearchResult, CoverageResult, SingleCellResult,
                             CoexpressionResult, UMAPCoordinates)
        return locals()[name]

    elif name == 'Config':
        from .config import Config
        return Config
    
    elif name in ('search_gene', 'search_sequence'):
        from .client import search_gene, search_sequence
        return locals()[name]

    elif name in ('login', 'logout', 'get_stored_token'):
        from .auth import login, logout, get_stored_token
        return locals()[name]

    else:
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'")

# For backwards compatibility and to make imports work in older Python versions
def __dir__():
    """Return the list of available attributes"""
    return __all__
