__version__ = "0.1.0"
__author__ = "Malva Team"
__email__ = "hello@malva.bio"

from .client import MalvaClient
from .exceptions import MalvaAPIError, AuthenticationError, SearchError, QuotaExceededError
from .config import Config
from .models import SearchResult, CoverageResult

__all__ = [
    'MalvaClient',
    'MalvaAPIError',
    'AuthenticationError', 
    'SearchError',
    'QuotaExceededError',
    'SearchResult',
    'CoverageResult',
    'Config'
]
