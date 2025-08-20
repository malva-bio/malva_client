from typing import Dict, Any, Optional


class MalvaAPIError(Exception):
    """Base exception for Malva API errors"""
    
    def __init__(self, message: str, response_data: Optional[Dict[str, Any]] = None):
        self.message = message
        self.response_data = response_data or {}
        super().__init__(message)


class AuthenticationError(MalvaAPIError):
    """Authentication failed"""
    pass


class SearchError(MalvaAPIError):
    """Search operation failed"""
    
    def __init__(self, message: str, job_id: Optional[str] = None, 
                 response_data: Optional[Dict[str, Any]] = None):
        self.job_id = job_id
        super().__init__(message, response_data)


class QuotaExceededError(MalvaAPIError):
    """API quota exceeded"""
    
    def __init__(self, message: str, quota_info: Optional[Dict[str, Any]] = None):
        self.quota_info = quota_info or {}
        super().__init__(message, quota_info)


class ValidationError(MalvaAPIError):
    """Input validation failed"""
    pass


class ConfigurationError(MalvaAPIError):
    """Configuration error"""
    pass


class StorageError(MalvaAPIError):
    """Local storage error"""
    pass


class NetworkError(MalvaAPIError):
    """Network connectivity error"""
    pass


class TimeoutError(MalvaAPIError):
    """Operation timed out"""
    pass