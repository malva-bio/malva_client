class MalvaAPIError(Exception):
    """Base exception for Malva API errors"""
    pass

class AuthenticationError(MalvaAPIError):
    """Authentication related errors"""
    pass

class SearchError(MalvaAPIError):
    """Search related errors"""
    pass

class QuotaExceededError(MalvaAPIError):
    """Quota exceeded errors"""
    def __init__(self, message, quota_info=None):
        super().__init__(message)
        self.quota_info = quota_info