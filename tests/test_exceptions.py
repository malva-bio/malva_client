"""Tests for malva_client.exceptions."""

from malva_client.exceptions import (
    MalvaAPIError,
    AuthenticationError,
    SearchError,
    QuotaExceededError,
    ValidationError,
    ConfigurationError,
    StorageError,
    NetworkError,
    TimeoutError,
)


class TestMalvaAPIError:
    def test_message_and_response_data(self):
        err = MalvaAPIError("something broke", {"code": 500})
        assert err.message == "something broke"
        assert err.response_data == {"code": 500}
        assert str(err) == "something broke"

    def test_defaults_response_data_to_empty_dict(self):
        err = MalvaAPIError("oops")
        assert err.response_data == {}


class TestAuthenticationError:
    def test_inherits_from_base(self):
        err = AuthenticationError("bad token")
        assert isinstance(err, MalvaAPIError)
        assert err.message == "bad token"


class TestSearchError:
    def test_has_job_id(self):
        err = SearchError("failed", job_id="job-42")
        assert err.job_id == "job-42"
        assert isinstance(err, MalvaAPIError)

    def test_job_id_defaults_to_none(self):
        err = SearchError("failed")
        assert err.job_id is None


class TestQuotaExceededError:
    def test_has_quota_info(self):
        info = {"account_type": "free", "quota_exceeded": True}
        err = QuotaExceededError("limit reached", quota_info=info)
        assert err.quota_info == info
        assert isinstance(err, MalvaAPIError)

    def test_defaults_quota_info_to_empty_dict(self):
        err = QuotaExceededError("limit")
        assert err.quota_info == {}


class TestAllExceptionsInheritFromBase:
    def test_all_inherit(self):
        classes = [
            AuthenticationError,
            SearchError,
            QuotaExceededError,
            ValidationError,
            ConfigurationError,
            StorageError,
            NetworkError,
            TimeoutError,
        ]
        for cls in classes:
            assert issubclass(cls, MalvaAPIError), f"{cls.__name__} should inherit MalvaAPIError"


class TestExceptionStr:
    def test_str_returns_message(self):
        for cls in [MalvaAPIError, AuthenticationError, ValidationError,
                    ConfigurationError, StorageError, NetworkError, TimeoutError]:
            err = cls("test message")
            assert str(err) == "test message"
