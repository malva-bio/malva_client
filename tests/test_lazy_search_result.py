"""Lazy SearchResult accessors must trigger the deferred results fetch."""

from typing import TYPE_CHECKING, cast

import pandas as pd
import pytest

from malva_client.models import SearchResult

if TYPE_CHECKING:
    from malva_client.client import MalvaClient


class _StubResult(SearchResult):
    """SearchResult whose fetch is replaced by a counter."""

    calls = 0

    def _ensure_loaded(self):
        if not getattr(self, "_needs_lazy_fetch", False):
            return
        type(self).calls += 1
        self._needs_lazy_fetch = False
        self._df = pd.DataFrame({"a": [1, 2, 3]})


@pytest.fixture(autouse=True)
def _reset_counter():
    _StubResult.calls = 0


def _new_result() -> _StubResult:
    return _StubResult({"job_id": "job-1"}, cast("MalvaClient", None))


@pytest.mark.parametrize(
    "accessor",
    [
        lambda r: r.shape,
        lambda r: len(r),
        lambda r: r["a"],
        lambda r: r.head(),
        lambda r: r.columns,
    ],
)
def test_lazy_access_triggers_single_fetch(accessor):
    result = _new_result()
    assert _StubResult.calls == 0

    accessor(result)

    assert _StubResult.calls == 1


def test_lazy_fetch_happens_only_once():
    result = _new_result()

    assert result.shape == (3, 1)
    assert len(result) == 3
    assert list(result["a"]) == [1, 2, 3]

    assert _StubResult.calls == 1


def test_no_fetch_when_results_already_present():
    loaded = SearchResult.__new__(SearchResult)
    loaded._needs_lazy_fetch = False
    loaded._df = pd.DataFrame({"a": [1]})
    loaded._enriched = True

    assert loaded.shape == (1, 1)
