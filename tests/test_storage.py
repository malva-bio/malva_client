"""Tests for malva_client.storage."""

import json
import time
from datetime import datetime, timedelta
from pathlib import Path

import pytest

from malva_client.storage import MalvaStorage


@pytest.fixture
def storage(tmp_path):
    """Create a MalvaStorage using a temporary directory."""
    return MalvaStorage(storage_dir=tmp_path)


class TestInit:
    def test_creates_database(self, tmp_path):
        s = MalvaStorage(storage_dir=tmp_path)
        assert (tmp_path / "malva.db").exists()


class TestSaveSearch:
    def test_inserts_and_returns_id(self, storage):
        search_id = storage.save_search(query="BRCA1", job_id="j1", server_url="https://s.io")
        assert isinstance(search_id, int)
        assert search_id >= 1

    def test_multiple_inserts_increment(self, storage):
        id1 = storage.save_search(query="TP53")
        id2 = storage.save_search(query="EGFR")
        assert id2 > id1


class TestGetSearchHistory:
    def test_returns_ordered_by_timestamp_desc(self, storage):
        storage.save_search(query="first")
        storage.save_search(query="second")

        history = storage.get_search_history(limit=10)
        assert len(history) == 2
        assert history[0]["query"] == "second"
        assert history[1]["query"] == "first"

    def test_filter_by_server_url(self, storage):
        storage.save_search(query="a", server_url="https://s1.io")
        storage.save_search(query="b", server_url="https://s2.io")

        history = storage.get_search_history(limit=10, server_url="https://s1.io")
        assert len(history) == 1
        assert history[0]["query"] == "a"

    def test_limit_works(self, storage):
        for i in range(5):
            storage.save_search(query=f"q{i}")
        history = storage.get_search_history(limit=2)
        assert len(history) == 2


class TestJobStatus:
    def test_save_and_get(self, storage):
        storage.save_job_status(
            job_id="j1", status="pending",
            query="BRCA1", server_url="https://s.io",
        )
        result = storage.get_job_status("j1")
        assert result is not None
        assert result["status"] == "pending"
        assert result["query"] == "BRCA1"

    def test_upsert_replaces(self, storage):
        storage.save_job_status(job_id="j1", status="pending")
        storage.save_job_status(job_id="j1", status="completed",
                                results_data={"foo": "bar"})
        result = storage.get_job_status("j1")
        assert result["status"] == "completed"
        assert result["results_data"] == {"foo": "bar"}

    def test_not_found_returns_none(self, storage):
        assert storage.get_job_status("nonexistent") is None

    def test_parses_json_fields(self, storage):
        storage.save_job_status(
            job_id="j1", status="completed",
            results_data={"key": [1, 2, 3]},
            metadata={"extra": True},
        )
        result = storage.get_job_status("j1")
        assert result["results_data"] == {"key": [1, 2, 3]}
        assert result["metadata"] == {"extra": True}


class TestGetPendingJobs:
    def test_filters_by_status(self, storage):
        storage.save_job_status(job_id="j1", status="pending")
        storage.save_job_status(job_id="j2", status="completed")
        storage.save_job_status(job_id="j3", status="running")

        pending = storage.get_pending_jobs()
        job_ids = [j["job_id"] for j in pending]
        assert "j1" in job_ids
        assert "j3" in job_ids
        assert "j2" not in job_ids

    def test_filters_by_server_url(self, storage):
        storage.save_job_status(job_id="j1", status="pending", server_url="https://a.io")
        storage.save_job_status(job_id="j2", status="pending", server_url="https://b.io")

        pending = storage.get_pending_jobs(server_url="https://a.io")
        assert len(pending) == 1
        assert pending[0]["job_id"] == "j1"


class TestCleanup:
    def test_deletes_old_entries(self, storage):
        import sqlite3

        # Insert a row with a very old timestamp directly
        old_ts = (datetime.now() - timedelta(days=60)).isoformat()
        with sqlite3.connect(storage.db_path) as conn:
            conn.execute(
                "INSERT INTO search_history (timestamp, query, status) VALUES (?, ?, ?)",
                (old_ts, "old_query", "completed"),
            )
            conn.execute(
                "INSERT INTO job_status (job_id, status, timestamp) VALUES (?, ?, ?)",
                ("old-job", "completed", old_ts),
            )
            conn.commit()

        # Also add a recent entry
        storage.save_search(query="recent")
        storage.save_job_status(job_id="new-job", status="pending")

        storage.cleanup_old_data(days=30)

        history = storage.get_search_history(limit=100)
        assert len(history) == 1
        assert history[0]["query"] == "recent"

        # old completed job should be gone, new pending job should remain
        assert storage.get_job_status("old-job") is None
        assert storage.get_job_status("new-job") is not None


class TestGetStats:
    def test_returns_counts(self, storage):
        storage.save_search(query="a")
        storage.save_search(query="b")
        storage.save_job_status(job_id="j1", status="completed")

        stats = storage.get_stats()
        assert stats["total_searches"] == 2
        assert stats["database_size"] > 0
        assert "completed" in stats["job_counts"]


class TestExportHistory:
    def test_export_json(self, storage, tmp_path):
        storage.save_search(query="export_me", job_id="e1")

        output = tmp_path / "export.json"
        storage.export_history(output, format="json")

        data = json.loads(output.read_text())
        assert len(data) == 1
        assert data[0]["query"] == "export_me"
