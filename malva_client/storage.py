"""
Storage module for Malva client
Manages search history, job status tracking, and local data persistence
"""

import json
import sqlite3
import time
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Any

from rich.console import Console

console = Console()


class MalvaStorage:
    """Manages local storage for Malva client data"""
    
    def __init__(self, storage_dir: Optional[Path] = None):
        """
        Initialize storage manager
        
        Args:
            storage_dir: Directory for storage files. Defaults to ~/.malva/
        """
        if storage_dir is None:
            storage_dir = Path.home() / '.malva'
        
        self.storage_dir = Path(storage_dir)
        self.storage_dir.mkdir(parents=True, exist_ok=True)
        
        self.db_path = self.storage_dir / 'malva.db'
        self._init_database()
    
    def _init_database(self):
        """Initialize SQLite database with required tables"""
        with sqlite3.connect(self.db_path) as conn:
            conn.execute('''
                CREATE TABLE IF NOT EXISTS search_history (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    timestamp TEXT NOT NULL,
                    query TEXT NOT NULL,
                    job_id TEXT,
                    status TEXT NOT NULL DEFAULT 'submitted',
                    server_url TEXT,
                    results_count INTEGER DEFAULT 0,
                    error_message TEXT,
                    metadata TEXT
                )
            ''')
            
            conn.execute('''
                CREATE TABLE IF NOT EXISTS job_status (
                    job_id TEXT PRIMARY KEY,
                    status TEXT NOT NULL,
                    timestamp TEXT NOT NULL,
                    query TEXT,
                    server_url TEXT,
                    results_data TEXT,
                    error_message TEXT,
                    metadata TEXT
                )
            ''')
            
            # Create indexes for performance
            conn.execute('CREATE INDEX IF NOT EXISTS idx_search_timestamp ON search_history(timestamp)')
            conn.execute('CREATE INDEX IF NOT EXISTS idx_job_timestamp ON job_status(timestamp)')
            
            conn.commit()
    
    def save_search(self, query: str, job_id: Optional[str] = None, 
                   server_url: Optional[str] = None, status: str = 'submitted',
                   results_count: int = 0, error_message: Optional[str] = None,
                   metadata: Optional[Dict] = None) -> int:
        """
        Save a search to history
        
        Args:
            query: Search query
            job_id: Job ID if available
            server_url: Server URL
            status: Search status
            results_count: Number of results found
            error_message: Error message if failed
            metadata: Additional metadata
            
        Returns:
            Search history ID
        """
        timestamp = datetime.now().isoformat()
        metadata_json = json.dumps(metadata) if metadata else None
        
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.execute('''
                INSERT INTO search_history 
                (timestamp, query, job_id, status, server_url, results_count, error_message, metadata)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            ''', (timestamp, query, job_id, status, server_url, results_count, error_message, metadata_json))
            
            conn.commit()
            return cursor.lastrowid
    
    def update_search_status(self, search_id: int, status: str, 
                           results_count: int = 0, error_message: Optional[str] = None):
        """
        Update search status
        
        Args:
            search_id: Search history ID
            status: New status
            results_count: Number of results
            error_message: Error message if any
        """
        with sqlite3.connect(self.db_path) as conn:
            conn.execute('''
                UPDATE search_history 
                SET status = ?, results_count = ?, error_message = ?
                WHERE id = ?
            ''', (status, results_count, error_message, search_id))
            conn.commit()
    
    def get_search_history(self, limit: int = 20, server_url: Optional[str] = None) -> List[Dict]:
        """
        Get search history
        
        Args:
            limit: Maximum number of results
            server_url: Filter by server URL
            
        Returns:
            List of search history entries
        """
        with sqlite3.connect(self.db_path) as conn:
            conn.row_factory = sqlite3.Row
            
            if server_url:
                cursor = conn.execute('''
                    SELECT * FROM search_history 
                    WHERE server_url = ?
                    ORDER BY timestamp DESC 
                    LIMIT ?
                ''', (server_url, limit))
            else:
                cursor = conn.execute('''
                    SELECT * FROM search_history 
                    ORDER BY timestamp DESC 
                    LIMIT ?
                ''', (limit,))
            
            return [dict(row) for row in cursor.fetchall()]
    
    def save_job_status(self, job_id: str, status: str, query: Optional[str] = None,
                       server_url: Optional[str] = None, results_data: Optional[Dict] = None,
                       error_message: Optional[str] = None, metadata: Optional[Dict] = None):
        """
        Save or update job status
        
        Args:
            job_id: Job ID
            status: Job status
            query: Original query
            server_url: Server URL
            results_data: Results data if completed
            error_message: Error message if failed
            metadata: Additional metadata
        """
        timestamp = datetime.now().isoformat()
        results_json = json.dumps(results_data) if results_data else None
        metadata_json = json.dumps(metadata) if metadata else None
        
        with sqlite3.connect(self.db_path) as conn:
            conn.execute('''
                INSERT OR REPLACE INTO job_status 
                (job_id, status, timestamp, query, server_url, results_data, error_message, metadata)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            ''', (job_id, status, timestamp, query, server_url, results_json, error_message, metadata_json))
            conn.commit()
    
    def get_job_status(self, job_id: str) -> Optional[Dict]:
        """
        Get job status
        
        Args:
            job_id: Job ID
            
        Returns:
            Job status data or None if not found
        """
        with sqlite3.connect(self.db_path) as conn:
            conn.row_factory = sqlite3.Row
            cursor = conn.execute('''
                SELECT * FROM job_status WHERE job_id = ?
            ''', (job_id,))
            
            row = cursor.fetchone()
            if row:
                data = dict(row)
                # Parse JSON fields
                if data['results_data']:
                    data['results_data'] = json.loads(data['results_data'])
                if data['metadata']:
                    data['metadata'] = json.loads(data['metadata'])
                return data
            return None
    
    def get_pending_jobs(self, server_url: Optional[str] = None) -> List[Dict]:
        """
        Get pending jobs
        
        Args:
            server_url: Filter by server URL
            
        Returns:
            List of pending jobs
        """
        with sqlite3.connect(self.db_path) as conn:
            conn.row_factory = sqlite3.Row
            
            if server_url:
                cursor = conn.execute('''
                    SELECT * FROM job_status 
                    WHERE status IN ('pending', 'running') AND server_url = ?
                    ORDER BY timestamp DESC
                ''', (server_url,))
            else:
                cursor = conn.execute('''
                    SELECT * FROM job_status 
                    WHERE status IN ('pending', 'running')
                    ORDER BY timestamp DESC
                ''')
            
            return [dict(row) for row in cursor.fetchall()]
    
    def cleanup_old_data(self, days: int = 30):
        """
        Clean up old data
        
        Args:
            days: Number of days to keep data
        """
        cutoff_date = (datetime.now() - timedelta(days=days)).isoformat()
        
        with sqlite3.connect(self.db_path) as conn:
            # Clean up old search history
            cursor = conn.execute('''
                DELETE FROM search_history WHERE timestamp < ?
            ''', (cutoff_date,))
            history_deleted = cursor.rowcount
            
            # Clean up old completed jobs
            cursor = conn.execute('''
                DELETE FROM job_status 
                WHERE timestamp < ? AND status IN ('completed', 'error', 'failed')
            ''', (cutoff_date,))
            jobs_deleted = cursor.rowcount
            
            conn.commit()
            
            if history_deleted > 0 or jobs_deleted > 0:
                console.print(f"Cleaned up {history_deleted} old searches and {jobs_deleted} old jobs")
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get storage statistics
        
        Returns:
            Dictionary with storage statistics
        """
        with sqlite3.connect(self.db_path) as conn:
            # Count searches
            cursor = conn.execute('SELECT COUNT(*) FROM search_history')
            total_searches = cursor.fetchone()[0]
            
            # Count jobs by status
            cursor = conn.execute('''
                SELECT status, COUNT(*) FROM job_status GROUP BY status
            ''')
            job_counts = dict(cursor.fetchall())
            
            # Get recent activity (last 7 days)
            week_ago = (datetime.now() - timedelta(days=7)).isoformat()
            cursor = conn.execute('''
                SELECT COUNT(*) FROM search_history WHERE timestamp > ?
            ''', (week_ago,))
            recent_searches = cursor.fetchone()[0]
            
            return {
                'total_searches': total_searches,
                'recent_searches': recent_searches,
                'job_counts': job_counts,
                'database_size': self.db_path.stat().st_size if self.db_path.exists() else 0
            }
    
    def export_history(self, output_path: Path, format: str = 'json'):
        """
        Export search history
        
        Args:
            output_path: Output file path
            format: Export format ('json' or 'csv')
        """
        history = self.get_search_history(limit=1000)  # Get more for export
        
        if format == 'json':
            with open(output_path, 'w') as f:
                json.dump(history, f, indent=2)
        elif format == 'csv':
            try:
                import pandas as pd
                df = pd.DataFrame(history)
                df.to_csv(output_path, index=False)
            except ImportError:
                raise ImportError("pandas required for CSV export. Install with: pip install pandas")
        else:
            raise ValueError(f"Unsupported export format: {format}")
        
        console.print(f"Exported {len(history)} entries to {output_path}")


# Global storage instance
_storage = None

def get_storage(storage_dir: Optional[Path] = None) -> MalvaStorage:
    """
    Get global storage instance
    
    Args:
        storage_dir: Storage directory
        
    Returns:
        MalvaStorage instance
    """
    global _storage
    if _storage is None:
        _storage = MalvaStorage(storage_dir)
    return _storage


# Convenience functions
def save_search(query: str, **kwargs) -> int:
    """Save search to history"""
    return get_storage().save_search(query, **kwargs)


def get_search_history(limit: int = 20, **kwargs) -> List[Dict]:
    """Get search history"""
    return get_storage().get_search_history(limit, **kwargs)


def save_job_status(job_id: str, status: str, **kwargs):
    """Save job status"""
    get_storage().save_job_status(job_id, status, **kwargs)


def get_job_status(job_id: str) -> Optional[Dict]:
    """Get job status"""
    return get_storage().get_job_status(job_id)


def cleanup_old_data(days: int = 30):
    """Clean up old data"""
    get_storage().cleanup_old_data(days)