import json
import os
from pathlib import Path
from typing import Optional

class Config:
    """Configuration management for Malva client"""
    
    def __init__(self, config_path: Optional[str] = None):
        self.config_path = Path(config_path) if config_path else Path.home() / '.malva' / 'config.json'
        self.server_url = None
        self.api_token = None
        self.verify_ssl = True
        
    @classmethod
    def load(cls, config_path: Optional[str] = None):
        """Load configuration from file or environment"""
        config = cls(config_path)
        
        # Load from file if exists
        if config.config_path.exists():
            try:
                with open(config.config_path) as f:
                    data = json.load(f)
                config.server_url = data.get('server_url')
                config.api_token = data.get('api_token')
                config.verify_ssl = data.get('verify_ssl', True)
            except Exception:
                pass
        
        # Override with environment variables
        config.server_url = os.environ.get('MALVA_API_URL', config.server_url)
        config.api_token = os.environ.get('MALVA_API_TOKEN', config.api_token)
        if 'MALVA_VERIFY_SSL' in os.environ:
            config.verify_ssl = os.environ.get('MALVA_VERIFY_SSL', 'true').lower() == 'true'
            
        return config
    
    def save(self):
        """Save configuration to file"""
        self.config_path.parent.mkdir(parents=True, exist_ok=True)
        data = {
            'server_url': self.server_url,
            'api_token': self.api_token,
            'verify_ssl': self.verify_ssl
        }
        with open(self.config_path, 'w') as f:
            json.dump(data, f, indent=2)