import json
import os
from pathlib import Path
from typing import Optional, Dict, List, Any

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

    def set_value(self, key: str, value: str):
        """
        Set configuration value
        
        Args:
            key: Configuration key
            value: Configuration value
        """
        valid_keys = {
            'server_url': 'server_url',
            'server': 'server_url',
            'url': 'server_url',
            'api_token': 'api_token',
            'token': 'api_token',
            'verify_ssl': 'verify_ssl',
            'ssl': 'verify_ssl'
        }
        
        if key not in valid_keys:
            raise ValueError(f"Invalid configuration key: {key}. Valid keys: {list(valid_keys.keys())}")
        
        attr_name = valid_keys[key]
        
        # Convert value types
        if attr_name == 'verify_ssl':
            value = value.lower() in ('true', '1', 'yes', 'on')
        
        setattr(self, attr_name, value)

    def get_value(self, key: str) -> Optional[str]:
        """
        Get configuration value
        
        Args:
            key: Configuration key
            
        Returns:
            Configuration value or None
        """
        valid_keys = {
            'server_url': 'server_url',
            'server': 'server_url',
            'url': 'server_url',
            'api_token': 'api_token',
            'token': 'api_token',
            'verify_ssl': 'verify_ssl',
            'ssl': 'verify_ssl'
        }
        
        if key not in valid_keys:
            return None
        
        attr_name = valid_keys[key]
        return getattr(self, attr_name)

    def show(self) -> Dict[str, Any]:
        """
        Show current configuration
        
        Returns:
            Dictionary with current configuration
        """
        return {
            'server_url': self.server_url,
            'api_token': '***' if self.api_token else None,
            'verify_ssl': self.verify_ssl,
            'config_path': str(self.config_path)
        }

    def validate(self) -> List[str]:
        """
        Validate current configuration
        
        Returns:
            List of validation errors (empty if valid)
        """
        errors = []
        
        if not self.server_url:
            errors.append("server_url is required")
        elif not self.server_url.startswith(('http://', 'https://')):
            errors.append("server_url must start with http:// or https://")
        
        if not self.api_token:
            errors.append("api_token is required")
        
        return errors

    def reset_to_defaults(self):
        """Reset configuration to defaults"""
        self.server_url = None
        self.api_token = None
        self.verify_ssl = True

    def delete_config_file(self) -> bool:
        """
        Delete configuration file
        
        Returns:
            True if file was deleted
        """
        try:
            if self.config_path.exists():
                self.config_path.unlink()
                return True
            return False
        except Exception:
            return False