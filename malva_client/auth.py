"""
Authentication module for Malva client
Handles interactive login and secure token storage
"""

import getpass
import sys
from typing import Optional
from urllib.parse import urljoin

import requests
from rich.console import Console
from rich.prompt import Prompt

from malva_client.exceptions import AuthenticationError, MalvaAPIError

console = Console()

try:
    import keyring
    HAS_KEYRING = True
except ImportError:
    HAS_KEYRING = False
    console.print("[yellow]Warning: keyring not available. Tokens will not be stored securely.[/yellow]")


class AuthManager:
    """Manages authentication and token storage"""
    
    SERVICE_NAME = "malva_client"
    USERNAME = "api-token"
    
    def __init__(self, server_url: str):
        self.server_url = server_url.rstrip('/')
    
    def login_interactive(self) -> str:
        """
        Perform interactive login and return API token
        
        Returns:
            API token string
        """
        console.print(f"[bold]Malva Client Login[/bold]")
        console.print(f"Server: {self.server_url}")
        console.print()
        
        # Check if we have a stored token first
        stored_token = self.get_stored_token()
        if stored_token:
            use_stored = Prompt.ask(
                "Found stored API token. Use it?", 
                choices=["y", "n"], 
                default="y"
            )
            if use_stored.lower() == "y":
                if self._validate_token(stored_token):
                    console.print("[green]Using stored token - login successful![/green]")
                    return stored_token
                else:
                    console.print("[yellow]Stored token is invalid, continuing with new login...[/yellow]")
                    self.logout()  # Remove invalid token
        
        console.print("Please get your API token from:")
        console.print(f"1. Go to {self.server_url}")
        console.print("2. Click the menu button (⋯) in the upper left")
        console.print("3. Select 'Profile' from the dropdown")
        console.print("4. Click 'Generate API Token'")
        console.print("5. Copy the token and paste it below")
        console.print()
        
        # Get token from user
        token = getpass.getpass("Enter your API token: ").strip()
        
        if not token:
            raise AuthenticationError("No token provided")
        
        # Validate token
        if not self._validate_token(token):
            raise AuthenticationError("Invalid token provided")
        
        # Store token securely
        if self._store_token(token):
            console.print("[green]Token stored securely for future use[/green]")
        
        console.print("[green]Login successful![/green]")
        return token
    
    def logout(self) -> bool:
        """
        Remove stored token
        
        Returns:
            True if token was removed successfully
        """
        if not HAS_KEYRING:
            console.print("[yellow]No secure storage available - nothing to logout from[/yellow]")
            return True
        
        try:
            keyring.delete_password(self.SERVICE_NAME, self.USERNAME)
            console.print("[green]Logged out successfully[/green]")
            return True
        except keyring.errors.PasswordDeleteError:
            console.print("[yellow]No stored token found[/yellow]")
            return True
        except Exception as e:
            console.print(f"[red]Error during logout: {e}[/red]")
            return False
    
    def get_stored_token(self) -> Optional[str]:
        """
        Get stored API token
        
        Returns:
            Stored token or None if not found
        """
        if not HAS_KEYRING:
            return None
        
        try:
            token = keyring.get_password(self.SERVICE_NAME, self.USERNAME)
            return token
        except Exception:
            return None
    
    def _store_token(self, token: str) -> bool:
        """
        Store token securely
        
        Args:
            token: API token to store
            
        Returns:
            True if stored successfully
        """
        if not HAS_KEYRING:
            console.print("[yellow]Keyring not available - token not stored[/yellow]")
            return False
        
        try:
            keyring.set_password(self.SERVICE_NAME, self.USERNAME, token)
            return True
        except Exception as e:
            console.print(f"[yellow]Could not store token securely: {e}[/yellow]")
            return False
    
    def _validate_token(self, token: str) -> bool:
        """
        Validate token with the server
        
        Args:
            token: Token to validate
            
        Returns:
            True if token is valid
        """
        try:
            session = requests.Session()
            session.headers.update({'Authorization': f'Bearer {token}'})
            
            # Test with quota endpoint
            response = session.get(
                urljoin(self.server_url, '/api/quota-status'),
                timeout=10
            )
            
            if response.status_code == 200:
                return True
            elif response.status_code in [401, 403]:
                return False
            else:
                # If quota endpoint doesn't exist, try health endpoint
                response = session.get(
                    urljoin(self.server_url, '/health'),
                    timeout=10
                )
                return response.status_code == 200
                
        except Exception:
            return False


def login(server_url: str) -> str:
    """
    Convenience function for interactive login
    
    Args:
        server_url: Malva server URL
        
    Returns:
        API token
    """
    auth_manager = AuthManager(server_url)
    return auth_manager.login_interactive()


def logout(server_url: str) -> bool:
    """
    Convenience function for logout
    
    Args:
        server_url: Malva server URL
        
    Returns:
        True if successful
    """
    auth_manager = AuthManager(server_url)
    return auth_manager.logout()


def get_stored_token(server_url: str) -> Optional[str]:
    """
    Convenience function to get stored token
    
    Args:
        server_url: Malva server URL
        
    Returns:
        Stored token or None
    """
    auth_manager = AuthManager(server_url)
    return auth_manager.get_stored_token()