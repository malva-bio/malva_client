import sys
import time
import json
from pathlib import Path
from typing import Optional

import click
from rich.console import Console
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn
from rich.prompt import Prompt, Confirm
from rich.panel import Panel

console = Console()

@click.group()
@click.version_option()
@click.option('--config', type=click.Path(), help='Path to configuration file')
@click.option('--server', help='Malva server URL')
@click.option('--token', help='API token for authentication')
@click.option('--no-ssl-verify', is_flag=True, help='Disable SSL verification')
@click.pass_context
def main(ctx, config, server, token, no_ssl_verify):
    """Malva CLI - Search genomic datasets from the command line."""
    from malva_client.config import Config
    
    ctx.ensure_object(dict)
    
    config_obj = Config.load(config_path=config)

    if server:
        config_obj.server_url = server
    if token:
        config_obj.api_token = token
    if no_ssl_verify:
        config_obj.verify_ssl = False
    
    ctx.obj['config'] = config_obj

@main.command()
@click.argument('query')
@click.option('--output', '-o', type=click.Path(), help='Output file path')
@click.option('--format', 'output_format', type=click.Choice(['csv', 'json', 'excel', 'table']),
              default='table', help='Output format')
@click.option('--no-wait', is_flag=True, help='Submit search without waiting for completion')
@click.option('--max-wait', default=300, help='Maximum time to wait for completion (seconds)')
@click.option('--no-enrich', is_flag=True, help='Skip metadata enrichment')
@click.option('--window-size', '-w', type=int, default=None,
              help='Sliding window size (k-mers per window)')
@click.option('--threshold', '-t', type=float, default=None,
              help='Match threshold (0.0-1.0)')
@click.option('--stranded', is_flag=True, default=False,
              help='Restrict to single-strand search')
@click.pass_context
def search(ctx, query, output, output_format, no_wait, max_wait, no_enrich,
           window_size, threshold, stranded):
    """Search genomic datasets."""
    config = ctx.obj['config']
    from malva_client.client import MalvaClient
    from malva_client.exceptions import MalvaAPIError, AuthenticationError, QuotaExceededError
    from malva_client.storage import save_search
    
    try:
        client = MalvaClient(
            base_url=config.server_url,
            api_token=config.api_token,
            verify_ssl=config.verify_ssl
        )
        
        console.print(f"Searching for: [bold]{query}[/bold]")
        
        search_kwargs = {}
        if window_size is not None:
            search_kwargs['window_size'] = window_size
        if threshold is not None:
            search_kwargs['threshold'] = threshold
        if stranded:
            search_kwargs['stranded'] = True

        if no_wait:
            job_id = client.submit_search(query, **search_kwargs)
            console.print(f"Search submitted with job ID: [bold]{job_id}[/bold]")
            console.print(f"Check status with: [bold]malva_client status {job_id}[/bold]")
            return

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console
        ) as progress:
            task = progress.add_task("Searching...", total=None)
            results = client.search(query, max_wait=max_wait, **search_kwargs)

        if not no_enrich and hasattr(results, 'enrich_with_metadata'):
            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                console=console
            ) as progress:
                task = progress.add_task("Enriching with metadata...", total=None)
                results.enrich_with_metadata()

        if output_format == 'table':
            _display_results_table(results)

        if output:
            _save_results_to_file(results, output, output_format)
            
    except AuthenticationError as e:
        console.print(f"[red]Authentication error: {e}[/red]")
        console.print("Run [bold]malva_client login[/bold] or get an API token from your profile page.")
        sys.exit(1)
    except QuotaExceededError as e:
        console.print(f"[red]Quota exceeded: {e}[/red]")
        if e.quota_info:
            console.print(f"Account type: {e.quota_info.get('account_type', 'unknown')}")
        sys.exit(1)
    except MalvaAPIError as e:
        console.print(f"[red]Error: {e}[/red]")
        sys.exit(1)

@main.command()
@click.pass_context
def login(ctx):
    """Interactive login to Malva"""
    config = ctx.obj['config']
    from malva_client.auth import login as auth_login
    from malva_client.exceptions import AuthenticationError
    
    try:
        if not config.server_url:
            server_url = Prompt.ask(
                "Enter Malva server URL",
                default="https://malva.mdc-berlin.de"
            )
            config.server_url = server_url
        
        console.print(f"Logging in to: [bold]{config.server_url}[/bold]")
        
        token = auth_login(config.server_url)
        config.api_token = token
        config.save()
        
        console.print("[green]Login successful![/green]")
        console.print(f"Configuration saved to: {config.config_path}")
        
    except AuthenticationError as e:
        console.print(f"[red]Login failed: {e}[/red]")
        sys.exit(1)
    except KeyboardInterrupt:
        console.print("\n[yellow]Login cancelled[/yellow]")
        sys.exit(1)

@main.command()
@click.option('--key', help='Configuration key to set')
@click.option('--value', help='Configuration value')
@click.option('--server', help='Server URL')
@click.option('--token', help='API token')
@click.option('--show', is_flag=True, help='Show current configuration')
@click.option('--reset', is_flag=True, help='Reset configuration to defaults')
@click.pass_context  
def config_cmd(ctx, key, value, server, token, show, reset):
    """Manage configuration"""
    config = ctx.obj['config']
    
    if reset:
        if Confirm.ask("Reset configuration to defaults?"):
            config.reset_to_defaults()
            config.save()
            console.print("[green]Configuration reset to defaults[/green]")
        return
    
    if show or (not key and not server and not token):
        _show_config(config)
        return

    changed = False
    
    if server:
        config.server_url = server
        changed = True
        console.print(f"Set server_url to: [bold]{server}[/bold]")
    
    if token:
        config.api_token = token
        changed = True
        console.print("Set api_token: [bold]***[/bold]")
    
    if key and value:
        try:
            config.set_value(key, value)
            changed = True
            console.print(f"Set {key} to: [bold]{value}[/bold]")
        except ValueError as e:
            console.print(f"[red]Error: {e}[/red]")
            return
    
    if changed:
        config.save()
        console.print(f"Configuration saved to: {config.config_path}")
    
    # Validate configuration
    errors = config.validate()
    if errors:
        console.print("\n[yellow]Configuration warnings:[/yellow]")
        for error in errors:
            console.print(f"  - {error}")

main.add_command(config_cmd, name='config')

@main.command()
@click.argument('job_id')
@click.pass_context
def status(ctx, job_id):
    """Check search status"""
    config = ctx.obj['config']
    from malva_client.client import MalvaClient
    from malva_client.exceptions import MalvaAPIError, AuthenticationError
    
    try:
        client = MalvaClient(
            base_url=config.server_url,
            api_token=config.api_token,
            verify_ssl=config.verify_ssl
        )
        
        status_data = client.get_job_status(job_id)

        table = Table(title=f"Job Status: {job_id}")
        table.add_column("Property", style="cyan")
        table.add_column("Value", style="white")
        
        table.add_row("Job ID", job_id)
        table.add_row("Status", _colorize_status(status_data.get('status', 'unknown')))
        
        if 'query' in status_data:
            table.add_row("Query", status_data['query'])
        
        if status_data.get('status') == 'completed':
            results = status_data.get('results', {})
            total_cells = sum(r.get('ncells', 0) for r in results.values())
            table.add_row("Total Cells", f"{total_cells:,}")
            
        if status_data.get('error'):
            table.add_row("Error", f"[red]{status_data['error']}[/red]")
        
        console.print(table)
        
        if status_data.get('status') == 'completed':
            console.print(f"\nGet results with: [bold]malva_client results {job_id}[/bold]")
        
    except AuthenticationError as e:
        console.print(f"[red]Authentication error: {e}[/red]")
        sys.exit(1)
    except MalvaAPIError as e:
        console.print(f"[red]Error: {e}[/red]")
        sys.exit(1)

@main.command()
@click.argument('job_id')
@click.option('--output', '-o', type=click.Path(), help='Output file path')
@click.option('--format', 'output_format', type=click.Choice(['csv', 'json', 'excel', 'table']),
              default='table', help='Output format')
@click.pass_context
def results(ctx, job_id, output, output_format):
    """Get results from previous search"""
    config = ctx.obj['config']
    from malva_client.client import MalvaClient
    from malva_client.exceptions import MalvaAPIError, AuthenticationError, SearchError
    
    try:
        client = MalvaClient(
            base_url=config.server_url,
            api_token=config.api_token,
            verify_ssl=config.verify_ssl
        )
        
        console.print(f"Retrieving results for job: [bold]{job_id}[/bold]")
        results = client.get_job_results(job_id)
        
        if output_format == 'table' and not output:
            _display_results_table(results)
        
        if output:
            _save_results_to_file(results, output, output_format)
        elif output_format != 'table':
            extension = 'csv' if output_format == 'csv' else 'json' if output_format == 'json' else 'xlsx'
            output = f"malva_results_{job_id}.{extension}"
            _save_results_to_file(results, output, output_format)
        
    except AuthenticationError as e:
        console.print(f"[red]Authentication error: {e}[/red]")
        sys.exit(1)
    except SearchError as e:
        console.print(f"[red]Search error: {e}[/red]")
        sys.exit(1)
    except MalvaAPIError as e:
        console.print(f"[red]Error: {e}[/red]")
        sys.exit(1)

@main.command()
@click.option('--limit', default=20, help='Number of searches to show')
@click.option('--export', type=click.Path(), help='Export history to file')
@click.pass_context
def history(ctx, limit, export):
    """View search history"""
    from malva_client.storage import get_search_history
    
    config = ctx.obj['config']
    
    try:
        history = get_search_history(limit=limit, server_url=config.server_url)
        
        if not history:
            console.print("[yellow]No search history found[/yellow]")
            return
        
        table = Table(title="Search History")
        table.add_column("Date", style="cyan")
        table.add_column("Query", style="white", max_width=40)
        table.add_column("Status", style="white")
        table.add_column("Results", style="white")
        table.add_column("Job ID", style="dim", max_width=20)
        
        for entry in history:
            timestamp = entry['timestamp'][:19].replace('T', ' ')  # Format datetime
            query = entry['query'][:37] + "..." if len(entry['query']) > 40 else entry['query']
            status = _colorize_status(entry['status'])
            results = f"{entry.get('results_count', 0):,}" if entry.get('results_count') else "-"
            job_id = entry.get('job_id', '')[:17] + "..." if entry.get('job_id') and len(entry.get('job_id', '')) > 20 else entry.get('job_id', '')
            
            table.add_row(timestamp, query, status, results, job_id)
        
        console.print(table)
        
        if export:
            from malva_client.storage import get_storage
            storage = get_storage()
            format_type = 'json' if export.endswith('.json') else 'csv'
            storage.export_history(Path(export), format_type)
        
    except Exception as e:
        console.print(f"[red]Error retrieving history: {e}[/red]")
        sys.exit(1)

@main.command()
@click.pass_context
def quota(ctx):
    """Check API quota status"""
    config = ctx.obj['config']
    from malva_client.client import MalvaClient
    from malva_client.exceptions import MalvaAPIError, AuthenticationError
    
    try:
        client = MalvaClient(
            base_url=config.server_url,
            api_token=config.api_token,
            verify_ssl=config.verify_ssl
        )
        
        quota_data = client.get_quota_status()
        
        table = Table(title="API Quota Status")
        table.add_column("Property", style="cyan")
        table.add_column("Value", style="white")
        
        table.add_row("Account Type", quota_data.get('account_type', 'unknown'))
        table.add_row("Searches Used", f"{quota_data.get('searches_used', 0):,}")
        table.add_row("Search Limit", f"{quota_data.get('search_limit', 0):,}")
        
        remaining = quota_data.get('search_limit', 0) - quota_data.get('searches_used', 0)
        table.add_row("Remaining", f"{remaining:,}")
        
        if 'reset_time' in quota_data:
            table.add_row("Reset Time", quota_data['reset_time'])
        
        console.print(table)

        if remaining <= 2:
            console.print("[yellow]Warning: Low quota remaining![/yellow]")
        
    except AuthenticationError as e:
        console.print(f"[red]Authentication error: {e}[/red]")
        console.print("Run [bold]malva_client login[/bold] to authenticate.")
        sys.exit(1)
    except MalvaAPIError as e:
        console.print(f"[red]Error: {e}[/red]")
        sys.exit(1)

@main.command()
@click.pass_context
def logout(ctx):
    """Logout and remove stored credentials"""
    config = ctx.obj['config']
    from malva_client.auth import logout as auth_logout
    
    if auth_logout(config.server_url or ""):
        # Also clear config file
        config.api_token = None
        config.save()
        console.print("[green]Logged out successfully[/green]")
    else:
        console.print("[yellow]No stored credentials found[/yellow]")

@main.command()
@click.option('--days', default=30, help='Days of data to keep')
@click.pass_context
def cleanup(ctx, days):
    """Clean up old search data"""
    from malva_client.storage import cleanup_old_data
    
    if Confirm.ask(f"Delete search data older than {days} days?"):
        cleanup_old_data(days)
        console.print(f"[green]Cleaned up data older than {days} days[/green]")

def _display_results_table(results):
    """Display search results in a table"""
    if hasattr(results, 'results') and results.results:
        table = Table(title="Search Results")
        table.add_column("Gene/Sequence", style="cyan")
        table.add_column("Cells Found", style="white")
        
        total_cells = 0
        for key, result in results.results.items():
            if key != "_sample_metadata":
                ncells = result.get('ncells', 0)
                table.add_row(key, f"{ncells:,}")
                total_cells += ncells
        
        console.print(table)
        console.print(f"\nTotal cells found: [bold]{total_cells:,}[/bold]")
    else:
        console.print("[yellow]No results found[/yellow]")

def _save_results_to_file(results, output_path, format_type):
    """Save results to file"""
    output_path = Path(output_path)
    
    try:
        df = results.to_pandas()
        
        if format_type == 'csv':
            df.to_csv(output_path, index=False)
        elif format_type == 'json':
            df.to_json(output_path, orient='records', indent=2)
        elif format_type == 'excel':
            try:
                df.to_excel(output_path, index=False)
            except ImportError:
                console.print("[red]Excel export requires openpyxl. Install with: pip install openpyxl[/red]")
                return
        
        console.print(f"Results saved to: [bold]{output_path}[/bold] ({len(df)} rows)")
        
    except Exception as e:
        console.print(f"[red]Error saving file: {e}[/red]")

def _show_config(config):
    """Display current configuration"""
    config_data = config.show()
    
    table = Table(title="Current Configuration")
    table.add_column("Setting", style="cyan")
    table.add_column("Value", style="white")
    
    for key, value in config_data.items():
        if key == 'api_token' and value:
            value = "***" if value != '***' else value
        table.add_row(key.replace('_', ' ').title(), str(value) if value else "[dim]Not set[/dim]")
    
    console.print(table)

    errors = config.validate()
    if errors:
        console.print("\n[red]Configuration Issues:[/red]")
        for error in errors:
            console.print(f"  - {error}")

def _colorize_status(status):
    """Add color to status text"""
    status_colors = {
        'completed': '[green]completed[/green]',
        'pending': '[yellow]pending[/yellow]',
        'running': '[blue]running[/blue]',
        'error': '[red]error[/red]',
        'failed': '[red]failed[/red]',
        'cancelled': '[dim]cancelled[/dim]'
    }
    return status_colors.get(status, status)

if __name__ == '__main__':
    main()