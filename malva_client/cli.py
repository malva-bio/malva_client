import sys
import click
from rich.console import Console
from rich.table import Table

console = Console()

@click.group()
@click.version_option()
@click.option('--config', type=click.Path(), help='Path to configuration file')
@click.option('--server', help='Malva server URL')
@click.option('--token', help='API token for authentication')
@click.pass_context
def main(ctx, config, server, token):
    from malva_client.config import Config
    """Malva CLI - Search genomic datasets from the command line."""
    ctx.ensure_object(dict)
    
    # Load configuration
    config_obj = Config.load(config_path=config)
    
    # Override config with CLI arguments
    if server:
        config_obj.server_url = server
    if token:
        config_obj.api_token = token
    
    ctx.obj['config'] = config_obj

@main.command()
@click.argument('query')
@click.option('--output', '-o', type=click.Path(), help='Output file path')
@click.option('--format', 'output_format', type=click.Choice(['csv', 'json', 'table']), 
              default='table', help='Output format')
@click.pass_context
def search(ctx, query, output, output_format):
    """Search genomic datasets."""
    config = ctx.obj['config']
    from malva_client.client import MalvaClient
    from malva_client.exceptions import MalvaAPIError, AuthenticationError, QuotaExceededError
    
    try:
        client = MalvaClient(
            base_url=config.server_url,
            api_token=config.api_token,
            verify_ssl=config.verify_ssl
        )
        
        console.print(f"Searching for: [bold]{query}[/bold]")
        results = client.search(query)
        
        # Display results based on format
        if output_format == 'table':
            table = Table(title="Search Results")
            table.add_column("Gene/Sequence", style="cyan")
            table.add_column("Cells Found", style="white")
            
            for key, result in results.results.items():
                table.add_row(key, f"{result.get('ncells', 0):,}")
            
            console.print(table)
        
        if output:
            # Save results to file
            df = results.to_dataframe()
            if output_format == 'csv':
                df.to_csv(output, index=False)
            elif output_format == 'json':
                df.to_json(output, orient='records', indent=2)
            console.print(f"Results saved to: [bold]{output}[/bold]")
            
    except AuthenticationError as e:
        console.print(f"[red]Authentication error: {e}[/red]")
        console.print("Get an API token from your profile page at the Malva website.")
        sys.exit(1)
    except QuotaExceededError as e:
        console.print(f"[red]Quota exceeded: {e}[/red]")
        if e.quota_info:
            console.print(f"Account type: {e.quota_info.get('account_type', 'unknown')}")
        sys.exit(1)
    except MalvaAPIError as e:
        console.print(f"[red]Error: {e}[/red]")
        sys.exit(1)