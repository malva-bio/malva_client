import json
import time
import gzip
import io
import zipfile
import requests
import numpy as np
import pandas as pd
from typing import Dict, List, Any, Optional, Union, Tuple
from pathlib import Path
import logging
from urllib.parse import urljoin
import h5py
import anndata as ad

from malva_client.exceptions import MalvaAPIError, AuthenticationError, SearchError, QuotaExceededError

logger = logging.getLogger(__name__)


def _expr_data_from_columnar(expr_col: dict) -> dict:
    """Expand columnar (_fmt: col) expression_data to row-major format.

    The server encodes aggregated expression data in a compact columnar layout
    (si, ci, v0-v4 as parallel arrays) to minimise JSON payload. This function
    reconstructs the row-major 'data' list expected by _convert_to_dataframe.
    """
    si = expr_col.get('si', [])
    n = len(si)
    if n == 0:
        return {
            'data': [],
            'samples': expr_col.get('samples', []),
            'cell_types': expr_col.get('cell_types', []),
            'columns': expr_col.get('columns', []),
            'celltype_sample_counts': expr_col.get('celltype_sample_counts', {}),
        }
    ci = expr_col.get('ci', [])
    v0 = expr_col.get('v0', [])
    v1 = expr_col.get('v1', [])
    v2 = expr_col.get('v2', [])
    v3 = expr_col.get('v3', [])
    v4 = expr_col.get('v4', [])
    return {
        'data': [[si[j], ci[j], v0[j], v1[j], v2[j], v3[j], v4[j]] for j in range(n)],
        'samples': expr_col.get('samples', []),
        'cell_types': expr_col.get('cell_types', []),
        'columns': expr_col.get('columns', []),
        'celltype_sample_counts': expr_col.get('celltype_sample_counts', {}),
    }


class MalvaDataFrame:
    """
    Wrapper around pandas DataFrame with sample metadata enrichment and analysis methods
    """

    FIELD_MAPPING = {
        'specimen_from_organism.organ': 'organ',
        'specimen_from_organism.organ_part': 'organ_part', 
        'specimen_from_organism.diseases': 'disease',
        'project.project_core.project_short_name': 'study',
        'project.project_core.project_title': 'study_title',
        'project.contributors.laboratory': 'laboratory',
        'genus_species': 'species',
        'development_stage': 'development_stage',
        'protocol_harmonized': 'protocol',
        'sex': 'sex',
        'age': 'age',
        'diseases': 'disease_alt',
        'uuid': 'sample_uuid'
    }
    
    def __init__(self, expr_df: pd.DataFrame, client: 'MalvaClient', sample_metadata: Dict[int, Dict] = None):
        self.client = client
        self._df = expr_df.copy()
        self._sample_metadata = sample_metadata or {}
        self._enrich_with_metadata()
    
    def _enrich_with_metadata(self):
        """Enrich the dataframe with sample metadata using simplified field names"""
        if not self._sample_metadata:
            return
            
        # Create metadata DataFrame with renamed columns
        metadata_records = []
        for sample_id, metadata in self._sample_metadata.items():
            record = {'sample_id': sample_id}
            
            # Add metadata with both original and mapped field names
            for original_field, value in metadata.items():
                # Add original field
                record[original_field] = value
                
                # Add mapped field if mapping exists
                if original_field in self.FIELD_MAPPING:
                    mapped_field = self.FIELD_MAPPING[original_field]
                    record[mapped_field] = value
            
            metadata_records.append(record)
        
        if metadata_records:
            metadata_df = pd.DataFrame(metadata_records)
            
            # Merge with expression data
            self._df = self._df.merge(
                metadata_df, 
                on='sample_id', 
                how='left',
                suffixes=('', '_metadata')
            )
            
            self._convert_to_categorical()

        self._extract_celltype_sample_counts()

    def _extract_celltype_sample_counts(self):
        """Extract celltype_sample_counts from the raw data if available"""
        try:
            if hasattr(self, 'raw_data') and self.raw_data:
                results = self.raw_data.get('results', {})
                for gene_key, gene_data in results.items():
                    if isinstance(gene_data, dict) and 'expression_data' in gene_data:
                        expr_data = gene_data['expression_data']
                        if 'celltype_sample_counts' in expr_data:
                            self._celltype_sample_counts = expr_data['celltype_sample_counts']
                            print(f"✓ Extracted sample counts per cell type: {len(self._celltype_sample_counts)} cell types")
                            return
            
            # If we can't find it, we'll have to live with the fallback method
            print("ℹ️  Could not extract true sample counts per cell type")
            
        except Exception as e:
            print(f"⚠️  Error extracting sample counts: {e}")

    def _convert_to_categorical(self):
        """Convert string/object columns to categorical for faster filtering"""
        categorical_candidates = [
            'organ', 'organ_part', 'disease', 'species', 'development_stage', 
            'sex', 'study', 'study_title', 'laboratory', 'protocol', 'cell_type'
        ]
        
        for col in categorical_candidates:
            if col in self._df.columns and self._df[col].dtype == 'object':
                # Only convert if the column has a reasonable number of unique values
                # (avoid converting columns with mostly unique values)
                unique_ratio = self._df[col].nunique() / len(self._df)
                if unique_ratio < 0.5:  # Less than 50% unique values
                    self._df[col] = self._df[col].astype('category')
    
    @property 
    def df(self) -> pd.DataFrame:
        """Access the underlying pandas DataFrame"""
        return self._df
    
    def __getattr__(self, name):
        """Delegate pandas DataFrame methods"""
        return getattr(self._df, name)
    
    def __getitem__(self, key):
        """Delegate pandas DataFrame indexing"""
        return self._df[key]
    
    def __len__(self):
        return len(self._df)
    
    def filter_by(self, **kwargs) -> 'MalvaDataFrame':
        """
        Filter data by any combination of metadata fields (optimized for large datasets)
        Uses simplified field names (e.g., 'organ' instead of 'specimen_from_organism.organ')
        
        Args:
            **kwargs: Field=value pairs for filtering
            
        Returns:
            New MalvaDataFrame with filtered data
            
        Example:
            df.filter_by(organ='brain', disease='normal', cell_type='neuron')
            df.filter_by(species='Homo sapiens', study='BrainDepressiveDisorder')
        """
        if not kwargs:
            return MalvaDataFrame(self._df.copy(), self.client, self._sample_metadata)
        
        # Start with all rows as True
        mask = pd.Series(True, index=self._df.index)
        
        # Create reverse mapping for convenience
        reverse_mapping = {v: k for k, v in self.FIELD_MAPPING.items()}
        
        for field, value in kwargs.items():
            # Check if field exists directly
            actual_field = field
            if field not in self._df.columns:
                # Try mapped field
                if field in reverse_mapping and reverse_mapping[field] in self._df.columns:
                    actual_field = reverse_mapping[field]
                else:
                    logging.warning(f"Field '{field}' not found. Available fields: {self.available_filter_fields()}")
                    continue
            
            # Apply filter using boolean indexing (much faster)
            if isinstance(value, (list, tuple)):
                # For categorical columns, this is very fast
                field_mask = self._df[actual_field].isin(value)
            else:
                # For categorical columns, this is also very fast
                field_mask = self._df[actual_field] == value
            
            # Combine with existing mask using AND operation
            mask = mask & field_mask
        
        # Apply the final mask in one operation
        filtered_df = self._df[mask].copy()
        
        return MalvaDataFrame(filtered_df, self.client, self._sample_metadata)
    
    EXPRESSION_COLUMN_ALIASES = {
        'norm_expr': 'rel',
        'kpt_expr': 'exp',
        'raw_expr': 'exp',
        'fraction_positive': 'pct',
        'pct_positive': 'pct',
        'raw_kmer_mean': 'raw_kmers',
    }

    def _resolve_expression_column(self, expr_column: str) -> str:
        """Resolve legacy expression column names to the canonical schema."""
        if expr_column in self._df.columns:
            return expr_column
        alias = self.EXPRESSION_COLUMN_ALIASES.get(expr_column)
        if alias and alias in self._df.columns:
            return alias
        raise ValueError(f"Expression column '{expr_column}' not found")

    def aggregate_by(self, group_by: Union[str, List[str]],
                    agg_func: str = 'mean',
                    expr_column: str = 'rel') -> pd.DataFrame:
        """
        Aggregate expression data by specified grouping variables
        
        Args:
            group_by: Column name(s) to group by
            agg_func: Aggregation function ('mean', 'median', 'sum', 'count', 'std')
            expr_column: Expression column to aggregate
            
        Returns:
            DataFrame with aggregated results
            
        Example:
            df.aggregate_by('cell_type')
            df.aggregate_by(['organ', 'cell_type'])
        """
        expr_column = self._resolve_expression_column(expr_column)
        
        group_cols = [group_by] if isinstance(group_by, str) else group_by
        
        # Resolve field names
        resolved_group_cols = []
        reverse_mapping = {v: k for k, v in self.FIELD_MAPPING.items()}
        
        for col in group_cols:
            if col in self._df.columns:
                resolved_group_cols.append(col)
            elif col in reverse_mapping and reverse_mapping[col] in self._df.columns:
                resolved_group_cols.append(reverse_mapping[col])
            else:
                raise ValueError(f"Grouping column '{col}' not found")
        
        # Simple aggregation approach
        grouped = self._df.groupby(resolved_group_cols)
        
        # Create result DataFrame step by step
        if agg_func == 'mean':
            expr_agg = grouped[expr_column].mean()
        elif agg_func == 'median':
            expr_agg = grouped[expr_column].median()
        elif agg_func == 'sum':
            expr_agg = grouped[expr_column].sum()
        elif agg_func == 'count':
            expr_agg = grouped[expr_column].count()
        elif agg_func == 'std':
            expr_agg = grouped[expr_column].std()
        elif agg_func == 'min':
            expr_agg = grouped[expr_column].min()
        elif agg_func == 'max':
            expr_agg = grouped[expr_column].max()
        else:
            raise ValueError(f"Unsupported aggregation function: {agg_func}")
        
        # Build result DataFrame
        result = pd.DataFrame({
            f"{agg_func}_{expr_column}": expr_agg,
            'n_observations': grouped.size(),
        })
        
        # Add additional metrics if columns exist
        if 'cell_count' in self._df.columns:
            result['total_cells'] = grouped['cell_count'].sum()
        else:
            result['total_cells'] = grouped.size()  # Use observation count as proxy
        
        if 'sample_id' in self._df.columns:
            result['n_samples'] = grouped['sample_id'].nunique()
        else:
            result['n_samples'] = 1  # Default
        
        # Reset index and sort
        result = result.reset_index()
        result = result.sort_values(f"{agg_func}_{expr_column}", ascending=False)
        
        return result.round(6)
    
    def plot_expression_by(self, group_by: str, limit: int = None, sort_by: str = 'mean', 
                      ascending: bool = False, **kwargs):
        """
        Plot expression levels grouped by a metadata field
        
        Args:
            group_by: Column to group by for plotting
            limit: Maximum number of categories to show (shows top N)
            sort_by: How to sort categories ('mean', 'median', 'count', 'alphabetical')
            ascending: Whether to sort in ascending order (False shows highest first)
            **kwargs: Additional arguments passed to matplotlib
        """
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
        except ImportError:
            raise ImportError("matplotlib and seaborn required for plotting. Install with: pip install matplotlib seaborn")
        
        if group_by not in self._df.columns:
            raise ValueError(f"Column '{group_by}' not found")
        
        expr_column = self._resolve_expression_column('rel')
        
        # Prepare data for plotting
        plot_data = self._df.copy()
        
        # Calculate sorting metrics for each category
        if sort_by != 'alphabetical':
            if sort_by == 'mean':
                sort_values = plot_data.groupby(group_by)[expr_column].mean()
            elif sort_by == 'median':
                sort_values = plot_data.groupby(group_by)[expr_column].median()
            elif sort_by == 'count':
                sort_values = plot_data.groupby(group_by)[expr_column].count()
            else:
                raise ValueError("sort_by must be one of: 'mean', 'median', 'count', 'alphabetical'")
            
            # Sort categories
            sorted_categories = sort_values.sort_values(ascending=ascending).index.tolist()
        else:
            # Alphabetical sorting
            sorted_categories = sorted(plot_data[group_by].dropna().unique())
            if not ascending:
                sorted_categories = sorted_categories[::-1]
        
        # Apply limit if specified
        if limit is not None and limit > 0:
            sorted_categories = sorted_categories[:limit]
            plot_data = plot_data[plot_data[group_by].isin(sorted_categories)]
        
        # Create the plot
        plt.figure(figsize=kwargs.get('figsize', (12, 6)))
        
        # Remove kwargs that shouldn't go to seaborn
        plot_kwargs = {k: v for k, v in kwargs.items() if k not in ['figsize']}
        
        # Create box plot with ordered categories
        sns.boxplot(data=plot_data, x=group_by, y=expr_column, 
                    order=sorted_categories, **plot_kwargs)
        
        # Customize plot
        plt.xticks(rotation=45, ha='right')
        plt.xlabel(group_by.replace('_', ' ').title())
        plt.ylabel('Normalized Expression')
        
        # Create informative title
        title_parts = [f'Expression Distribution by {group_by.replace("_", " ").title()}']
        if limit:
            title_parts.append(f'(Top {limit}')
            if sort_by != 'alphabetical':
                title_parts.append(f'by {sort_by})')
            else:
                title_parts.append('alphabetically)')
        title = ' '.join(title_parts)
        
        plt.title(title)
        plt.tight_layout()
        
        # Add summary statistics as text if there are few categories
        if len(sorted_categories) <= 10:
            stats_text = []
            for i, category in enumerate(sorted_categories):
                cat_data = plot_data[plot_data[group_by] == category][expr_column]
                if not cat_data.empty:
                    mean_val = cat_data.mean()
                    count = len(cat_data)
                    stats_text.append(f'{category}: μ={mean_val:.3f} (n={count})')
            
            # Add text box with statistics
            if stats_text and len(stats_text) <= 8:  # Only if not too many
                textstr = '\n'.join(stats_text)
                props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
                plt.gca().text(0.02, 0.98, textstr, transform=plt.gca().transAxes, 
                            fontsize=8, verticalalignment='top', bbox=props)
        
        return plt.gcf()
    
    def plot_expression_summary(self, group_by: str, limit: int = 10, 
                                plot_type: str = 'box', **kwargs):
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
            import numpy as np
        except ImportError:
            raise ImportError("matplotlib and seaborn required for plotting")
        
        if group_by not in self._df.columns:
            raise ValueError(f"Column '{group_by}' not found")
        
        # Get top categories by mean expression
        expr_column = self._resolve_expression_column('rel')
        top_categories = (self._df.groupby(group_by)[expr_column]
                        .mean().nlargest(limit).index.tolist())
        
        plot_data = self._df[self._df[group_by].isin(top_categories)]
        
        # Create 1x3 subplot layout
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        fig.suptitle(f'Expression Analysis by {group_by.replace("_", " ").title()}', 
                    fontsize=16)
        
        # 1. Box plot (or other plot type)
        if plot_type == 'box':
            sns.boxplot(data=plot_data, x=group_by, y=expr_column, 
                        order=top_categories, ax=axes[0])
            axes[0].set_title('Expression Distribution (Box Plot)')
        elif plot_type == 'violin':
            sns.violinplot(data=plot_data, x=group_by, y=expr_column, 
                        order=top_categories, ax=axes[0])
            axes[0].set_title('Expression Distribution (Violin Plot)')
        elif plot_type == 'strip':
            sns.stripplot(data=plot_data, x=group_by, y=expr_column, 
                        order=top_categories, ax=axes[0])
            axes[0].set_title('Expression Distribution (Strip Plot)')
        else:  # bar plot
            means = plot_data.groupby(group_by)[expr_column].agg(['mean', 'std'])
            means = means.reindex(top_categories)
            means.plot(kind='bar', y='mean', yerr='std', ax=axes[0], 
                    color='skyblue', capsize=4)
            axes[0].set_title('Mean Expression ± SD')
            axes[0].legend().remove()
        
        axes[0].tick_params(axis='x', rotation=45)
        axes[0].set_xlabel(group_by.replace('_', ' ').title())
        axes[0].set_ylabel('Normalized Expression')
        
        # 2. Number of samples per category - FIXED TO USE SERVER DATA
        sample_counts_per_category = self._get_true_sample_counts(group_by, top_categories)
        sample_counts = pd.Series(sample_counts_per_category)
        
        sample_counts.plot(kind='bar', ax=axes[1], color='lightcoral')
        axes[1].set_title('Number of Samples')
        axes[1].set_xlabel(group_by.replace('_', ' ').title())
        axes[1].set_ylabel('Sample Count')
        axes[1].tick_params(axis='x', rotation=45)
        
        # Add value labels on bars
        for i, v in enumerate(sample_counts.values):
            axes[1].text(i, v + 0.01 * max(sample_counts.values), str(v), 
                        ha='center', va='bottom', fontsize=9)
        
        # 3. Total number of cells per category
        if 'cell_count' in plot_data.columns:
            cell_counts = plot_data.groupby(group_by)['cell_count'].sum().reindex(top_categories)
            ylabel = 'Total Cells'
            title = 'Number of Cells'
        else:
            cell_counts = plot_data.groupby(group_by).size().reindex(top_categories)
            ylabel = 'Observations'
            title = 'Number of Observations'
        
        cell_counts.plot(kind='bar', ax=axes[2], color='lightgreen')
        axes[2].set_title(title)
        axes[2].set_xlabel(group_by.replace('_', ' ').title())
        axes[2].set_ylabel(ylabel)
        axes[2].tick_params(axis='x', rotation=45)
        
        # Add value labels on bars (format large numbers)
        for i, v in enumerate(cell_counts.values):
            if v >= 1000:
                label = f'{v/1000:.1f}K' if v < 1000000 else f'{v/1000000:.1f}M'
            else:
                label = str(int(v))
            axes[2].text(i, v + 0.01 * max(cell_counts.values), label, 
                        ha='center', va='bottom', fontsize=9)
        
        # Format y-axis for large numbers
        if max(cell_counts.values) >= 1000:
            axes[2].yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1000:.0f}K' if x >= 1000 else f'{x:.0f}'))
        
        plt.tight_layout()
        
        # Print summary statistics - FIXED WITH CORRECT SAMPLE COUNTS
        print(f"\n📊 Summary for {group_by.replace('_', ' ').title()}:")
        print("-" * 50)
        for category in top_categories:
            cat_data = plot_data[plot_data[group_by] == category]
            n_samples = sample_counts_per_category.get(category, cat_data['sample_id'].nunique())
            n_cells = cat_data['cell_count'].sum() if 'cell_count' in cat_data.columns else len(cat_data)
            mean_expr = cat_data[expr_column].mean()
            print(f"{category}: {n_samples} samples, {n_cells:,} cells, μ={mean_expr:.3f}")
        
        return fig

    def _get_true_sample_counts(self, group_by: str, categories: List[str]) -> Dict[str, int]:
        """Extract true sample counts from server-provided data"""
        
        # Method 1: Check if we have celltype_sample_counts in the data structure
        if hasattr(self, '_sample_metadata') and group_by == 'cell_type':
            # Look for the celltype_sample_counts in the original search results
            # This is a bit hacky but necessary given the current data structure
            
            # Try to find expression_data with celltype_sample_counts
            sample_counts = {}
            
            # Check if any columns contain serialized data with our counts
            for col in self._df.columns:
                if 'expression_data' in str(col).lower():
                    continue  # Skip for now, this is complex to extract
            
            # Alternative: Check if we stored this somewhere during enrichment
            if hasattr(self, '_celltype_sample_counts'):
                for category in categories:
                    sample_counts[category] = self._celltype_sample_counts.get(category, 0)
                return sample_counts
        
        # Method 2: Fallback - count unique samples per category (this gives wrong results but is better than nothing)
        sample_counts = {}
        for category in categories:
            cat_data = self._df[self._df[group_by] == category]
            sample_counts[category] = cat_data['sample_id'].nunique()
        
        return sample_counts

    def to_pandas(self) -> pd.DataFrame:
        """Convert to regular pandas DataFrame"""
        return self._df.copy()
    
    def available_fields(self) -> Dict[str, List[str]]:
        """
        Get categorized list of available fields for filtering/grouping
        
        Returns:
            Dictionary with categorized field names
        """
        all_columns = list(self._df.columns)
        
        # Categorize fields
        categories = {
            'expression': [col for col in all_columns if col in ['rel', 'exp', 'pct', 'raw_kmers', 'cell_count', 'sample_idx', 'cell_type_idx']],
            'basic': [col for col in all_columns if col in ['sample_id', 'cell_type', 'gene_sequence']],
            'biological': [col for col in all_columns if col in ['organ', 'organ_part', 'disease', 'species', 'development_stage', 'sex', 'age']],
            'study': [col for col in all_columns if col in ['study', 'study_title', 'laboratory', 'protocol']],
            'technical': [col for col in all_columns if col in ['chunk_id', 'internal_sample_id', 'num_cells', 'num_data', 'num_kmers', 'sample_uuid']],
            'original_fields': [col for col in all_columns if '.' in col]  # Keep original nested field names
        }
        
        # Remove empty categories
        return {k: v for k, v in categories.items() if v}
    
    def available_filter_fields(self) -> List[str]:
        """Get simplified list of commonly used fields for filtering"""
        common_fields = [
            'organ', 'disease', 'species', 'development_stage', 'sex', 'age',
            'study', 'laboratory', 'protocol', 'cell_type'
        ]
        return [field for field in common_fields if field in self._df.columns]
    
    def field_info(self) -> pd.DataFrame:
        """
        Get detailed information about available fields
        
        Returns:
            DataFrame with field names, types, unique values count, and examples
        """
        info_data = []
        
        for col in self._df.columns:
            unique_vals = self._df[col].dropna().unique()
            
            info_data.append({
                'field_name': col,
                'data_type': str(self._df[col].dtype),
                'non_null_count': self._df[col].count(),
                'unique_values': len(unique_vals),
                'example_values': ', '.join(map(str, unique_vals[:3])) if len(unique_vals) > 0 else 'No data',
                'category': self._categorize_field(col)
            })
        
        return pd.DataFrame(info_data).sort_values(['category', 'field_name'])
    
    def _categorize_field(self, field_name: str) -> str:
        """Categorize a field for better organization"""
        if field_name in ['rel', 'exp', 'pct', 'raw_kmers', 'cell_count', 'sample_idx', 'cell_type_idx']:
            return 'expression'
        elif field_name in ['sample_id', 'cell_type', 'gene_sequence']:
            return 'basic'
        elif field_name in ['organ', 'organ_part', 'disease', 'species', 'development_stage', 'sex', 'age']:
            return 'biological'
        elif field_name in ['study', 'study_title', 'laboratory', 'protocol']:
            return 'study'
        elif '.' in field_name:
            return 'original_nested'
        else:
            return 'technical'
    
    def unique_values(self, field: str, limit: int = 20) -> List:
        """
        Get unique values for a specific field
        
        Args:
            field: Field name to get unique values for
            limit: Maximum number of values to return
            
        Returns:
            Sorted list of unique values
        """
        if field not in self._df.columns:
            # Try to find field in mapping
            reverse_mapping = {v: k for k, v in self.FIELD_MAPPING.items()}
            if field in reverse_mapping:
                field = reverse_mapping[field]
            else:
                available = self.available_filter_fields()
                raise ValueError(f"Field '{field}' not found. Available fields: {available}")
        
        unique_vals = self._df[field].dropna().unique()
        sorted_vals = sorted(unique_vals, key=lambda x: str(x))
        
        if len(sorted_vals) > limit:
            return sorted_vals[:limit]
        return sorted_vals
    
    def value_counts(self, field: str, limit: int = 10) -> pd.Series:
        """Get value counts for a field"""
        if field not in self._df.columns:
            reverse_mapping = {v: k for k, v in self.FIELD_MAPPING.items()}
            if field in reverse_mapping:
                field = reverse_mapping[field]
            else:
                raise ValueError(f"Field '{field}' not found")
        
        return self._df[field].value_counts().head(limit)

class SearchResult(MalvaDataFrame):
    """Container for search results that extends MalvaDataFrame for direct analysis.

    Data loading is lazy: if raw_data contains no 'results' (e.g. it is the
    initial POST response or a lightweight status response), the DataFrame is
    populated on first access to .df by fetching
    /api/expression-data/<job_id>/results from the server.
    """

    def __init__(self, raw_data: Dict[str, Any], client: 'MalvaClient'):
        self.raw_data = raw_data
        self.client = client
        self._enriched = False
        self._job_id = raw_data.get('job_id') if isinstance(raw_data, dict) else None

        # Try to build DataFrame immediately; mark for lazy fetch if no data yet.
        expr_df = self._convert_to_dataframe()
        self._needs_lazy_fetch = expr_df.empty and bool(self._job_id)

        super().__init__(expr_df, client, sample_metadata={})

    def _ensure_loaded(self):
        """Lazily fetch expression data from the server if not yet loaded."""
        if not self._needs_lazy_fetch:
            return
        self._needs_lazy_fetch = False
        aggregate_expression = (
            self.raw_data.get('aggregate_expression')
            if isinstance(self.raw_data, dict)
            else None
        )
        try:
            import time as _time
            logger.info(f"Fetching results for job {self._job_id} ...")
            t0 = _time.time()
            response = self.client._request(
                'GET', f'/api/expression-data/{self._job_id}/results'
            )
            # response may be a dict or a SearchResults wrapper
            if hasattr(response, '_data'):
                response = response._data
            n_genes = len(response.get('results', {})) if isinstance(response, dict) else 0
            logger.info(
                f"Results parsed in {_time.time() - t0:.1f}s: {n_genes} genes"
            )
            if aggregate_expression is not None and isinstance(response, dict):
                response.setdefault('aggregate_expression', aggregate_expression)
            if isinstance(response, dict):
                response.setdefault('job_id', self._job_id)
            self.raw_data = response
            if isinstance(response, dict) and response.get('job_id'):
                self._job_id = response.get('job_id')
            self._df = self._convert_to_dataframe()
            sample_meta = response.get('_sample_metadata', {}) if isinstance(response, dict) else {}
            self._sample_metadata = sample_meta
            self._enrich_with_metadata()
        except Exception as exc:
            logger.error(f"Failed to fetch results for job {self._job_id}: {exc}")

    @property
    def df(self) -> 'pd.DataFrame':
        """Return the expression DataFrame, fetching lazily if needed."""
        self._ensure_loaded()
        return self._df
    
    def _convert_to_dataframe(self) -> pd.DataFrame:
        """Convert raw search results to DataFrame - keep sample-level aggregates"""
        if not self.raw_data or not isinstance(self.raw_data, dict):
            return pd.DataFrame()
        
        results = self.raw_data.get('results', {})
        if not results:
            return pd.DataFrame()
        
        all_data = []
        for gene_seq, result in results.items():
            if gene_seq.startswith('_') or not isinstance(result, dict):
                continue

            # Handle the new compact expression_data format
            if 'expression_data' in result:
                expression_data = result['expression_data']
                # Expand columnar format (_fmt: col) to row-major before reading 'data'.
                # The /api/expression-data/<job_id>/results endpoint returns columnar layout.
                if isinstance(expression_data, dict) and expression_data.get('_fmt') == 'col':
                    expression_data = _expr_data_from_columnar(expression_data)
                data = expression_data.get('data', [])
                columns = expression_data.get('columns', [])
                samples = expression_data.get('samples', [])
                cell_types = expression_data.get('cell_types', [])
                
                if data and columns:
                    # Keep sample x cell-type aggregate rows and expose the same
                    # metric columns used by Expression Explorer display modes.
                    col_idx = {str(name): idx for idx, name in enumerate(columns)}

                    def _value(row, name, default_idx=None, default=0.0):
                        idx = col_idx.get(name, default_idx)
                        if idx is None or idx >= len(row):
                            return default
                        return row[idx]

                    sample_rows = []
                    for row in data:
                        if len(row) >= 5:
                            sample_idx = int(_value(row, 'sample_idx', 0, 0))
                            cell_type_idx = int(_value(row, 'cell_type_idx', 1, 0))
                            norm_expr = float(_value(row, 'norm_expr', 2, 0.0))
                            kpt_expr = float(_value(row, 'kpt_expr', 3, norm_expr))
                            cell_count = int(_value(row, 'cell_count', 4, 0))
                            fraction_positive = float(_value(row, 'fraction_positive', 5, 0.0))
                            raw_kmer_mean = float(_value(row, 'raw_kmer_mean', col_idx.get('raw_expr', 6), kpt_expr))
                            sample_id = samples[sample_idx] if sample_idx < len(samples) else f"sample_{sample_idx}"
                            cell_type = cell_types[cell_type_idx] if cell_type_idx < len(cell_types) else f"celltype_{cell_type_idx}"

                            sample_rows.append({
                                'sample_idx': sample_idx,
                                'cell_type_idx': cell_type_idx,
                                'sample_id': sample_id,
                                'cell_type': cell_type,
                                'gene_sequence': gene_seq,
                                'rel': norm_expr,
                                'exp': kpt_expr,
                                'pct': fraction_positive * 100.0,
                                'raw_kmers': raw_kmer_mean,
                                'cell_count': cell_count,
                                'sample_celltype_id': f"{sample_id}_{cell_type}",
                            })
                    
                    if sample_rows:
                        df = pd.DataFrame(sample_rows)
                        all_data.append(df)
                    continue
            
            # FALLBACK: Handle old format for backward compatibility
            expression_data = result.get('expression_data', {})
            if not expression_data:
                # Try old direct format
                if 'cell' in result and 'expression' in result and 'sample' in result:
                    df = pd.DataFrame({
                        'cell_id': result['cell'],
                        'rel': result['expression'],
                        'sample_id': result['sample'],
                        'gene_sequence': gene_seq
                    })
                    all_data.append(df)
                continue
            
            # Handle old expression_data format
            data = expression_data.get('data', [])
            columns = expression_data.get('columns', [])
            samples = expression_data.get('samples', [])
            cell_types = expression_data.get('cell_types', [])
            
            if not data:
                continue
            
            # Create DataFrame from compact data
            df = pd.DataFrame(data, columns=columns)
            
            # Replace indices with actual values
            if 'sample_idx' in df.columns:
                df['sample_id'] = df['sample_idx'].map(lambda x: samples[x] if x < len(samples) else x)
            if 'cell_type_idx' in df.columns:
                df['cell_type'] = df['cell_type_idx'].map(lambda x: cell_types[x] if x < len(cell_types) else 'Unknown')
            
            df['gene_sequence'] = gene_seq
            all_data.append(df)
        
        if all_data:
            result_df = pd.concat(all_data, ignore_index=True)
            if 'norm_expr' in result_df.columns and 'rel' not in result_df.columns:
                result_df['rel'] = result_df['norm_expr']
            if 'kpt_expr' in result_df.columns and 'exp' not in result_df.columns:
                result_df['exp'] = result_df['kpt_expr']
            elif 'raw_expr' in result_df.columns and 'exp' not in result_df.columns:
                result_df['exp'] = result_df['raw_expr']
            if 'fraction_positive' in result_df.columns and 'pct' not in result_df.columns:
                result_df['pct'] = result_df['fraction_positive'].astype(float) * 100.0
            elif 'pct_positive' in result_df.columns and 'pct' not in result_df.columns:
                result_df['pct'] = result_df['pct_positive']
            if 'raw_kmer_mean' in result_df.columns and 'raw_kmers' not in result_df.columns:
                result_df['raw_kmers'] = result_df['raw_kmer_mean']

            legacy_columns = [
                'norm_expr', 'kpt_expr', 'raw_expr', 'fraction_positive',
                'pct_positive', 'raw_kmer_mean',
            ]
            result_df = result_df.drop(columns=[c for c in legacy_columns if c in result_df.columns])
            if 'cell_type' in result_df.columns:
                result_df['cell_type'] = result_df['cell_type'].astype('category')
            return result_df
        else:
            return pd.DataFrame()
    
    @property
    def job_id(self) -> str:
        """Get the job ID for this search"""
        if isinstance(self.raw_data, dict) and self.raw_data.get('job_id'):
            return str(self.raw_data.get('job_id'))
        return str(self._job_id or '')
    
    @property
    def status(self) -> str:
        """Get the current status of the search"""
        return self.raw_data.get('status', 'unknown')
    
    @property
    def results(self) -> Dict[str, Any]:
        """Get the raw search results"""
        self._ensure_loaded()
        return self.raw_data.get('results', {})

    def __len__(self) -> int:
        self._ensure_loaded()
        return len(self._df)

    @property
    def total_cells(self) -> int:
        """Get total number of cells found"""
        self._ensure_loaded()
        if 'cell_count' in self._df.columns:
            return int(self._df['cell_count'].sum())
        else:
            # Fallback to counting combinations
            return len(self._df)
    
    def enrich_with_metadata(self) -> 'SearchResult':
        """
        Enrich the search results with sample metadata
        
        Returns:
            Self with enriched metadata
        """
        if self._enriched:
            return self

        # Ensure lazy-loaded data is available before checking emptiness
        self._ensure_loaded()

        if self._df.empty:
            return self

        # Get unique sample IDs
        sample_ids = self._df['sample_id'].unique().tolist()
        
        # Fetch sample metadata
        try:
            sample_metadata = self.client.get_sample_metadata(sample_ids)
            total_cells = self._df['cell_count'].sum() if 'cell_count' in self._df.columns else len(self._df)
            print(f"✓ Enriched with metadata for {len(sample_metadata)} samples ({total_cells:,} total cells)")
        except Exception as e:
            logging.warning(f"Could not fetch sample metadata: {e}")
            sample_metadata = {}
        
        # Update the metadata and re-enrich
        self._sample_metadata = sample_metadata
        self._enrich_with_metadata()
        self._enriched = True
        
        return self
    
    def __repr__(self) -> str:
        self._ensure_loaded()
        if self._df.empty:
            return "SearchResult: No data found"
        
        lines = []
        lines.append("🔬 Malva Search Results")
        lines.append("=" * 50)
        
        # FIXED: Always try to get actual cell count first
        if 'cell_count' in self._df.columns and not self._df['cell_count'].isna().all():
            total_cells = int(self._df['cell_count'].sum())
            lines.append(f"📊 Total cells: {total_cells:,}")
            lines.append(f"📊 Sample/cell_type combinations: {len(self._df)}")
        else:
            # Fallback: count rows but label correctly
            lines.append(f"📊 Data points: {len(self._df)}")
            lines.append("📊 Individual cell counts not available")
            
        lines.append(f"🧬 Genes/sequences: {self._df['gene_sequence'].nunique() if 'gene_sequence' in self._df.columns else 'N/A'}")
        lines.append(f"🧪 Samples: {self._df['sample_id'].nunique() if 'sample_id' in self._df.columns else 'N/A'}")
        lines.append(f"🔬 Cell types: {self._df['cell_type'].nunique() if 'cell_type' in self._df.columns else 'N/A'}")
        
        # Expression stats remain the same (these are per sample/cell_type aggregate)
        if 'rel' in self._df.columns:
            lines.append(f"📈 Expression range: {self._df['rel'].min():.3f} - {self._df['rel'].max():.3f}")
            lines.append(f"📊 Mean expression: {self._df['rel'].mean():.3f}")
        
        lines.append("")
        
        # Metadata status
        if self._enriched:
            lines.append("✅ Enriched with sample metadata")
            biological_fields = [f for f in ['organ', 'disease', 'species', 'study'] if f in self._df.columns]
            if biological_fields:
                lines.append(f"🏷️  Available metadata: {', '.join(biological_fields)}")
        else:
            lines.append("ℹ️  Basic expression data only")
            lines.append("💡 Run .enrich_with_metadata() to add sample metadata for filtering by:")
            lines.append("   • Organ, disease, species")
            lines.append("   • Study, laboratory, protocol")
            lines.append("   • Age, sex, development stage")
        
        lines.append("")
        lines.append("🔍 Available methods:")
        lines.append("   • .filter_by(organ='brain', disease='normal')")
        lines.append("   • .aggregate_by('cell_type')")
        lines.append("   • .plot_expression_by('organ')")
        lines.append("   • .available_filter_fields()")
        
        return "\n".join(lines)
    
    def __str__(self) -> str:
        """String representation"""
        return self.__repr__()
    
    def _repr_html_(self) -> str:
        """HTML representation for Jupyter notebooks"""
        self._ensure_loaded()
        if self._df.empty:
            return "<div><strong>SearchResult:</strong> No data found</div>"
        
        html = []
        html.append('<div style="border: 1px solid #ddd; padding: 10px; border-radius: 5px; font-family: monospace;">')
        html.append('<h3 style="margin-top: 0;">🔬 Malva Search Results</h3>')
        
        # Stats table
        html.append('<table style="border-collapse: collapse; margin: 10px 0;">')
        html.append(f'<tr><td><strong>Total cells:</strong></td><td>{len(self._df):,}</td></tr>')
        html.append(f'<tr><td><strong>Genes/sequences:</strong></td><td>{self._df["gene_sequence"].nunique() if "gene_sequence" in self._df.columns else "N/A"}</td></tr>')
        html.append(f'<tr><td><strong>Samples:</strong></td><td>{self._df["sample_id"].nunique() if "sample_id" in self._df.columns else "N/A"}</td></tr>')
        html.append(f'<tr><td><strong>Cell types:</strong></td><td>{self._df["cell_type"].nunique() if "cell_type" in self._df.columns else "N/A"}</td></tr>')
        
        if 'rel' in self._df.columns:
            html.append(f'<tr><td><strong>Expression range:</strong></td><td>{self._df["rel"].min():.3f} - {self._df["rel"].max():.3f}</td></tr>')
            html.append(f'<tr><td><strong>Mean expression:</strong></td><td>{self._df["rel"].mean():.3f}</td></tr>')
        
        html.append('</table>')
        
        # Metadata status
        if self._enriched:
            html.append('<div style="color: green;">✅ <strong>Enriched with sample metadata</strong></div>')
        else:
            html.append('<div style="color: orange;">ℹ️ <strong>Basic expression data only</strong></div>')
            html.append('<div style="margin: 5px 0; font-size: 0.9em;">Run <code>.enrich_with_metadata()</code> to add sample metadata for advanced filtering</div>')
        
        html.append('</div>')
        return ''.join(html)

import pandas as pd
import numpy as np
from typing import Dict, List, Any, Optional, Union, TYPE_CHECKING

if TYPE_CHECKING:
    from malva_client.client import MalvaClient


class CellExpressionMatrixResult:
    """
    Per-cell result returned by MalvaClient.retrieve_cells().

    New results are built directly from search and metadata endpoints. The
    older ZIP-backed constructor remains supported for compatibility.
    """

    def __init__(self, archive_path: Optional[Union[str, Path]] = None,
                 export_record: Optional[Dict[str, Any]] = None,
                 client: Optional['MalvaClient'] = None,
                 *,
                 cells: Optional[pd.DataFrame] = None,
                 features: Optional[pd.DataFrame] = None,
                 matrix_entries: Optional[pd.DataFrame] = None,
                 normalization_factors: Optional[pd.DataFrame] = None,
                 sample_metadata: Optional[pd.DataFrame] = None,
                 barcodes: Optional[pd.DataFrame] = None,
                 job_id: Optional[str] = None,
                 source: str = 'direct'):
        self.archive_path = Path(archive_path) if archive_path is not None else None
        self.export_record = export_record or {}
        self.client = client
        self.export_id = self.export_record.get('export_id')
        self.job_id = job_id or self.export_record.get('job_id')
        self.source = source or self.export_record.get('source', 'direct')
        self.status = self.export_record.get('status', 'completed')

        self._cells_df: Optional[pd.DataFrame] = cells.copy() if cells is not None else None
        self._features_df: Optional[pd.DataFrame] = features.copy() if features is not None else None
        self._normalization_df: Optional[pd.DataFrame] = (
            normalization_factors.copy() if normalization_factors is not None else None
        )
        self._sample_metadata_df: Optional[pd.DataFrame] = (
            sample_metadata.copy() if sample_metadata is not None else None
        )
        self._barcodes_df: Optional[pd.DataFrame] = barcodes.copy() if barcodes is not None else None
        self._matrix_entries_df: Optional[pd.DataFrame] = (
            matrix_entries.copy() if matrix_entries is not None else None
        )

    def _read_gz_tsv(self, member: str) -> pd.DataFrame:
        if self.archive_path is None:
            return pd.DataFrame()

        if not self.archive_path.exists():
            raise FileNotFoundError(f"Archive not found: {self.archive_path}")
        with zipfile.ZipFile(self.archive_path) as zf:
            if member not in zf.namelist():
                return pd.DataFrame()
            with zf.open(member) as raw:
                with gzip.GzipFile(fileobj=raw, mode='rb') as gz:
                    return pd.read_csv(gz, sep='\t')

    @property
    def cells(self) -> pd.DataFrame:
        """Rows of the matrix: row_index, sample_id, cell_id."""
        if self._cells_df is None:
            self._cells_df = self._read_gz_tsv('cells.tsv.gz')
        return self._cells_df.copy()

    @property
    def features(self) -> pd.DataFrame:
        """Columns of the matrix: feature_index, job_id, feature, label, and source."""
        if self._features_df is None:
            self._features_df = self._read_gz_tsv('features.tsv.gz')
        return self._features_df.copy()

    @property
    def normalization_factors(self) -> pd.DataFrame:
        """Per-cell size factors aligned by row_index, when available."""
        if self._normalization_df is None:
            self._normalization_df = self._read_gz_tsv('normalization_factors.tsv.gz')
        return self._normalization_df.copy()

    @property
    def sample_metadata(self) -> pd.DataFrame:
        """Sample metadata table keyed by sample_id, when available."""
        if self._sample_metadata_df is None:
            self._sample_metadata_df = self._read_gz_tsv('sample_metadata.tsv.gz')
        return self._sample_metadata_df.copy()

    @property
    def barcodes(self) -> pd.DataFrame:
        """Optional barcode table aligned by row_index."""
        if self._barcodes_df is None:
            self._barcodes_df = self._read_gz_tsv('barcodes.tsv.gz')
        return self._barcodes_df.copy()

    @property
    def matrix_entries(self) -> pd.DataFrame:
        """
        Sparse matrix entries as row_index, feature_index, value.

        Values are raw per-cell expression or k-mer hit counts. Missing
        row/feature pairs are zero.
        """
        if self._matrix_entries_df is not None:
            return self._matrix_entries_df.copy()

        if self.archive_path is None:
            self._matrix_entries_df = pd.DataFrame(columns=['row_index', 'feature_index', 'value'])
            return self._matrix_entries_df.copy()

        if not self.archive_path.exists():
            raise FileNotFoundError(f"Archive not found: {self.archive_path}")

        with zipfile.ZipFile(self.archive_path) as zf:
            if 'matrix.mtx.gz' not in zf.namelist():
                self._matrix_entries_df = pd.DataFrame(columns=['row_index', 'feature_index', 'value'])
                return self._matrix_entries_df.copy()

            with zf.open('matrix.mtx.gz') as raw:
                with gzip.GzipFile(fileobj=raw, mode='rb') as gz:
                    text = io.TextIOWrapper(gz, encoding='utf-8')
                    for line in text:
                        if not line.startswith('%'):
                            break
                    try:
                        entries = pd.read_csv(
                            text,
                            sep=r'\s+',
                            names=['row_index', 'feature_index', 'value'],
                            engine='python',
                        )
                    except pd.errors.EmptyDataError:
                        entries = pd.DataFrame(columns=['row_index', 'feature_index', 'value'])

        if entries.empty:
            entries = pd.DataFrame(columns=['row_index', 'feature_index', 'value'])
        else:
            entries['row_index'] = entries['row_index'].astype('int64')
            entries['feature_index'] = entries['feature_index'].astype('int64')
            entries['value'] = entries['value'].astype('float32')

        self._matrix_entries_df = entries
        return self._matrix_entries_df.copy()

    def get_cell_ids(self, sample_ids: Optional[Union[int, List[int]]] = None) -> pd.DataFrame:
        """
        Return sample_id/cell_id pairs from the matrix rows.

        Args:
            sample_ids: Optional encoded sample ID or list of sample IDs.
        """
        cells = self.cells
        if sample_ids is None or cells.empty:
            return cells
        if isinstance(sample_ids, (int, np.integer)):
            sample_set = {int(sample_ids)}
        else:
            sample_set = {int(s) for s in sample_ids}
        return cells[cells['sample_id'].isin(sample_set)].copy()

    def positive_cells(self, feature: Optional[Union[str, int]] = None,
                       sample_ids: Optional[Union[int, List[int]]] = None) -> pd.DataFrame:
        """
        Return cells with non-zero expression in any retrieved feature or one feature.

        Args:
            feature: Feature label/name or feature_index. If omitted, returns
                all cells that are present in the retrieved matrix rows.
            sample_ids: Optional encoded sample ID or list of sample IDs.
        """
        cells = self.get_cell_ids(sample_ids)
        if feature is None or cells.empty:
            return cells

        features = self.features
        if isinstance(feature, (int, np.integer)):
            feature_indices = {int(feature)}
        else:
            feature_text = str(feature)
            mask = (
                features.get('feature', pd.Series(dtype=object)).astype(str).eq(feature_text)
                | features.get('label', pd.Series(dtype=object)).astype(str).eq(feature_text)
            )
            feature_indices = set(features.loc[mask, 'feature_index'].astype(int).tolist())
        if not feature_indices:
            raise ValueError(f"Feature not found in result: {feature}")

        entries = self.matrix_entries
        row_ids = set(entries.loc[entries['feature_index'].isin(feature_indices), 'row_index'].astype(int).tolist())
        return cells[cells['row_index'].isin(row_ids)].copy()

    def to_dataframe(self, normalized: bool = False,
                     include_sample_metadata: bool = False) -> pd.DataFrame:
        """
        Convert the sparse matrix to a long DataFrame.

        Args:
            normalized: Add normalized_value = value / size_factor when
                normalization factors are available.
            include_sample_metadata: Merge sample metadata by sample_id.
        """
        entries = self.matrix_entries
        cells = self.cells
        features = self.features

        if entries.empty:
            columns = ['row_index', 'sample_id', 'cell_id', 'feature_index', 'feature', 'label', 'value']
            return pd.DataFrame(columns=columns)

        df = entries.merge(cells, on='row_index', how='left')
        df = df.merge(features, on='feature_index', how='left')

        if normalized:
            norm = self.normalization_factors
            if not norm.empty and 'size_factor' in norm.columns:
                norm = norm[['row_index', 'total_counts', 'size_factor']].copy()
                df = df.merge(norm, on='row_index', how='left')
                sf = pd.to_numeric(df['size_factor'], errors='coerce')
                df['normalized_value'] = df['value'] / sf.replace(0, np.nan)

        if include_sample_metadata:
            meta = self.sample_metadata
            if not meta.empty and 'sample_id' in meta.columns:
                df = df.merge(meta, on='sample_id', how='left', suffixes=('', '_sample'))

        return df

    def for_sample(self, sample_id: int, normalized: bool = False,
                   include_sample_metadata: bool = False) -> pd.DataFrame:
        """Return long matrix entries for one encoded sample ID."""
        cells = self.get_cell_ids(int(sample_id))
        if cells.empty:
            return pd.DataFrame()
        row_ids = set(cells['row_index'].astype(int).tolist())
        entries = self.matrix_entries
        subset = entries[entries['row_index'].isin(row_ids)]
        if subset.empty:
            return pd.DataFrame()

        original_entries = self._matrix_entries_df
        try:
            self._matrix_entries_df = subset
            return self.to_dataframe(
                normalized=normalized,
                include_sample_metadata=include_sample_metadata,
            )
        finally:
            self._matrix_entries_df = original_entries

    def to_single_cell_result(self, feature: Optional[Union[str, int]] = None,
                              sample_ids: Optional[Union[int, List[int]]] = None,
                              normalized: bool = False) -> 'SingleCellResult':
        """
        Convert one retrieved feature to the legacy SingleCellResult shape.

        This is useful for existing downstream code that expects columns
        cell_id, expression, and sample_id.
        """
        features = self.features
        if features.empty:
            return SingleCellResult({'status': 'completed', 'results': {}}, self.client)

        if feature is None:
            if len(features) != 1:
                raise ValueError("feature is required when the result has multiple features")
            feature_index = int(features.iloc[0]['feature_index'])
            feature_label = str(features.iloc[0].get('label') or features.iloc[0].get('feature'))
        elif isinstance(feature, (int, np.integer)):
            feature_index = int(feature)
            row = features[features['feature_index'] == feature_index]
            feature_label = str(row.iloc[0].get('label') or row.iloc[0].get('feature')) if not row.empty else str(feature)
        else:
            feature_text = str(feature)
            mask = (
                features.get('feature', pd.Series(dtype=object)).astype(str).eq(feature_text)
                | features.get('label', pd.Series(dtype=object)).astype(str).eq(feature_text)
            )
            rows = features.loc[mask]
            if rows.empty:
                raise ValueError(f"Feature not found in result: {feature}")
            feature_index = int(rows.iloc[0]['feature_index'])
            feature_label = str(rows.iloc[0].get('label') or rows.iloc[0].get('feature'))

        df = self.to_dataframe(normalized=normalized)
        df = df[df['feature_index'] == feature_index]
        if sample_ids is not None:
            if isinstance(sample_ids, (int, np.integer)):
                sample_set = {int(sample_ids)}
            else:
                sample_set = {int(s) for s in sample_ids}
            df = df[df['sample_id'].isin(sample_set)]

        expr_col = 'normalized_value' if normalized and 'normalized_value' in df.columns else 'value'
        result_data = {
            'status': 'completed',
            'results': {
                feature_label: {
                    'cell': df['cell_id'].astype(int).tolist(),
                    'sample': df['sample_id'].astype(int).tolist(),
                    'expression': df[expr_col].astype(float).tolist(),
                    'ncells': int(len(df)),
                }
            }
        }
        return SingleCellResult(result_data, self.client)

    def project(self, dataset_id: str,
                sample_ids: Optional[Union[int, List[int]]] = None,
                feature: Optional[Union[str, int]] = None,
                **kwargs) -> 'CoexpressionResult':
        """
        Project retrieved positive cells onto a coexpression index.

        Args:
            dataset_id: Coexpression index or dataset identifier.
            sample_ids: Optional encoded sample ID or IDs to restrict cells.
            feature: Optional feature to restrict to cells positive for that
                feature. If omitted, uses all retrieved positive cells.
            **kwargs: Additional coexpression parameters, such as top_n_genes.
        """
        if self.client is None:
            raise MalvaAPIError("A MalvaClient is required for projection")
        cells = self.positive_cells(feature=feature, sample_ids=sample_ids)
        return self.client.project_cells(
            cells['cell_id'].astype(int).tolist(),
            cells['sample_id'].astype(int).tolist(),
            dataset_id=dataset_id,
            **kwargs,
        )

    def __repr__(self) -> str:
        n_cells = self.export_record.get('cell_count')
        n_features = self.export_record.get('feature_count')
        if n_cells is None:
            try:
                n_cells = len(self.cells)
            except Exception:
                n_cells = 'unknown'
        if n_features is None:
            try:
                n_features = len(self.features)
            except Exception:
                n_features = 'unknown'
        if self.archive_path is not None:
            location = f"archive_path='{self.archive_path}'"
        else:
            location = f"job_id='{self.job_id}', source='{self.source}'"
        return f"CellExpressionMatrixResult(cells={n_cells}, features={n_features}, {location})"


class SingleCellResult:
    """
    Represents search results at the single cell level (not aggregated by cell type)
    """
    
    def __init__(self, results_data: Dict[str, Any], client: Optional['MalvaClient'] = None):
        """
        Initialize SingleCellResult
        
        Args:
            results_data: Raw results from the API
            client: MalvaClient instance for metadata enrichment
        """
        self.raw_data = results_data
        self.client = client
        self._processed_data = None  # <-- This starts as None
        
        # Extract basic info
        self.job_id = results_data.get('job_id')
        self.status = results_data.get('status', 'unknown')
        
        # Process results if available and completed
        if self.status == 'completed' and 'results' in results_data:
            self._process_results()
    
    def _process_results(self):
        """Process raw results into structured data"""
        # Check if we have a nested results structure
        if 'results' in self.raw_data and 'results' in self.raw_data['results']:
            # We have nested results - this is the actual data
            results = self.raw_data['results']['results']
            print(f"DEBUG: Found nested results structure")
        elif 'results' in self.raw_data:
            # Try the direct results first
            results = self.raw_data.get('results', {})
            print(f"DEBUG: Using direct results structure")
        else:
            results = {}
            print(f"DEBUG: No results found")
        
        # Now look for seq_1 or similar keys
        result_keys = [k for k in results.keys() if k.startswith('seq') and not k.startswith('_')]
        
        if not result_keys:
            # Try to find any dict with 'cell' key
            for key, value in results.items():
                if isinstance(value, dict) and 'cell' in value:
                    result_keys = [key]
                    break
        
        if not result_keys:
            print(f"DEBUG: No valid result keys found in: {list(results.keys())}")
            self._processed_data = {
                'cells': [],
                'expression': [],
                'samples': [],
                'query_gene': None,
                'metadata': {'error': 'No sequence data found'}
            }
            return
        
        first_key = result_keys[0]
        result_data = results[first_key]
        
        print(f"DEBUG: Processing key '{first_key}' with type {type(result_data)}")
        
        if not isinstance(result_data, dict):
            print(f"DEBUG: Result data is not a dict: {type(result_data)}")
            self._processed_data = {
                'cells': [],
                'expression': [],
                'samples': [],
                'query_gene': None,
                'metadata': {'error': f'Invalid result data type: {type(result_data)}'}
            }
            return
        
        # Extract the arrays
        cells = result_data.get('cell', [])
        expressions = result_data.get('expression', [])
        samples = result_data.get('sample', [])
        
        print(f"DEBUG: Found {len(cells)} cells, {len(expressions)} expressions, {len(samples)} samples")
        
        # Store processed data
        self._processed_data = {
            'cells': cells,
            'expression': expressions,
            'samples': samples,
            'query_gene': result_data.get('sequence', first_key),
            'metadata': {
                'total_cells': result_data.get('ncells', len(cells)),
                'unique_samples': len(set(samples)) if samples else 0,
                'sequence_length': result_data.get('sequence_length'),
                'query_info': self.raw_data.get('query', {}),
                'searches_remaining': self.raw_data.get('searches_remaining')
            }
        }
        
        print(f"✓ Processed {len(cells)} cells from {len(set(samples))} samples")

    @property
    def is_completed(self) -> bool:
        """Check if the search is completed"""
        return self.status == 'completed'
    
    @property
    def is_pending(self) -> bool:
        """Check if the search is still pending"""
        return self.status == 'pending'
    
    @property
    def has_results(self) -> bool:
        """Check if results are available"""
        return self._processed_data is not None and len(self._processed_data['cells']) > 0
    
    @property
    def cell_count(self) -> int:
        """Get total number of cells in results"""
        if not self._processed_data:
            return 0
        return len(self._processed_data['cells'])
    
    @property
    def sample_count(self) -> int:
        """Get number of unique samples in results"""
        if not self._processed_data:
            return 0
        return len(set(self._processed_data['samples']))
    
    @property
    def query_gene(self) -> Optional[str]:
        """Get the queried gene symbol or sequence"""
        if not self._processed_data:
            return None
        return self._processed_data.get('query_gene')
    
    def get_cell_ids(self) -> List[int]:
        """Get list of all cell IDs"""
        if not self.has_results:
            return []
        return self._processed_data['cells']
    
    def get_expression_values(self) -> List[float]:
        """Get list of all expression values"""
        if not self.has_results:
            return []
        return self._processed_data['expression']
    
    def get_sample_ids(self) -> List[int]:
        """Get list of all sample IDs"""
        if not self.has_results:
            return []
        return self._processed_data['samples']
    
    def to_dataframe(self) -> pd.DataFrame:
        """
        Convert results to a pandas DataFrame
        
        Returns:
            DataFrame with columns: cell_id, expression, sample_id
        """
        if not self.has_results:
            return pd.DataFrame(columns=['cell_id', 'expression', 'sample_id'])
        
        return pd.DataFrame({
            'cell_id': self._processed_data['cells'],
            'expression': self._processed_data['expression'],
            'sample_id': self._processed_data['samples']
        })
    
    def filter_by_expression(self, min_expression: float = 0, max_expression: float = float('inf')) -> 'SingleCellResult':
        """
        Filter results by expression thresholds
        
        Args:
            min_expression: Minimum expression value
            max_expression: Maximum expression value
            
        Returns:
            New SingleCellResult with filtered data
        """
        if not self.has_results:
            return self
        
        df = self.to_dataframe()
        filtered_df = df[(df['expression'] >= min_expression) & (df['expression'] <= max_expression)]
        
        # Create new result object with filtered data
        filtered_data = {
            'status': 'completed',
            'results': {
                'seq_1': {
                    'cell': filtered_df['cell_id'].tolist(),
                    'expression': filtered_df['expression'].tolist(),
                    'sample': filtered_df['sample_id'].tolist(),
                    'ncells': len(filtered_df)
                }
            }
        }
        
        return SingleCellResult(filtered_data, self.client)
    
    def filter_by_samples(self, sample_ids: List[int]) -> 'SingleCellResult':
        """
        Filter results to specific samples
        
        Args:
            sample_ids: List of sample IDs to keep
            
        Returns:
            New SingleCellResult with filtered data
        """
        if not self.has_results:
            return self
        
        df = self.to_dataframe()
        filtered_df = df[df['sample_id'].isin(sample_ids)]
        
        # Create new result object with filtered data
        filtered_data = {
            'status': 'completed',
            'results': {
                'seq_1': {
                    'cell': filtered_df['cell_id'].tolist(),
                    'expression': filtered_df['expression'].tolist(),
                    'sample': filtered_df['sample_id'].tolist(),
                    'ncells': len(filtered_df)
                }
            }
        }
        
        return SingleCellResult(filtered_data, self.client)
    
    def get_top_expressing_cells(self, n: int = 100) -> 'SingleCellResult':
        """
        Get top N expressing cells
        
        Args:
            n: Number of top cells to return
            
        Returns:
            New SingleCellResult with top expressing cells
        """
        if not self.has_results:
            return self
        
        df = self.to_dataframe()
        top_df = df.nlargest(n, 'expression')
        
        # Create new result object with top cells
        filtered_data = {
            'status': 'completed',
            'results': {
                'seq_1': {
                    'cell': top_df['cell_id'].tolist(),
                    'expression': top_df['expression'].tolist(),
                    'sample': top_df['sample_id'].tolist(),
                    'ncells': len(top_df)
                }
            }
        }
        
        return SingleCellResult(filtered_data, self.client)
    
    def aggregate_by_sample(self) -> pd.DataFrame:
        """
        Aggregate expression data by sample
        
        Returns:
            DataFrame with sample-level statistics
        """
        if not self.has_results:
            return pd.DataFrame(columns=['sample_id', 'mean_expression', 'cell_count', 'max_expression'])
        
        df = self.to_dataframe()
        
        aggregated = df.groupby('sample_id').agg({
            'expression': ['mean', 'max', 'std', 'count'],
            'cell_id': 'count'
        }).round(4)
        
        # Flatten column names
        aggregated.columns = ['mean_expression', 'max_expression', 'std_expression', 'expression_count', 'cell_count']
        aggregated = aggregated.reset_index()
        
        return aggregated
    
    def get_expression_stats(self) -> Dict[str, float]:
        """
        Get basic statistics about expression values
        
        Returns:
            Dictionary with expression statistics
        """
        if not self.has_results:
            return {}
        
        expression_values = np.array(self._processed_data['expression'])
        
        return {
            'mean': float(np.mean(expression_values)),
            'median': float(np.median(expression_values)),
            'std': float(np.std(expression_values)),
            'min': float(np.min(expression_values)),
            'max': float(np.max(expression_values)),
            'q25': float(np.percentile(expression_values, 25)),
            'q75': float(np.percentile(expression_values, 75)),
            'non_zero_cells': int(np.sum(expression_values > 0)),
            'total_cells': len(expression_values)
        }
    
    def enrich_with_metadata(self, sample_metadata: bool = True) -> pd.DataFrame:
        """
        Enrich results with metadata from the client
        
        Args:
            sample_metadata: Whether to include sample metadata
            
        Returns:
            DataFrame with enriched metadata
        """
        if not self.has_results or not self.client:
            return self.to_dataframe()
        
        df = self.to_dataframe()
        
        try:
            if sample_metadata:
                # Get unique sample IDs
                unique_samples = list(set(df['sample_id'].tolist()))
                
                # Get sample metadata
                sample_meta = self.client.get_sample_metadata(unique_samples)
                
                # Create sample metadata DataFrame
                meta_records = []
                for sample_id, metadata in sample_meta.items():
                    record = {'sample_id': int(sample_id)}
                    record.update(metadata)
                    meta_records.append(record)
                
                if meta_records:
                    sample_df = pd.DataFrame(meta_records)
                    # Merge with results
                    df = df.merge(sample_df, on='sample_id', how='left')
                    print(f"✓ Enriched with metadata for {len(sample_df)} samples")
                
        except Exception as e:
            print(f"Warning: Could not enrich with metadata: {e}")
        
        return df
    
    def save_to_csv(self, filename: str, include_metadata: bool = True):
        """
        Save results to CSV file
        
        Args:
            filename: Output filename
            include_metadata: Whether to include metadata enrichment
        """
        if include_metadata and self.client:
            df = self.enrich_with_metadata()
        else:
            df = self.to_dataframe()
        
        df.to_csv(filename, index=False)
        print(f"✓ Results saved to {filename} ({len(df)} cells)")
    
    def __repr__(self) -> str:
        if not self.has_results:
            if self.status == 'pending':
                return f"SingleCellResult(status='pending', job_id='{self.job_id}')"
            return "SingleCellResult(status='no results', cells=0)"
        
        return (f"SingleCellResult(sequence_length={self._processed_data['metadata'].get('sequence_length', 'N/A')}, "
                f"cells={self.cell_count:,}, samples={self.sample_count})")
    
    def __len__(self) -> int:
        return self.cell_count


class CoverageResult:
    """
    Represents genomic coverage data from the Malva genome browser.

    Coverage data is organized as a matrix with genomic positions as rows
    and cell types as columns. Each cell contains a coverage value.
    """

    def __init__(self, raw_data: Dict[str, Any], client: Optional['MalvaClient'] = None):
        """
        Initialize CoverageResult

        Args:
            raw_data: Raw coverage data from the API
            client: MalvaClient instance for follow-up requests
        """
        self.raw_data = raw_data
        self.client = client

        self.job_id = raw_data.get('job_id', '')
        self.region = raw_data.get('region', '')
        self.positions = raw_data.get('positions', [])
        self.cell_types = raw_data.get('cell_types', [])
        self.coverage_matrix = raw_data.get('coverage_matrix', [])
        self.total_windows = raw_data.get('total_windows', len(self.positions))

    def to_dataframe(self) -> pd.DataFrame:
        """
        Convert coverage data to a pandas DataFrame.

        Returns:
            DataFrame with positions as index and cell types as columns
        """
        if not self.coverage_matrix or not self.cell_types:
            return pd.DataFrame()

        # Server returns coverage_matrix as (cell_types x positions);
        # transpose to (positions x cell_types) for the DataFrame
        matrix = np.array(self.coverage_matrix)
        if matrix.ndim == 2 and matrix.shape[1] != len(self.cell_types):
            matrix = matrix.T

        df = pd.DataFrame(
            matrix,
            columns=self.cell_types
        )

        if self.positions:
            df.index = self.positions
            df.index.name = 'position'

        return df

    def get_filter_options(self) -> Dict[str, Any]:
        """
        Get available filter options for this coverage result.

        Returns:
            Dictionary with filter options
        """
        if not self.client or not self.job_id:
            return {}
        return self.client.get_coverage_filter_options(self.job_id)

    def download_wig(self, output_path: str, **filters) -> str:
        """
        Download coverage data as a WIG file.

        Args:
            output_path: Path to save the WIG file
            **filters: Optional metadata filters

        Returns:
            Path to the saved file
        """
        if not self.client or not self.job_id:
            raise MalvaAPIError("Client and job_id required to download WIG files")
        return self.client.download_coverage_wig(self.job_id, output_path, **filters)

    def plot(self, cell_types: Optional[List[str]] = None, **kwargs):
        """
        Plot coverage across the genomic region.

        Args:
            cell_types: Specific cell types to plot (default: all)
            **kwargs: Additional arguments passed to matplotlib
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("matplotlib required for plotting. Install with: pip install matplotlib")

        df = self.to_dataframe()
        if df.empty:
            print("No coverage data to plot")
            return

        if cell_types:
            available = [ct for ct in cell_types if ct in df.columns]
            if not available:
                raise ValueError(f"None of {cell_types} found. Available: {list(df.columns)}")
            df = df[available]

        fig, ax = plt.subplots(figsize=kwargs.get('figsize', (14, 6)))
        for col in df.columns:
            ax.plot(range(len(df)), df[col].values, label=col, alpha=0.8)

        ax.set_xlabel('Genomic Position')
        ax.set_ylabel('Coverage')
        ax.set_title(f'Coverage: {self.region}' if self.region else 'Genomic Coverage')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
        plt.tight_layout()
        return fig

    def __repr__(self) -> str:
        n_positions = len(self.positions)
        n_cell_types = len(self.cell_types)
        return (f"CoverageResult(region='{self.region}', "
                f"positions={n_positions}, cell_types={n_cell_types}, "
                f"job_id='{self.job_id}')")

    def __len__(self) -> int:
        return len(self.positions)


class UMAPCoordinates:
    """
    Lightweight container for UMAP coordinates from the coexpression API.

    Wraps the compact parallel-array format returned by
    ``GET /api/coexpression/umap/<dataset_id>`` and provides conversion to
    a pandas DataFrame and a simple scatter-plot method.
    """

    def __init__(self, raw_data: Dict[str, Any], client: Optional['MalvaClient'] = None):
        """
        Initialize UMAPCoordinates.

        Args:
            raw_data: Raw response from the UMAP endpoint
            client: MalvaClient instance for follow-up requests
        """
        self.raw_data = raw_data
        self.client = client

        self.x = raw_data.get('x', [])
        self.y = raw_data.get('y', [])
        self.metacell_ids = raw_data.get('metacell_id', [])
        self.n_cells = raw_data.get('n_cells', [])
        self.samples = raw_data.get('sample', [])
        self.clusters = raw_data.get('cluster', [])

    def to_dataframe(self) -> pd.DataFrame:
        """
        Convert to a pandas DataFrame.

        Returns:
            DataFrame with columns: x, y, metacell_id, n_cells, sample, cluster
        """
        data: Dict[str, Any] = {'x': self.x, 'y': self.y}

        if self.metacell_ids:
            data['metacell_id'] = self.metacell_ids
        if self.n_cells:
            data['n_cells'] = self.n_cells
        if self.samples:
            data['sample'] = self.samples
        if self.clusters:
            data['cluster'] = self.clusters

        return pd.DataFrame(data)

    def plot(self, color_by: str = 'cluster', point_size: Optional[float] = None,
             cmap: str = 'tab20', figsize: Tuple[int, int] = (10, 8)):
        """
        Scatter plot of UMAP coordinates.

        Args:
            color_by: Column to color points by (default ``'cluster'``)
            point_size: Marker size (auto-scaled if ``None``)
            cmap: Matplotlib colormap name
            figsize: Figure size as ``(width, height)``

        Returns:
            matplotlib Figure
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("matplotlib required for plotting. Install with: pip install matplotlib")

        df = self.to_dataframe()
        if df.empty:
            print("No UMAP data to plot")
            return

        fig, ax = plt.subplots(figsize=figsize)

        if point_size is None:
            point_size = max(1, 200 / (len(df) ** 0.5))

        if color_by in df.columns:
            col = df[color_by]
            if col.dtype == 'object' or hasattr(col, 'cat'):
                categories = col.unique()
                color_map = {cat: i for i, cat in enumerate(categories)}
                colors = col.map(color_map)
                scatter = ax.scatter(df['x'], df['y'], c=colors, s=point_size,
                                     cmap=cmap, alpha=0.7)
                handles = [plt.Line2D([0], [0], marker='o', color='w',
                                       markerfacecolor=plt.get_cmap(cmap)(color_map[cat] / max(len(categories) - 1, 1)),
                                       markersize=8, label=str(cat))
                           for cat in categories]
                ax.legend(handles=handles, bbox_to_anchor=(1.05, 1),
                          loc='upper left', fontsize='small')
            else:
                scatter = ax.scatter(df['x'], df['y'], c=col, s=point_size,
                                     cmap=cmap, alpha=0.7)
                plt.colorbar(scatter, ax=ax, label=color_by)
        else:
            ax.scatter(df['x'], df['y'], s=point_size, alpha=0.7)

        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')
        ax.set_title(f'UMAP ({color_by})')
        plt.tight_layout()
        return fig

    def __repr__(self) -> str:
        return f"UMAPCoordinates(n_points={len(self.x)}, clusters={len(set(self.clusters)) if self.clusters else 'N/A'})"

    def __len__(self) -> int:
        return len(self.x)


class CoexpressionResult:
    """
    Full coexpression analysis result from the Malva coexpression API.

    Wraps the response from ``POST /api/coexpression/query-by-job`` and
    provides DataFrame conversions, top-gene retrieval, and plotting
    helpers for correlated genes, GO enrichment, and UMAP scores.
    """

    def __init__(self, raw_data: Dict[str, Any], client: Optional['MalvaClient'] = None):
        """
        Initialize CoexpressionResult.

        Args:
            raw_data: Raw response from the coexpression endpoint
            client: MalvaClient instance for follow-up requests
        """
        self.raw_data = raw_data
        self.client = client

        self.dataset_id = raw_data.get('dataset_id', '')
        self.correlated_genes = raw_data.get('correlated_genes', [])
        self.umap_scores = raw_data.get('umap_scores', {})
        self.go_enrichment = raw_data.get('go_enrichment', [])
        self.cell_type_enrichment = raw_data.get('cell_type_enrichment', [])
        self.tissue_breakdown = raw_data.get('tissue_breakdown', [])
        self.n_query_cells = raw_data.get('n_query_cells', 0)
        self.n_mapped_metacells = raw_data.get('n_mapped_metacells', 0)

    # -- DataFrame conversions ------------------------------------------------

    def genes_to_dataframe(self) -> pd.DataFrame:
        """
        Convert correlated genes to a DataFrame.

        Returns:
            DataFrame with columns: gene, correlation, p_value (plus any
            extra fields returned by the server)
        """
        if not self.correlated_genes:
            return pd.DataFrame(columns=['gene', 'correlation', 'p_value'])
        return pd.DataFrame(self.correlated_genes)

    def scores_to_dataframe(self) -> pd.DataFrame:
        """
        Convert UMAP scores to a DataFrame.

        Returns:
            DataFrame with metacell-level score data
        """
        scores = self.umap_scores
        if not scores:
            return pd.DataFrame()

        # Handle parallel-array (compact) format
        if isinstance(scores, dict) and 'metacell_ids' in scores:
            data: Dict[str, Any] = {}
            for key, values in scores.items():
                if isinstance(values, list):
                    data[key] = values
            return pd.DataFrame(data)

        # Handle list-of-dicts format
        if isinstance(scores, list):
            return pd.DataFrame(scores)

        return pd.DataFrame()

    def umap_to_dataframe(self) -> pd.DataFrame:
        """
        Convert UMAP score data to a DataFrame with x/y coordinates.

        Falls back to :meth:`scores_to_dataframe` when coordinates are
        embedded in the scores payload.

        Returns:
            DataFrame with UMAP coordinates and scores
        """
        return self.scores_to_dataframe()

    def go_to_dataframe(self) -> pd.DataFrame:
        """
        Convert GO enrichment results to a DataFrame.

        Returns:
            DataFrame with columns such as go_id, name, fdr, etc.
        """
        if not self.go_enrichment:
            return pd.DataFrame(columns=['go_id', 'name', 'fdr'])
        return pd.DataFrame(self.go_enrichment)

    def cell_type_enrichment_to_dataframe(self) -> pd.DataFrame:
        """
        Convert cell-type enrichment to a DataFrame.

        Returns:
            DataFrame with cell-type enrichment data
        """
        if not self.cell_type_enrichment:
            return pd.DataFrame()
        return pd.DataFrame(self.cell_type_enrichment)

    def tissue_breakdown_to_dataframe(self) -> pd.DataFrame:
        """
        Convert tissue breakdown to a DataFrame.

        Returns:
            DataFrame with tissue breakdown data
        """
        if not self.tissue_breakdown:
            return pd.DataFrame()
        return pd.DataFrame(self.tissue_breakdown)

    # -- Convenience ----------------------------------------------------------

    def get_top_genes(self, n: int = 20, min_correlation: float = 0.0) -> List[str]:
        """
        Get the top *n* correlated gene names.

        Args:
            n: Number of genes to return
            min_correlation: Minimum correlation to include

        Returns:
            List of gene names
        """
        df = self.genes_to_dataframe()
        if df.empty:
            return []

        if 'correlation' in df.columns:
            df = df[df['correlation'] >= min_correlation]
            df = df.sort_values('correlation', ascending=False)

        gene_col = 'gene' if 'gene' in df.columns else df.columns[0]
        return df[gene_col].head(n).tolist()

    # -- Plotting -------------------------------------------------------------

    def plot_umap(self, color_by: str = 'positive_fraction', point_size: Optional[float] = None,
                  cmap: str = 'viridis', figsize: Tuple[int, int] = (10, 8)):
        """
        Scatter plot of UMAP coordinates coloured by a score column.

        Args:
            color_by: Column in the scores data to use for colouring
            point_size: Marker size (auto-scaled if ``None``)
            cmap: Matplotlib colormap name
            figsize: Figure size as ``(width, height)``

        Returns:
            matplotlib Figure
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("matplotlib required for plotting. Install with: pip install matplotlib")

        df = self.scores_to_dataframe()
        if df.empty:
            print("No UMAP score data to plot")
            return

        fig, ax = plt.subplots(figsize=figsize)

        if point_size is None:
            point_size = max(1, 200 / (len(df) ** 0.5))

        x_col = 'umap_x' if 'umap_x' in df.columns else df.columns[0]
        y_col = 'umap_y' if 'umap_y' in df.columns else df.columns[1]

        if color_by in df.columns:
            # resort if the data is numeric, to plot properly
            if pd.api.types.is_any_real_numeric_dtype(df[color_by]):
                df = df.sort_values(by=color_by, ascending=True)
            scatter = ax.scatter(df[x_col], df[y_col], c=df[color_by],
                                 s=point_size, cmap=cmap, alpha=0.7)
            plt.colorbar(scatter, ax=ax, label=color_by)
        else:
            ax.scatter(df[x_col], df[y_col], s=point_size, alpha=0.7)

        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')
        ax.set_title(f'Coexpression UMAP ({color_by})')
        plt.tight_layout()
        return fig

    def plot_top_genes(self, n: int = 20, figsize: Tuple[int, int] = (8, 6)):
        """
        Horizontal bar chart of the top correlated genes.

        Args:
            n: Number of genes to show
            figsize: Figure size as ``(width, height)``

        Returns:
            matplotlib Figure
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("matplotlib required for plotting. Install with: pip install matplotlib")

        df = self.genes_to_dataframe()
        if df.empty:
            print("No correlated gene data to plot")
            return

        if 'correlation' in df.columns:
            df = df.sort_values('correlation', ascending=False).head(n)
        else:
            df = df.head(n)

        gene_col = 'gene' if 'gene' in df.columns else df.columns[0]
        corr_col = 'correlation' if 'correlation' in df.columns else df.columns[1]

        fig, ax = plt.subplots(figsize=figsize)
        y_pos = range(len(df))
        ax.barh(y_pos, df[corr_col].values, color='steelblue')
        ax.set_yticks(y_pos)
        ax.set_yticklabels(df[gene_col].values)
        ax.invert_yaxis()
        ax.set_xlabel('Correlation')
        ax.set_title(f'Top {len(df)} Correlated Genes')
        plt.tight_layout()
        return fig

    def plot_go_enrichment(self, n: int = 15, figsize: Tuple[int, int] = (8, 6)):
        """
        Bar chart of GO enrichment results (−log10 FDR).

        Args:
            n: Number of GO terms to show
            figsize: Figure size as ``(width, height)``

        Returns:
            matplotlib Figure
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("matplotlib required for plotting. Install with: pip install matplotlib")

        df = self.go_to_dataframe()
        if df.empty:
            print("No GO enrichment data to plot")
            return

        fdr_col = 'fdr' if 'fdr' in df.columns else 'p_value' if 'p_value' in df.columns else None
        name_col = 'name' if 'name' in df.columns else df.columns[0]

        if fdr_col is None:
            print("No FDR or p-value column found in GO enrichment data")
            return

        df = df.sort_values(fdr_col, ascending=True).head(n)
        df['neg_log10_fdr'] = -np.log10(df[fdr_col].clip(lower=1e-300))

        fig, ax = plt.subplots(figsize=figsize)
        y_pos = range(len(df))
        ax.barh(y_pos, df['neg_log10_fdr'].values, color='darkorange')
        ax.set_yticks(y_pos)
        ax.set_yticklabels(df[name_col].values)
        ax.invert_yaxis()
        ax.set_xlabel('-log10(FDR)')
        ax.set_title(f'Top {len(df)} GO Enrichment Terms')
        plt.tight_layout()
        return fig

    def __repr__(self) -> str:
        parts = [f"CoexpressionResult(dataset='{self.dataset_id}'"]
        parts.append(f"genes={len(self.correlated_genes)}")
        parts.append(f"query_cells={self.n_query_cells}")
        parts.append(f"metacells={self.n_mapped_metacells}")
        if self.go_enrichment:
            parts.append(f"go_terms={len(self.go_enrichment)}")
        return ', '.join(parts) + ')'

    def __len__(self) -> int:
        return len(self.correlated_genes)
