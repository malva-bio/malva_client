import json
import time
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
    
    def aggregate_by(self, group_by: Union[str, List[str]], 
                    agg_func: str = 'mean',
                    expr_column: str = 'norm_expr') -> pd.DataFrame:
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
        if expr_column not in self._df.columns:
            raise ValueError(f"Expression column '{expr_column}' not found")
        
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
        
        if 'norm_expr' not in self._df.columns:
            raise ValueError("Expression column 'norm_expr' not found")
        
        # Prepare data for plotting
        plot_data = self._df.copy()
        
        # Calculate sorting metrics for each category
        if sort_by != 'alphabetical':
            if sort_by == 'mean':
                sort_values = plot_data.groupby(group_by)['norm_expr'].mean()
            elif sort_by == 'median':
                sort_values = plot_data.groupby(group_by)['norm_expr'].median()
            elif sort_by == 'count':
                sort_values = plot_data.groupby(group_by)['norm_expr'].count()
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
        sns.boxplot(data=plot_data, x=group_by, y='norm_expr', 
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
                cat_data = plot_data[plot_data[group_by] == category]['norm_expr']
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
        """
        Create comprehensive expression plots with three visualizations
        
        Args:
            group_by: Column to group by for plotting
            limit: Maximum number of categories to show
            plot_type: Type of plot for main plot ('box', 'violin', 'bar', 'strip')
            **kwargs: Additional arguments
        """
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
            import numpy as np
        except ImportError:
            raise ImportError("matplotlib and seaborn required for plotting")
        
        if group_by not in self._df.columns:
            raise ValueError(f"Column '{group_by}' not found")
        
        # Get top categories by mean expression
        top_categories = (self._df.groupby(group_by)['norm_expr']
                        .mean().nlargest(limit).index.tolist())
        
        plot_data = self._df[self._df[group_by].isin(top_categories)]
        
        # Create 1x3 subplot layout
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        fig.suptitle(f'Expression Analysis by {group_by.replace("_", " ").title()}', 
                    fontsize=16)
        
        # 1. Box plot (or other plot type)
        if plot_type == 'box':
            sns.boxplot(data=plot_data, x=group_by, y='norm_expr', 
                        order=top_categories, ax=axes[0])
            axes[0].set_title('Expression Distribution (Box Plot)')
        elif plot_type == 'violin':
            sns.violinplot(data=plot_data, x=group_by, y='norm_expr', 
                        order=top_categories, ax=axes[0])
            axes[0].set_title('Expression Distribution (Violin Plot)')
        elif plot_type == 'strip':
            sns.stripplot(data=plot_data, x=group_by, y='norm_expr', 
                        order=top_categories, ax=axes[0])
            axes[0].set_title('Expression Distribution (Strip Plot)')
        else:  # bar plot
            means = plot_data.groupby(group_by)['norm_expr'].agg(['mean', 'std'])
            means = means.reindex(top_categories)
            means.plot(kind='bar', y='mean', yerr='std', ax=axes[0], 
                    color='skyblue', capsize=4)
            axes[0].set_title('Mean Expression ± SD')
            axes[0].legend().remove()
        
        axes[0].tick_params(axis='x', rotation=45)
        axes[0].set_xlabel(group_by.replace('_', ' ').title())
        axes[0].set_ylabel('Normalized Expression')
        
        # 2. Number of samples per category
        sample_counts = plot_data.groupby(group_by)['sample_id'].nunique().reindex(top_categories)
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
        # Check if cell_count column exists, otherwise use observation count
        if 'cell_count' in plot_data.columns:
            cell_counts = plot_data.groupby(group_by)['cell_count'].sum().reindex(top_categories)
            ylabel = 'Total Cells'
            title = 'Number of Cells'
        else:
            # Fallback: count observations (rows) as proxy for cells
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
        
        # Print summary statistics
        print(f"\n📊 Summary for {group_by.replace('_', ' ').title()}:")
        print("-" * 50)
        for category in top_categories:
            cat_data = plot_data[plot_data[group_by] == category]
            n_samples = cat_data['sample_id'].nunique()
            n_cells = cat_data['cell_count'].sum() if 'cell_count' in cat_data.columns else len(cat_data)
            mean_expr = cat_data['norm_expr'].mean()
            print(f"{category}: {n_samples} samples, {n_cells:,} cells, μ={mean_expr:.3f}")
        
        return fig

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
            'expression': [col for col in all_columns if col in ['norm_expr', 'cell_count', 'sample_idx', 'cell_type_idx']],
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
        if field_name in ['norm_expr', 'cell_count', 'sample_idx', 'cell_type_idx']:
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
    """Container for search results that extends MalvaDataFrame for direct analysis"""
    
    def __init__(self, raw_data: Dict[str, Any], client: 'MalvaClient'):
        self.raw_data = raw_data
        self.client = client
        self._enriched = False
        
        # Convert to DataFrame immediately
        expr_df = self._convert_to_dataframe()
        
        # Initialize as MalvaDataFrame (without metadata initially)
        super().__init__(expr_df, client, sample_metadata={})
    
    def _convert_to_dataframe(self) -> pd.DataFrame:
        """Convert raw search results to DataFrame"""
        if not self.raw_data or not isinstance(self.raw_data, dict):
            return pd.DataFrame()
        
        results = self.raw_data.get('results', {})
        if not results:
            return pd.DataFrame()
        
        all_data = []
        for gene_seq, result in results.items():
            if gene_seq == "_sample_metadata":
                continue
                
            expression_data = result.get('expression_data', {})
            if not expression_data:
                continue
                
            # Convert compact format to DataFrame
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
            # Convert cell_type to categorical for faster filtering
            if 'cell_type' in result_df.columns:
                result_df['cell_type'] = result_df['cell_type'].astype('category')
            return result_df
        else:
            return pd.DataFrame()
    
    @property
    def job_id(self) -> str:
        """Get the job ID for this search"""
        return self.raw_data.get('job_id', '')
    
    @property
    def status(self) -> str:
        """Get the current status of the search"""
        return self.raw_data.get('status', 'unknown')
    
    @property
    def results(self) -> Dict[str, Any]:
        """Get the raw search results"""
        return self.raw_data.get('results', {})
    
    @property
    def total_cells(self) -> int:
        """Get total number of cells found"""
        return sum(result.get('ncells', 0) for result in self.results.values())
    
    def enrich_with_metadata(self) -> 'SearchResult':
        """
        Enrich the search results with sample metadata
        
        Returns:
            Self with enriched metadata
        """
        if self._enriched:
            return self
            
        if self._df.empty:
            return self
        
        # Get unique sample IDs
        sample_ids = self._df['sample_id'].unique().tolist()
        
        # Fetch sample metadata
        try:
            sample_metadata = self.client.get_sample_metadata(sample_ids)
            print(f"✓ Enriched with metadata for {len(sample_metadata)} samples")
        except Exception as e:
            logging.warning(f"Could not fetch sample metadata: {e}")
            sample_metadata = {}
        
        # Update the metadata and re-enrich
        self._sample_metadata = sample_metadata
        self._enrich_with_metadata()
        self._enriched = True
        
        return self
    
    def __repr__(self) -> str:
        """Nice representation for Jupyter notebooks"""
        if self._df.empty:
            return "SearchResult: No data found"
        
        lines = []
        lines.append("🔬 Malva Search Results")
        lines.append("=" * 50)
        
        # Basic stats
        lines.append(f"📊 Total cells: {len(self._df):,}")
        lines.append(f"🧬 Genes/sequences: {self._df['gene_sequence'].nunique() if 'gene_sequence' in self._df.columns else 'N/A'}")
        lines.append(f"🧪 Samples: {self._df['sample_id'].nunique() if 'sample_id' in self._df.columns else 'N/A'}")
        lines.append(f"🔬 Cell types: {self._df['cell_type'].nunique() if 'cell_type' in self._df.columns else 'N/A'}")
        
        # Expression stats
        if 'norm_expr' in self._df.columns:
            lines.append(f"📈 Expression range: {self._df['norm_expr'].min():.3f} - {self._df['norm_expr'].max():.3f}")
            lines.append(f"📊 Mean expression: {self._df['norm_expr'].mean():.3f}")
        
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
        
        if 'norm_expr' in self._df.columns:
            html.append(f'<tr><td><strong>Expression range:</strong></td><td>{self._df["norm_expr"].min():.3f} - {self._df["norm_expr"].max():.3f}</td></tr>')
            html.append(f'<tr><td><strong>Mean expression:</strong></td><td>{self._df["norm_expr"].mean():.3f}</td></tr>')
        
        html.append('</table>')
        
        # Metadata status
        if self._enriched:
            html.append('<div style="color: green;">✅ <strong>Enriched with sample metadata</strong></div>')
        else:
            html.append('<div style="color: orange;">ℹ️ <strong>Basic expression data only</strong></div>')
            html.append('<div style="margin: 5px 0; font-size: 0.9em;">Run <code>.enrich_with_metadata()</code> to add sample metadata for advanced filtering</div>')
        
        html.append('</div>')
        return ''.join(html)

class CoverageDataFrame(MalvaDataFrame):
    """
    Specialized MalvaDataFrame for coverage data with genomic position handling
    """
    
    def __init__(self, coverage_df: pd.DataFrame, client: 'MalvaClient', sample_metadata: Dict[int, Dict] = None):
        super().__init__(coverage_df, client, sample_metadata)
        
    def filter_by_position(self, start: int = None, end: int = None) -> 'CoverageDataFrame':
        """
        Filter coverage data by genomic position range
        
        Args:
            start: Start position (inclusive)
            end: End position (inclusive)
            
        Returns:
            New CoverageDataFrame with filtered positions
        """
        filtered_df = self._df.copy()
        
        if 'genomic_position' in filtered_df.columns:
            if start is not None:
                filtered_df = filtered_df[filtered_df['genomic_position'] >= start]
            if end is not None:
                filtered_df = filtered_df[filtered_df['genomic_position'] <= end]
        
        return CoverageDataFrame(filtered_df, self.client, self._sample_metadata)

    def filter_by(self, **kwargs) -> 'CoverageDataFrame':
        """
        Filter coverage data while preserving genomic positions
        Uses simplified field names (e.g., 'organ' instead of 'specimen_from_organism.organ')
        
        Args:
            **kwargs: Field=value pairs for filtering
            
        Returns:
            New CoverageDataFrame with filtered data (preserves genomic positions)
        """
        if not kwargs:
            return CoverageDataFrame(self._df.copy(), self.client, self._sample_metadata)
        
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
            
            # Apply filter using boolean indexing
            if isinstance(value, (list, tuple)):
                field_mask = self._df[actual_field].isin(value)
            else:
                field_mask = self._df[actual_field] == value
            
            mask = mask & field_mask
        
        # Apply the final mask in one operation
        filtered_df = self._df[mask].copy()
        
        return CoverageDataFrame(filtered_df, self.client, self._sample_metadata)
    
    def aggregate_by(self, group_by: Union[str, List[str]], 
                    agg_func: str = 'mean',
                    expr_column: str = 'coverage') -> 'CoverageDataFrame':
        """
        Aggregate coverage data by specified grouping variables while preserving genomic positions
        
        Args:
            group_by: Column name(s) to group by
            agg_func: Aggregation function ('mean', 'median', 'sum', 'count', 'std')
            expr_column: Expression column to aggregate (default: 'coverage')
            
        Returns:
            New CoverageDataFrame with aggregated results that can still be plotted
            
        Example:
            df.aggregate_by('cell_type')  # Returns aggregated coverage by position and cell_type
            df.aggregate_by(['organ', 'cell_type'])
        """
        if expr_column not in self._df.columns:
            raise ValueError(f"Expression column '{expr_column}' not found")
        
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
        
        # IMPORTANT: Always include genomic_position in grouping for coverage data
        if 'genomic_position' in self._df.columns and 'genomic_position' not in resolved_group_cols:
            final_group_cols = ['genomic_position'] + resolved_group_cols
        else:
            final_group_cols = resolved_group_cols
        
        # Simple aggregation approach
        grouped = self._df.groupby(final_group_cols)
        
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
        
        # Build result DataFrame with coverage-specific structure
        result_data = {
            'coverage': expr_agg,  # Always name it 'coverage' for plotting compatibility
            'n_observations': grouped.size(),
        }
        
        # Add additional metrics if columns exist
        if 'cell_count' in self._df.columns:
            result_data['total_cells'] = grouped['cell_count'].sum()
        
        if 'sample_id' in self._df.columns:
            result_data['n_samples'] = grouped['sample_id'].nunique()
        
        result_df = pd.DataFrame(result_data).reset_index()
        
        # Ensure genomic_position is preserved and coverage column exists for plotting
        if 'genomic_position' not in result_df.columns:
            logging.warning("genomic_position not preserved in aggregation - plotting may not work")
        
        # Don't pass sample_metadata to aggregated data
        # Convert to CoverageDataFrame to maintain plotting capabilities but without sample metadata
        return CoverageDataFrame(result_df, self.client, sample_metadata={})
    
    def available_fields(self) -> Dict[str, List[str]]:
        """
        Get categorized list of available fields for filtering/grouping (coverage-specific)
        
        Returns:
            Dictionary with categorized field names
        """
        all_columns = list(self._df.columns)
        
        # Categorize fields with coverage-specific categories
        categories = {
            'coverage': [col for col in all_columns if col in ['coverage', 'cell_type_coverage', 'genomic_position']],
            'basic': [col for col in all_columns if col in ['sample_id', 'cell_type', 'cell_count']],
            'biological': [col for col in all_columns if col in ['organ', 'organ_part', 'disease', 'species', 'development_stage', 'sex', 'age']],
            'study': [col for col in all_columns if col in ['study', 'study_title', 'laboratory', 'protocol']],
            'technical': [col for col in all_columns if col in ['chunk_id', 'internal_sample_id', 'num_cells', 'num_data', 'num_kmers', 'sample_uuid']],
            'original_fields': [col for col in all_columns if '.' in col]
        }
        
        # Remove empty categories
        return {k: v for k, v in categories.items() if v}
    
    def plot_coverage_profile(self, group_by: str = None, limit: int = 10, **kwargs):
        """
        Plot coverage profile across genomic positions
        
        Args:
            group_by: Group coverage by field ('sample_id', 'cell_type', etc.)
            limit: Maximum number of groups to plot
            **kwargs: Additional matplotlib arguments
        """
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
        except ImportError:
            raise ImportError("matplotlib and seaborn required for plotting")
        
        if 'genomic_position' not in self._df.columns or 'coverage' not in self._df.columns:
            raise ValueError("genomic_position and coverage columns required")
        
        plt.figure(figsize=kwargs.get('figsize', (15, 6)))
        
        if group_by and group_by in self._df.columns:
            # Plot by groups
            top_groups = self._df.groupby(group_by)['coverage'].mean().nlargest(limit).index
            
            for group in top_groups:
                group_data = self._df[self._df[group_by] == group]
                if not group_data.empty:
                    # Check if data is already aggregated by position
                    if group_data.groupby('genomic_position').size().max() == 1:
                        # Data is already aggregated - plot directly
                        group_data = group_data.sort_values('genomic_position')
                        plt.plot(group_data['genomic_position'], group_data['coverage'], 
                            alpha=0.7, label=f'{group_by}: {group}')
                    else:
                        # Aggregate by position to reduce noise
                        agg_data = group_data.groupby('genomic_position')['coverage'].mean().reset_index()
                        plt.plot(agg_data['genomic_position'], agg_data['coverage'], 
                            alpha=0.7, label=f'{group_by}: {group}')
        else:
            # Plot overall coverage
            if self._df.groupby('genomic_position').size().max() == 1:
                # Data is already aggregated
                plot_data = self._df.sort_values('genomic_position')
                plt.plot(plot_data['genomic_position'], plot_data['coverage'], alpha=0.8)
            else:
                # Aggregate by position first
                agg_data = self._df.groupby('genomic_position')['coverage'].mean().reset_index()
                plt.plot(agg_data['genomic_position'], agg_data['coverage'], alpha=0.8)
        
        plt.xlabel('Genomic Position')
        plt.ylabel('Coverage')
        plt.title('Coverage Profile')
        
        if group_by:
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        return plt.gcf()

class CoverageResult(CoverageDataFrame):
    """Container for genome browser coverage results that extends CoverageDataFrame for direct analysis"""
    
    def __init__(self, raw_data: Dict[str, Any]):
        self.raw_data = raw_data
        self._enriched = False
        self._client = None
        self._sample_metadata = {}
        
        # Convert to DataFrame immediately
        coverage_df = self._convert_to_dataframe()
        
        # Initialize as CoverageDataFrame (without metadata initially)
        super().__init__(coverage_df, client=None, sample_metadata={})
        
    def _convert_to_dataframe(self) -> pd.DataFrame:
        """Convert raw coverage results to DataFrame"""
        if not self.raw_data or not isinstance(self.raw_data, dict):
            return pd.DataFrame()
        
        coverage_data = self.raw_data.get('coverage_data', {})
        if not coverage_data:
            return pd.DataFrame()
        
        # Check if we have sample-level data
        if self._has_sample_data():
            return self._convert_sample_level_data()
        else:
            return self._convert_aggregated_data()
    
    def _has_sample_data(self) -> bool:
        """Check if sample-level data is available"""
        coverage_data = self.raw_data.get('coverage_data', {})
        window_info = coverage_data.get('window_info', {})
        return window_info.get('mode', 'aggregated') == 'sample_preserved'
    
    def _convert_sample_level_data(self) -> pd.DataFrame:
        """Convert sample-level coverage data to long-format DataFrame"""
        coverage_data = self.raw_data.get('coverage_data', {})
        sample_coverage_matrix = coverage_data.get('sample_coverage_matrix', [])
        positions = coverage_data.get('positions', [])
        samples = coverage_data.get('samples', [])
        cell_type_data = coverage_data.get('cell_type_data', {})
        
        if not sample_coverage_matrix or not positions:
            return pd.DataFrame()
        
        # Create long-format DataFrame efficiently
        rows = []
        
        # Pre-process cell type data for efficiency
        sample_cell_type_info = {}
        for sample_id, ct_data in cell_type_data.items():
            sample_id_int = int(sample_id)
            sample_cell_type_info[sample_id_int] = {}
            
            for cell_type, expressions in ct_data.items():
                if expressions:
                    total_expr = sum(expr['expression'] * expr['cell_count'] for expr in expressions)
                    total_cells = sum(expr['cell_count'] for expr in expressions)
                    avg_expr = total_expr / total_cells if total_cells > 0 else 0
                    sample_cell_type_info[sample_id_int][cell_type] = {
                        'cell_type_coverage': avg_expr,
                        'cell_count': total_cells
                    }
        
        # Build rows efficiently
        for i, sample_id in enumerate(samples):
            sample_coverage = sample_coverage_matrix[i] if i < len(sample_coverage_matrix) else [0] * len(positions)
            
            # Check if this sample has cell type data
            if sample_id in sample_cell_type_info and sample_cell_type_info[sample_id]:
                # Create one row per cell type
                for cell_type, ct_info in sample_cell_type_info[sample_id].items():
                    for j, position in enumerate(positions):
                        coverage = sample_coverage[j] if j < len(sample_coverage) else 0
                        
                        rows.append({
                            'sample_id': sample_id,
                            'genomic_position': position,
                            'coverage': coverage,
                            'cell_type': cell_type,
                            'cell_type_coverage': ct_info['cell_type_coverage'],
                            'cell_count': ct_info['cell_count']
                        })
            else:
                # No cell type data, create basic rows
                for j, position in enumerate(positions):
                    coverage = sample_coverage[j] if j < len(sample_coverage) else 0
                    
                    rows.append({
                        'sample_id': sample_id,
                        'genomic_position': position,
                        'coverage': coverage
                    })
        
        df = pd.DataFrame(rows)
        # Convert cell_type to categorical for faster filtering
        if 'cell_type' in df.columns and not df.empty:
            df['cell_type'] = df['cell_type'].astype('category')
        return df
    
    def _convert_aggregated_data(self) -> pd.DataFrame:
        """Convert aggregated coverage data to DataFrame"""
        coverage_data = self.raw_data.get('coverage_data', {})
        positions = coverage_data.get('positions', [])
        coverage_matrix = coverage_data.get('coverage_matrix', [])
        cell_types = coverage_data.get('cell_types', [])
        
        rows = []
        for i, cell_type in enumerate(cell_types):
            if i < len(coverage_matrix):
                for j, position in enumerate(positions):
                    if j < len(coverage_matrix[i]):
                        rows.append({
                            'cell_type': cell_type,
                            'genomic_position': position,
                            'coverage': coverage_matrix[i][j]
                        })
        
        return pd.DataFrame(rows)
    
    def set_client(self, client: 'MalvaClient'):
        """Set the client for metadata enrichment"""
        self._client = client
        self.client = client  # Update the inherited client reference
    
    @property
    def coverage_data(self) -> Dict[str, Any]:
        """Get the coverage data"""
        return self.raw_data.get('coverage_data', {})
    
    @property
    def mode(self) -> str:
        """Get the data mode (aggregated or sample_preserved)"""
        return self.coverage_data.get('window_info', {}).get('mode', 'aggregated')
    
    @property
    def has_sample_data(self) -> bool:
        """Check if sample-level data is available"""
        return self.mode == 'sample_preserved'
    
    @property
    def samples(self) -> List[int]:
        """Get sample IDs (only available if preserve_samples=True)"""
        return self.coverage_data.get('samples', [])
    
    @property
    def positions(self) -> List[int]:
        """Get genomic positions"""
        return self.coverage_data.get('positions', [])
    
    @property
    def region_info(self) -> Dict[str, Any]:
        """Get genomic region information"""
        return self.coverage_data.get('window_info', {}).get('region', {})
    
    def get_sample_coverage_dataframe(self) -> pd.DataFrame:
        """
        Get sample-level coverage data as DataFrame
        Returns: DataFrame with samples as rows, positions as columns
        """
        if not self.has_sample_data:
            raise ValueError("Sample-level data not available. Use preserve_samples=True when calling get_coverage()")
        
        sample_coverage_matrix = self.coverage_data.get('sample_coverage_matrix', [])
        positions = self.positions
        samples = self.samples
        
        if not sample_coverage_matrix or not positions:
            return pd.DataFrame()
        
        df = pd.DataFrame(sample_coverage_matrix, index=samples, columns=positions)
        df.index.name = 'sample_id'
        df.columns.name = 'genomic_position'
        return df
    
    def get_coverage_dataframe(self) -> 'CoverageDataFrame':
        """
        Get coverage data as a CoverageDataFrame for analysis
        This method is kept for backward compatibility, but now the CoverageResult itself is a CoverageDataFrame
        """
        return self
    
    def enrich_with_metadata(self) -> 'CoverageResult':
        """
        Enrich the coverage results with sample metadata
        
        Returns:
            Self with enriched metadata
        """
        if self._enriched or not self._client or not self.has_sample_data:
            return self
        
        try:
            sample_metadata = self._client.get_sample_metadata(self.samples)
            self._sample_metadata = sample_metadata
            
            # Re-enrich the DataFrame with new metadata
            self._enrich_with_metadata()
            self._enriched = True
            
            print(f"✓ Enriched with metadata for {len(sample_metadata)} samples")
        except Exception as e:
            logging.warning(f"Could not fetch sample metadata: {e}")
        
        return self
    
    # Remove the old filter_by, aggregate_by, plot_coverage_by methods since they're now inherited
    
    def plot_coverage_summary(self, group_by: str = 'cell_type', limit: int = 10, **kwargs):
        """
        Create comprehensive coverage plots
        
        Args:
            group_by: Field to group by
            limit: Maximum number of groups to show
            **kwargs: Additional arguments
        """
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
        except ImportError:
            raise ImportError("matplotlib and seaborn required for plotting")
        
        # Ensure we have enriched data for biological grouping
        biological_fields = {'organ', 'disease', 'species', 'study', 'laboratory', 'development_stage', 'sex', 'age'}
        if not self._enriched and group_by in biological_fields:
            self.enrich_with_metadata()
        
        # Create 1x3 subplot layout
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        fig.suptitle(f'Coverage Analysis by {group_by.replace("_", " ").title()}', fontsize=16)
        
        if group_by not in self._df.columns:
            available_fields = self.available_filter_fields()
            print(f"Warning: {group_by} not found. Available fields: {available_fields}")
            return fig
        
        # Get top groups by mean coverage
        top_groups = self._df.groupby(group_by)['coverage'].mean().nlargest(limit).index.tolist()
        plot_data = self._df[self._df[group_by].isin(top_groups)]
        
        # 1. Coverage distribution box plot
        sns.boxplot(data=plot_data, x=group_by, y='coverage', order=top_groups, ax=axes[0])
        axes[0].set_title('Coverage Distribution')
        axes[0].tick_params(axis='x', rotation=45)
        axes[0].set_ylabel('Coverage')
        
        # 2. Number of samples per group
        if 'sample_id' in plot_data.columns:
            sample_counts = plot_data.groupby(group_by)['sample_id'].nunique().reindex(top_groups)
            sample_counts.plot(kind='bar', ax=axes[1], color='lightcoral')
            axes[1].set_title('Number of Samples')
            axes[1].set_ylabel('Sample Count')
            axes[1].tick_params(axis='x', rotation=45)
            
            # Add value labels on bars
            for i, v in enumerate(sample_counts.values):
                axes[1].text(i, v + 0.01 * max(sample_counts.values), str(v), 
                            ha='center', va='bottom', fontsize=9)
        else:
            axes[1].text(0.5, 0.5, 'Sample data\nnot available', 
                        ha='center', va='center', transform=axes[1].transAxes)
        
        # 3. Average coverage by group
        avg_coverage = plot_data.groupby(group_by)['coverage'].mean().reindex(top_groups)
        avg_coverage.plot(kind='bar', ax=axes[2], color='lightgreen')
        axes[2].set_title('Average Coverage')
        axes[2].set_ylabel('Mean Coverage')
        axes[2].tick_params(axis='x', rotation=45)
        
        # Add value labels
        for i, v in enumerate(avg_coverage.values):
            axes[2].text(i, v + 0.01 * max(avg_coverage.values), f'{v:.3f}', 
                        ha='center', va='bottom', fontsize=9)
        
        plt.tight_layout()
        
        # Print summary
        print(f"\n📊 Coverage Summary by {group_by.replace('_', ' ').title()}:")
        print("-" * 50)
        for group in top_groups:
            group_data = plot_data[plot_data[group_by] == group]
            mean_cov = group_data['coverage'].mean()
            n_samples = group_data['sample_id'].nunique() if 'sample_id' in group_data.columns else len(group_data)
            n_positions = group_data['genomic_position'].nunique()
            print(f"{group}: {n_samples} samples, {n_positions} positions, μ={mean_cov:.3f}")
        
        return fig
    
    def __repr__(self) -> str:
        """Nice representation for Jupyter notebooks"""
        lines = []
        lines.append("🧬 Malva Coverage Results")
        lines.append("=" * 50)
        
        # Basic stats
        region = self.region_info
        chromosome = region.get('chromosome', 'Unknown')
        start = region.get('start', 0)
        end = region.get('end', 0)
        span = region.get('span', end - start)
        
        lines.append(f"📍 Region: {chromosome}:{start:,}-{end:,} ({span:,} bp)")
        lines.append(f"📊 Positions: {len(self.positions):,}")
        lines.append(f"🧪 Samples: {len(self.samples):,}")
        
        # Data info
        if not self._df.empty:
            lines.append(f"📈 Coverage range: {self._df['coverage'].min():.3f} - {self._df['coverage'].max():.3f}")
            lines.append(f"📊 Mean coverage: {self._df['coverage'].mean():.3f}")
            lines.append(f"🔬 Cell types: {self._df['cell_type'].nunique() if 'cell_type' in self._df.columns else 'N/A'}")
        
        # Mode and data type
        if self.has_sample_data:
            lines.append("✅ Sample-level data available")
            
            # Cell type info
            cell_type_data = self.coverage_data.get('cell_type_data', {})
            all_cell_types = set()
            for sample_ct_data in cell_type_data.values():
                all_cell_types.update(sample_ct_data.keys())
            
            if all_cell_types:
                lines.append(f"🔬 Cell types detected: {len(all_cell_types)}")
        else:
            lines.append("ℹ️  Aggregated data only")
        
        lines.append("")
        
        # Metadata status
        if self._enriched:
            lines.append("✅ Enriched with sample metadata")
            # Show available metadata fields
            biological_fields = [f for f in ['organ', 'disease', 'species', 'study'] if f in self._df.columns]
            if biological_fields:
                lines.append(f"🏷️  Available metadata: {', '.join(biological_fields)}")
        else:
            lines.append("ℹ️  Basic coverage data only")
            lines.append("💡 Run .enrich_with_metadata() to add sample metadata for filtering by:")
            lines.append("   • Organ, disease, species")
            lines.append("   • Study, laboratory, protocol")
            lines.append("   • Age, sex, development stage")
        
        lines.append("")
        lines.append("🔍 Available methods (inherited from CoverageDataFrame):")
        lines.append("   • .filter_by(organ='brain', cell_type='neuron')")
        lines.append("   • .aggregate_by('cell_type')")
        lines.append("   • .plot_expression_by('cell_type')  # plots coverage")
        lines.append("   • .plot_coverage_summary('organ')")
        lines.append("   • .available_fields()  # see all available fields")
        lines.append("   • .available_filter_fields()  # see filterable fields")
        lines.append("   • .filter_by_position(start=1000, end=2000)")
        
        return "\n".join(lines)
    
    def __str__(self) -> str:
        """String representation"""
        return self.__repr__()
    
    def _repr_html_(self) -> str:
        """HTML representation for Jupyter notebooks"""
        region = self.region_info
        chromosome = region.get('chromosome', 'Unknown')
        start = region.get('start', 0)
        end = region.get('end', 0)
        span = region.get('span', end - start)
        
        html = []
        html.append('<div style="border: 1px solid #ddd; padding: 10px; border-radius: 5px; font-family: monospace;">')
        html.append('<h3 style="margin-top: 0;">🧬 Malva Coverage Results</h3>')
        
        # Stats table
        html.append('<table style="border-collapse: collapse; margin: 10px 0;">')
        html.append(f'<tr><td><strong>Region:</strong></td><td>{chromosome}:{start:,}-{end:,} ({span:,} bp)</td></tr>')
        html.append(f'<tr><td><strong>Positions:</strong></td><td>{len(self.positions):,}</td></tr>')
        html.append(f'<tr><td><strong>Samples:</strong></td><td>{len(self.samples):,}</td></tr>')
        
        if not self._df.empty:
            html.append(f'<tr><td><strong>Coverage range:</strong></td><td>{self._df["coverage"].min():.3f} - {self._df["coverage"].max():.3f}</td></tr>')
            html.append(f'<tr><td><strong>Mean coverage:</strong></td><td>{self._df["coverage"].mean():.3f}</td></tr>')
            html.append(f'<tr><td><strong>Cell types:</strong></td><td>{self._df["cell_type"].nunique() if "cell_type" in self._df.columns else "N/A"}</td></tr>')
        
        if self.has_sample_data:
            html.append('<tr><td><strong>Data type:</strong></td><td>Sample-level data available</td></tr>')
        else:
            html.append('<tr><td><strong>Data type:</strong></td><td>Aggregated data only</td></tr>')
        
        html.append('</table>')
        
        # Metadata status
        if self._enriched:
            html.append('<div style="color: green;">✅ <strong>Enriched with sample metadata</strong></div>')
            biological_fields = [f for f in ['organ', 'disease', 'species', 'study'] if f in self._df.columns]
            if biological_fields:
                html.append(f'<div style="margin: 5px 0; font-size: 0.9em;">Available metadata: {", ".join(biological_fields)}</div>')
        else:
            html.append('<div style="color: orange;">ℹ️ <strong>Basic coverage data only</strong></div>')
            html.append('<div style="margin: 5px 0; font-size: 0.9em;">Run <code>.enrich_with_metadata()</code> to add sample metadata for advanced filtering</div>')
        
        html.append('<div style="margin-top: 10px; font-size: 0.9em; color: #666;">')
        html.append('<strong>Available methods:</strong> filter_by(), aggregate_by(), available_fields(), plot_coverage_summary(), and more...')
        html.append('</div>')
        
        html.append('</div>')
        return ''.join(html)


import pandas as pd
import numpy as np
from typing import Dict, List, Any, Optional, Union, TYPE_CHECKING

if TYPE_CHECKING:
    from malva_client.client import MalvaClient

class SingleCellResult:
    """
    Represents search results at the single cell level (not aggregated by cell type)
    
    This class provides access to individual cell data including:
    - Cell IDs
    - Expression values
    - Sample information
    - Metadata enrichment capabilities
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
        self._processed_data = None
        
        # Extract basic info
        if 'job_id' in results_data:
            self.job_id = results_data['job_id']
            self.status = results_data.get('status', 'unknown')
        else:
            self.job_id = None
            self.status = 'completed'
        
        # Process results if available
        if 'results' in results_data and results_data['results']:
            self._process_results()
    
    def _process_results(self):
        """Process raw results into structured data"""
        results = self.raw_data['results']
        
        # Get the first (and typically only) result key
        result_keys = list(results.keys())
        if not result_keys:
            self._processed_data = {
                'cells': [],
                'expression': [],
                'samples': [],
                'metadata': {}
            }
            return
        
        first_key = result_keys[0]
        result_data = results[first_key]
        
        self._processed_data = {
            'cells': result_data.get('cell', []),
            'expression': result_data.get('expression', []),
            'samples': result_data.get('sample', []),
            'query_gene': first_key,
            'metadata': {
                'total_cells': len(result_data.get('cell', [])),
                'unique_samples': len(set(result_data.get('sample', []))),
                'query_info': self.raw_data.get('query_info', {}),
                'search_params': self.raw_data.get('search_params', {})
            }
        }
    
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
        if not self.has_results:
            return 0
        return len(self._processed_data['cells'])
    
    @property
    def sample_count(self) -> int:
        """Get number of unique samples in results"""
        if not self.has_results:
            return 0
        return len(set(self._processed_data['samples']))
    
    @property
    def query_gene(self) -> Optional[str]:
        """Get the queried gene symbol"""
        if not self.has_results:
            return None
        return self._processed_data.get('query_gene')
    
    def get_cell_ids(self) -> List[str]:
        """Get list of all cell IDs"""
        if not self.has_results:
            return []
        return self._processed_data['cells']
    
    def get_expression_values(self) -> List[float]:
        """Get list of all expression values"""
        if not self.has_results:
            return []
        return self._processed_data['expression']
    
    def get_sample_ids(self) -> List[Union[str, int]]:
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
            'results': {
                self.query_gene: {
                    'cell': filtered_df['cell_id'].tolist(),
                    'expression': filtered_df['expression'].tolist(),
                    'sample': filtered_df['sample_id'].tolist()
                }
            }
        }
        
        return SingleCellResult(filtered_data, self.client)
    
    def filter_by_samples(self, sample_ids: List[Union[str, int]]) -> 'SingleCellResult':
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
            'results': {
                self.query_gene: {
                    'cell': filtered_df['cell_id'].tolist(),
                    'expression': filtered_df['expression'].tolist(),
                    'sample': filtered_df['sample_id'].tolist()
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
            'results': {
                self.query_gene: {
                    'cell': top_df['cell_id'].tolist(),
                    'expression': top_df['expression'].tolist(),
                    'sample': top_df['sample_id'].tolist()
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
    
    def enrich_with_metadata(self, sample_metadata: bool = True, cell_metadata: bool = False) -> pd.DataFrame:
        """
        Enrich results with metadata from the client
        
        Args:
            sample_metadata: Whether to include sample metadata
            cell_metadata: Whether to include cell metadata
            
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
                sample_df = pd.DataFrame.from_dict(sample_meta, orient='index')
                sample_df.index.name = 'sample_id'
                sample_df = sample_df.reset_index()
                
                # Merge with results
                df = df.merge(sample_df, on='sample_id', how='left')
            
            if cell_metadata:
                # Get cell metadata (this is more complex as we need sample_id, cell_id pairs)
                cell_sample_pairs = [(row['sample_id'], row['cell_id']) for _, row in df.iterrows()]
                
                # This would require the client to support cell metadata retrieval
                # Implementation depends on your specific API
                pass
                
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
        print(f"Results saved to {filename} ({len(df)} cells)")
    
    def __repr__(self) -> str:
        if not self.has_results:
            return f"SingleCellResult(status='{self.status}', cells=0)"
        
        return (f"SingleCellResult(gene='{self.query_gene}', "
                f"cells={self.cell_count}, samples={self.sample_count})")
    
    def __len__(self) -> int:
        return self.cell_count