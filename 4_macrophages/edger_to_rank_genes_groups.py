"""
Convert edgeR multi-group comparison differential expression results to scanpy rank_genes_groups format.

This module provides functions to convert edgeR output DataFrame to a format compatible with 
sc.tl.rank_genes_groups and store results in AnnData.uns['rank_genes_groups'].

edgeR output format:
- Index: gene names
- Columns: 'logCPM.<cond_i>_vs_<cond_j>', 'PValue.<cond_i>_vs_<cond_j>', 
           'FDR.<cond_i>_vs_<cond_j>', 'logFC.<cond_i>_vs_<cond_j>'
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional
import anndata as ad


def parse_comparison_columns(edger_df: pd.DataFrame) -> Dict[str, List[str]]:
    """Parse edgeR DataFrame column names to extract all comparison pairs.
    
    Parameters
    ----------
    edger_df : pd.DataFrame
        edgeR output DataFrame with gene names as index
        
    Returns
    -------
    Dict[str, List[str]]
        Dictionary mapping comparison names to their metric columns
    """
    comparisons = {}
    for col in edger_df.columns:
        if '.' in col:
            parts = col.split('.', 1)
            if len(parts) == 2:
                metric, comparison = parts
                if comparison not in comparisons:
                    comparisons[comparison] = []
                comparisons[comparison].append(metric)
    return comparisons


def parse_comparison_pair(comparison: str) -> Tuple[str, str]:
    """Parse comparison pair name to extract group name and reference group.
    
    Parameters
    ----------
    comparison : str
        Comparison name in format 'cond_i_vs_cond_j'
        
    Returns
    -------
    Tuple[str, str]
        (group_name, reference_group_name)
    """
    if '_vs_' in comparison:
        parts = comparison.split('_vs_', 1)
        if len(parts) == 2:
            return parts[0], parts[1]
    return comparison, ''


def create_rank_genes_groups_structure(
    edger_df: pd.DataFrame,
    comparisons: Optional[List[str]] = None,
    n_genes: Optional[int] = None,
    groupby: Optional[str] = None
) -> Dict:
    """Convert edgeR DataFrame to rank_genes_groups format structured arrays.
    
    Parameters
    ----------
    edger_df : pd.DataFrame
        edgeR output DataFrame
    comparisons : Optional[List[str]]
        List of comparison pairs to process. If None, extract all automatically
    n_genes : Optional[int]
        Number of genes to keep per comparison. If None, keep all genes
    groupby : Optional[str]
        Groupby parameter for params. If None, will be set to None (not required for edgeR results)
        
    Returns
    -------
    Dict
        Dictionary containing 'names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj'
    """
    comparison_dict = parse_comparison_columns(edger_df)
    
    if comparisons is None:
        comparisons = sorted(comparison_dict.keys())
    
    required_metrics = ['logCPM', 'PValue', 'FDR', 'logFC']
    valid_comparisons = []
    for comp in comparisons:
        if comp in comparison_dict:
            metrics = comparison_dict[comp]
            if all(metric in metrics for metric in required_metrics):
                valid_comparisons.append(comp)
            else:
                missing = [m for m in required_metrics if m not in metrics]
                print(f"Warning: Comparison '{comp}' missing metrics: {missing}")
        else:
            print(f"Warning: Comparison '{comp}' not found in DataFrame")
    
    if len(valid_comparisons) == 0:
        raise ValueError("No valid comparisons found")
    
    comparisons = valid_comparisons
    
    # Extract group names and reference groups
    # For edgeR pairwise comparisons, 'A_vs_B' means A relative to B
    # Use A as group name, B as reference group
    group_to_comparison = {}
    reference_groups = set()
    
    for comp in comparisons:
        group_name, ref_group = parse_comparison_pair(comp)
        if ref_group:
            reference_groups.add(ref_group)
            if group_name not in group_to_comparison:
                group_to_comparison[group_name] = comp
            else:
                print(f"Warning: Duplicate group name '{group_name}', using first comparison")
        else:
            group_to_comparison[comp] = comp
    
    # Determine reference group
    if len(reference_groups) == 1:
        reference = list(reference_groups)[0]
    elif len(reference_groups) > 1:
        reference = list(reference_groups)[0]
        print(f"Warning: Multiple reference groups {reference_groups}, using '{reference}'")
    else:
        reference = 'rest'
        print(f"Warning: Could not extract reference group, using 'rest'")
    
    group_names = list(group_to_comparison.keys())
    n_total_genes = len(edger_df.index)
    n_genes = n_genes or n_total_genes
    
    # Initialize result arrays
    names_dtype = np.dtype([(name, 'U100') for name in group_names])
    numeric_dtype = np.dtype([(name, 'f8') for name in group_names])
    
    names_array = np.empty(n_genes, dtype=names_dtype)
    scores_array = np.empty(n_genes, dtype=numeric_dtype)
    logfoldchanges_array = np.empty(n_genes, dtype=numeric_dtype)
    pvals_array = np.empty(n_genes, dtype=numeric_dtype)
    pvals_adj_array = np.empty(n_genes, dtype=numeric_dtype)
    
    # Fill arrays for each group
    for group_name in group_names:
        comp = group_to_comparison[group_name]
        logcpm_col = f'logCPM.{comp}'
        pvalue_col = f'PValue.{comp}'
        fdr_col = f'FDR.{comp}'
        logfc_col = f'logFC.{comp}'
        
        comp_data = edger_df[[logcpm_col, pvalue_col, fdr_col, logfc_col]].copy()
        comp_data['_abs_logfc'] = comp_data[logfc_col].abs()
        comp_data = comp_data.sort_values(
            by=[fdr_col, '_abs_logfc'],
            ascending=[True, False]
        ).drop(columns=['_abs_logfc'])
        
        top_genes = comp_data.head(n_genes)
        gene_names = top_genes.index.tolist()
        if len(gene_names) < n_genes:
            gene_names.extend([''] * (n_genes - len(gene_names)))
        
        names_array[group_name] = gene_names[:n_genes]
        scores_array[group_name] = top_genes[logcpm_col].fillna(0).values[:n_genes]
        logfoldchanges_array[group_name] = top_genes[logfc_col].fillna(0).values[:n_genes]
        pvals_array[group_name] = top_genes[pvalue_col].fillna(1.0).values[:n_genes]
        pvals_adj_array[group_name] = top_genes[fdr_col].fillna(1.0).values[:n_genes]
    
    return {
        'names': names_array,
        'scores': scores_array,
        'logfoldchanges': logfoldchanges_array,
        'pvals': pvals_array,
        'pvals_adj': pvals_adj_array,
        'params': {
            'groupby': groupby,
            'method': 'edger',
            'use_raw': False,
            'reference': reference,
            'n_genes': n_genes
        }
    }


def add_rank_genes_groups_to_anndata(
    adata: ad.AnnData,
    edger_df: pd.DataFrame,
    key: str = 'rank_genes_groups',
    comparisons: Optional[List[str]] = None,
    n_genes: Optional[int] = None,
    groupby: Optional[str] = None
) -> ad.AnnData:
    """Add edgeR results to AnnData.uns in rank_genes_groups format.
    
    Parameters
    ----------
    adata : ad.AnnData
        AnnData object to add results to
    edger_df : pd.DataFrame
        edgeR output DataFrame with gene names as index
    key : str
        Key name in uns, default 'rank_genes_groups'
    comparisons : Optional[List[str]]
        List of comparison pairs. If None, extract all automatically
    n_genes : Optional[int]
        Number of genes per comparison. If None, keep all genes
    groupby : Optional[str]
        Groupby parameter. If None, will try to auto-detect from adata.obs or set to None
        
    Returns
    -------
    ad.AnnData
        Modified AnnData object (modified in place)
    """
    # Auto-detect groupby if not specified
    if groupby is None and hasattr(adata, 'obs') and len(adata.obs.columns) > 0:
        # Try common column names
        for col in ['group', 'groups', 'condition', 'cluster']:
            if col in adata.obs.columns:
                groupby = col
                break
    
    rgg_dict = create_rank_genes_groups_structure(
        edger_df, 
        comparisons=comparisons,
        n_genes=n_genes,
        groupby=groupby
    )
    adata.uns[key] = rgg_dict
    return adata


def group_comparisons_by_reference(comparisons: List[str]) -> Dict[str, List[str]]:
    """Group comparison pairs by reference group.
    
    Parameters
    ----------
    comparisons : List[str]
        List of comparison pairs in format 'cond_i_vs_cond_j'
        
    Returns
    -------
    Dict[str, List[str]]
        Dictionary mapping reference group names to their comparison pairs
    """
    grouped = {}
    for comp in comparisons:
        _, ref_group = parse_comparison_pair(comp)
        ref_group = ref_group or 'unknown'
        if ref_group not in grouped:
            grouped[ref_group] = []
        grouped[ref_group].append(comp)
    return grouped


def create_multiple_rank_genes_groups(
    edger_df: pd.DataFrame,
    comparisons: Optional[List[str]] = None,
    n_genes: Optional[int] = None,
    group_by_reference: bool = True,
    groupby: Optional[str] = None
) -> Dict[str, Dict]:
    """Create separate rank_genes_groups results for multiple reference groups.
    
    Parameters
    ----------
    edger_df : pd.DataFrame
        edgeR output DataFrame
    comparisons : Optional[List[str]]
        List of comparison pairs. If None, extract all automatically
    n_genes : Optional[int]
        Number of genes per comparison. If None, keep all genes
    group_by_reference : bool
        If True, group by reference; if False, use single result for all
    groupby : Optional[str]
        Groupby parameter for params. If None, will be set to None
        
    Returns
    -------
    Dict[str, Dict]
        Dictionary mapping reference group names to rank_genes_groups dictionaries
    """
    comparison_dict = parse_comparison_columns(edger_df)
    
    if comparisons is None:
        comparisons = sorted(comparison_dict.keys())
    
    required_metrics = ['logCPM', 'PValue', 'FDR', 'logFC']
    valid_comparisons = []
    for comp in comparisons:
        if comp in comparison_dict:
            metrics = comparison_dict[comp]
            if all(metric in metrics for metric in required_metrics):
                valid_comparisons.append(comp)
            else:
                missing = [m for m in required_metrics if m not in metrics]
                print(f"Warning: Comparison '{comp}' missing metrics: {missing}")
        else:
            print(f"Warning: Comparison '{comp}' not found in DataFrame")
    
    if len(valid_comparisons) == 0:
        raise ValueError("No valid comparisons found")
    
    comparisons = valid_comparisons
    
    if not group_by_reference:
        result = create_rank_genes_groups_structure(
            edger_df,
            comparisons=comparisons,
            n_genes=n_genes,
            groupby=groupby
        )
        return {'all': result}
    
    grouped_comparisons = group_comparisons_by_reference(comparisons)
    results = {}
    for ref_group, comp_list in grouped_comparisons.items():
        print(f"Processing reference group '{ref_group}': {len(comp_list)} comparisons")
        result = create_rank_genes_groups_structure(
            edger_df,
            comparisons=comp_list,
            n_genes=n_genes,
            groupby=groupby
        )
        results[ref_group] = result
    
    return results


def add_multiple_rank_genes_groups_to_anndata(
    adata: ad.AnnData,
    edger_df: pd.DataFrame,
    key_prefix: str = 'rank_genes_groups',
    comparisons: Optional[List[str]] = None,
    n_genes: Optional[int] = None,
    group_by_reference: bool = True,
    groupby: Optional[str] = None
) -> ad.AnnData:
    """Add multiple rank_genes_groups results to AnnData with different keys per reference group.
    
    Parameters
    ----------
    adata : ad.AnnData
        AnnData object to add results to
    edger_df : pd.DataFrame
        edgeR output DataFrame with gene names as index
    key_prefix : str
        Prefix for keys in uns, default 'rank_genes_groups'
        Actual keys will be '{key_prefix}_{reference_group}'
    comparisons : Optional[List[str]]
        List of comparison pairs. If None, extract all automatically
    n_genes : Optional[int]
        Number of genes per comparison. If None, keep all genes
    group_by_reference : bool
        If True, group by reference; if False, use single result for all
    groupby : Optional[str]
        Groupby parameter. If None, will try to auto-detect from adata.obs or set to None
        
    Returns
    -------
    ad.AnnData
        Modified AnnData object (modified in place)
    """
    # Auto-detect groupby if not specified
    if groupby is None and hasattr(adata, 'obs') and len(adata.obs.columns) > 0:
        for col in ['group', 'groups', 'condition', 'cluster']:
            if col in adata.obs.columns:
                groupby = col
                break
    
    results = create_multiple_rank_genes_groups(
        edger_df,
        comparisons=comparisons,
        n_genes=n_genes,
        group_by_reference=group_by_reference,
        groupby=groupby
    )
    
    for ref_group, rgg_dict in results.items():
        if group_by_reference and ref_group != 'unknown':
            key = f"{key_prefix}_{ref_group}"
        else:
            key = key_prefix
        adata.uns[key] = rgg_dict
        print(f"Added results to adata.uns['{key}']")
    
    return adata


def create_separate_anndata_by_reference(
    base_adata: ad.AnnData,
    edger_df: pd.DataFrame,
    comparisons: Optional[List[str]] = None,
    n_genes: Optional[int] = None,
    groupby: Optional[str] = None
) -> Dict[str, ad.AnnData]:
    """Create separate AnnData objects for each reference group.
    
    Parameters
    ----------
    base_adata : ad.AnnData
        Base AnnData object to use as template
    edger_df : pd.DataFrame
        edgeR output DataFrame with gene names as index
    comparisons : Optional[List[str]]
        List of comparison pairs. If None, extract all automatically
    n_genes : Optional[int]
        Number of genes per comparison. If None, keep all genes
    groupby : Optional[str]
        Groupby parameter. If None, will try to auto-detect from base_adata.obs or set to None
        
    Returns
    -------
    Dict[str, ad.AnnData]
        Dictionary mapping reference group names to AnnData objects
    """
    # Auto-detect groupby if not specified
    if groupby is None and hasattr(base_adata, 'obs') and len(base_adata.obs.columns) > 0:
        for col in ['group', 'groups', 'condition', 'cluster']:
            if col in base_adata.obs.columns:
                groupby = col
                break
    
    results = create_multiple_rank_genes_groups(
        edger_df,
        comparisons=comparisons,
        n_genes=n_genes,
        group_by_reference=True,
        groupby=groupby
    )
    
    anndata_dict = {}
    for ref_group, rgg_dict in results.items():
        new_adata = base_adata.copy()
        new_adata.uns['rank_genes_groups'] = rgg_dict
        anndata_dict[ref_group] = new_adata
        print(f"Created AnnData object for reference group: '{ref_group}'")
    
    return anndata_dict
