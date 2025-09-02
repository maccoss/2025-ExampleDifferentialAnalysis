"""
Data Import Module for Proteomics Analysis Toolkit

Functions for loading and parsing Skyline quantitation files and metadata.
"""

import pandas as pd
import re
import os
from typing import Tuple, Dict, Any, Optional, List


def load_skyline_data(protein_file: str, metadata_file: str, 
                     peptide_file: Optional[str] = None) -> Tuple[pd.DataFrame, pd.DataFrame, Optional[pd.DataFrame]]:
    """
    Load Skyline protein quantitation data and metadata.
    
    Parameters:
    -----------
    protein_file : str
        Path to protein quantitation CSV file
    metadata_file : str  
        Path to sample metadata CSV file
    peptide_file : str, optional
        Path to peptide quantitation CSV file
    
    Returns:
    --------
    protein_data : pd.DataFrame
        Protein quantitation data
    metadata : pd.DataFrame
        Sample metadata
    peptide_data : pd.DataFrame or None
        Peptide quantitation data if provided
    """
    
    print("=== LOADING SKYLINE DATA ===\n")
    
    # Check if files exist
    for file_path, file_type in [(protein_file, "protein"), (metadata_file, "metadata")]:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"{file_type.title()} file not found: {file_path}")
    
    # Load protein data
    try:
        protein_data = pd.read_csv(protein_file)
        print(f"✓ Loaded protein data: {protein_data.shape}")
    except Exception as e:
        raise ValueError(f"Error loading protein file: {e}")
    
    # Load metadata
    try:
        metadata = pd.read_csv(metadata_file)
        print(f"✓ Loaded metadata: {metadata.shape}")
    except Exception as e:
        raise ValueError(f"Error loading metadata file: {e}")
    
    # Load peptide data if provided
    peptide_data = None
    if peptide_file:
        if os.path.exists(peptide_file):
            try:
                peptide_data = pd.read_csv(peptide_file)
                print(f"✓ Loaded peptide data: {peptide_data.shape}")
            except Exception as e:
                print(f"Warning: Could not load peptide file: {e}")
        else:
            print(f"Warning: Peptide file not found: {peptide_file}")
    
    print("\nData loading completed successfully!")
    return protein_data, metadata, peptide_data


def parse_uniprot_identifier(protein_id: str) -> Dict[str, str]:
    """
    Parse UniProt identifier from protein column.
    
    Handles formats like: sp|P12345|PROT_HUMAN -> P12345
    
    Parameters:
    -----------
    protein_id : str
        Protein identifier string
        
    Returns:
    --------
    dict with keys: accession, database, entry_name
    """
    if pd.isna(protein_id):
        return {'accession': '', 'database': '', 'entry_name': ''}
    
    protein_id = str(protein_id).strip()
    
    # Pattern: sp|P12345|PROT_HUMAN or tr|Q9ABC1|Q9ABC1_MOUSE
    match = re.match(r'^(sp|tr)\|([A-Z0-9]+)\|([A-Z0-9_]+)$', protein_id)
    if match:
        db = 'SwissProt' if match.group(1) == 'sp' else 'TrEMBL'
        return {
            'accession': match.group(2),
            'database': db,
            'entry_name': match.group(3)
        }
    
    # If no match, try to find accession pattern
    acc_match = re.search(r'([A-Z][A-Z0-9]{5,9})', protein_id)
    if acc_match:
        return {
            'accession': acc_match.group(1),
            'database': '',
            'entry_name': ''
        }
    
    return {'accession': '', 'database': '', 'entry_name': ''}


def parse_gene_from_description(protein_description: str) -> str:
    """
    Extract gene name from protein description using GN= pattern.
    
    Parameters:
    -----------
    protein_description : str
        Protein description text
        
    Returns:
    --------
    str : Gene name or empty string
    """
    if pd.isna(protein_description):
        return ''
    
    desc_str = str(protein_description).strip()
    
    # Extract gene name from GN= pattern
    gene_match = re.search(r'GN=([^=\s]+)', desc_str)
    return gene_match.group(1).strip() if gene_match else ''


def clean_description(protein_description: str) -> str:
    """
    Clean protein description by removing UniProt annotation fields.
    
    Parameters:
    -----------
    protein_description : str
        Raw protein description
        
    Returns:
    --------
    str : Cleaned description
    """
    if pd.isna(protein_description):
        return ''
    
    desc_str = str(protein_description).strip()
    
    # Remove patterns like OS=..., OX=..., GN=..., PE=..., SV=...
    clean_desc = re.sub(r'\s+(?:OS|OX|GN|PE|SV)=[^=]*(?=\s+(?:OS|OX|GN|PE|SV)=|$)', '', desc_str)
    return clean_desc.strip()


def identify_sample_columns(data: pd.DataFrame, metadata: pd.DataFrame) -> List[str]:
    """
    Identify sample columns in quantitation data based on metadata.
    
    Parameters:
    -----------
    data : pd.DataFrame
        Quantitation data
    metadata : pd.DataFrame
        Sample metadata
        
    Returns:
    --------
    List[str] : Sample column names
    """
    # Look for columns that match sample names in metadata
    replicate_col = None
    for col in ['Replicate', 'Sample', 'Sample_Name', 'SampleName']:
        if col in metadata.columns:
            replicate_col = col
            break
    
    if replicate_col is None:
        # Try first column as sample names
        replicate_col = metadata.columns[0]
        print(f"Warning: Using first metadata column '{replicate_col}' as sample names")
    
    metadata_samples = set(metadata[replicate_col].astype(str))
    
    # Find sample columns - look for columns ending with common intensity suffixes
    sample_columns = []
    intensity_suffixes = ['Normalized Area', 'Total Area Ratio', 'Area', 'Intensity', 'Abundance']
    
    for col in data.columns:
        # Check if any sample name appears in this column name
        for sample in metadata_samples:
            if str(sample) in str(col):
                sample_columns.append(col)
                break
        # Also check for columns ending with intensity suffixes
        else:
            for suffix in intensity_suffixes:
                if col.endswith(suffix):
                    sample_columns.append(col)
                    break
    
    print(f"Identified {len(sample_columns)} sample columns")
    return sample_columns


def clean_sample_names(sample_columns: List[str], common_prefix: Optional[str] = None, 
                      common_suffix: Optional[str] = None) -> Dict[str, str]:
    """
    Clean sample names by removing common prefixes/suffixes.
    
    Parameters:
    -----------
    sample_columns : List[str]
        List of sample column names
    common_prefix : str, optional
        Common prefix to remove
    common_suffix : str, optional  
        Common suffix to remove
        
    Returns:
    --------
    Dict[str, str] : Mapping from original to cleaned names
    """
    
    cleaned_names = {}
    
    # Auto-detect common prefix if not provided
    if common_prefix is None:
        if len(sample_columns) > 1:
            # Find common prefix
            prefix = os.path.commonprefix(sample_columns)
            # Remove trailing non-alphanumeric characters
            prefix = re.sub(r'[^a-zA-Z0-9]+$', '', prefix)
            common_prefix = prefix if len(prefix) > 0 else ""
        else:
            common_prefix = ""
    
    # Auto-detect common suffix if not provided  
    if common_suffix is None:
        if len(sample_columns) > 1:
            # Find common suffix by reversing strings
            reversed_names = [name[::-1] for name in sample_columns]
            suffix = os.path.commonprefix(reversed_names)[::-1]
            # Remove leading non-alphanumeric characters
            suffix = re.sub(r'^[^a-zA-Z0-9]+', '', suffix)
            common_suffix = suffix if len(suffix) > 0 else ""
        else:
            common_suffix = ""
    
    print(f"Removing common prefix: '{common_prefix}'")
    print(f"Removing common suffix: '{common_suffix}'")
    
    for original_name in sample_columns:
        cleaned_name = original_name
        
        # Remove prefix
        if common_prefix and cleaned_name.startswith(common_prefix):
            cleaned_name = cleaned_name[len(common_prefix):]
        
        # Remove suffix
        if common_suffix and cleaned_name.endswith(common_suffix):
            cleaned_name = cleaned_name[:-len(common_suffix)]
        
        # Clean up any remaining non-alphanumeric characters at start/end
        cleaned_name = re.sub(r'^[^a-zA-Z0-9]+', '', cleaned_name)
        cleaned_name = re.sub(r'[^a-zA-Z0-9]+$', '', cleaned_name)
        
        cleaned_names[original_name] = cleaned_name
    
    return cleaned_names


def match_samples_to_metadata(cleaned_sample_names: Dict[str, str], 
                             metadata: pd.DataFrame) -> Dict[str, Dict[str, Any]]:
    """
    Match cleaned sample names to metadata entries.
    
    Parameters:
    -----------
    cleaned_sample_names : Dict[str, str]
        Mapping from original to cleaned sample names
    metadata : pd.DataFrame
        Sample metadata
        
    Returns:
    --------
    Dict[str, Dict] : Sample metadata mapping
    """
    
    # Find replicate column
    replicate_col = None
    for col in ['Replicate', 'Sample', 'Sample_Name', 'SampleName']:
        if col in metadata.columns:
            replicate_col = col
            break
    
    if replicate_col is None:
        replicate_col = metadata.columns[0]
    
    sample_metadata = {}
    matched_count = 0
    
    for original_name, cleaned_name in cleaned_sample_names.items():
        # Try to find exact match first
        metadata_match = metadata[metadata[replicate_col] == cleaned_name]
        
        if metadata_match.empty:
            # Try partial matching
            for _, row in metadata.iterrows():
                if cleaned_name in str(row[replicate_col]) or str(row[replicate_col]) in cleaned_name:
                    metadata_match = pd.DataFrame([row])
                    break
        
        if not metadata_match.empty:
            matched_count += 1
            sample_metadata[original_name] = metadata_match.iloc[0].to_dict()
        else:
            # Create placeholder metadata
            sample_metadata[original_name] = {
                str(replicate_col): cleaned_name,
                'Group': 'Unknown',
                'matched': False
            }
    
    print(f"Matched {matched_count}/{len(cleaned_sample_names)} samples to metadata")
    
    return sample_metadata


def identify_and_classify_controls(sample_metadata: dict, metadata: Optional[pd.DataFrame] = None,
                                 control_keywords: Optional[Dict[str, List[str]]] = None,
                                 update_nan_only: bool = True) -> Tuple[dict, Dict[str, Any]]:
    """
    Automatically identify and reclassify control samples in sample metadata.
    
    This function searches through sample names and metadata columns for control-related 
    keywords and updates sample group assignments. It's designed to be flexible and work
    with various control sample naming conventions across different experiments.
    
    Parameters:
    -----------
    sample_metadata : dict
        Sample metadata dictionary to update
    metadata : pd.DataFrame, optional
        Original metadata DataFrame for additional context
    control_keywords : Dict[str, List[str]], optional
        Dictionary mapping control types to their keywords. Default includes:
        - pool: ['pool', 'pooled', 'poolsample']
        - qc: ['qc', 'quality', 'qualitycontrol'] 
        - reference: ['reference', 'ref', 'standard', 'std']
        - control: ['control', 'ctrl', 'cntrl']
        - blank: ['blank', 'buffer', 'negative']
        - spike: ['spike', 'spikein', 'positive']
    update_nan_only : bool, default True
        If True, only update samples with NaN/Unknown groups.
        If False, update all samples that match control patterns.
    
    Returns:
    --------
    updated_sample_metadata : dict
        Updated sample metadata with proper control classifications
    summary : Dict[str, Any]
        Summary including changes made, group counts, and identified patterns
    
    Example:
    --------
    >>> updated_metadata, summary = identify_and_classify_controls(sample_metadata)
    >>> print(f"Updated {summary['total_updated']} samples")
    >>> print(f"Found {len(summary['control_types'])} different control types")
    """
    
    if control_keywords is None:
        control_keywords = {
            'pool': ['pool', 'pooled', 'poolsample'],
            'qc': ['qc', 'quality', 'qualitycontrol'], 
            'reference': ['reference', 'ref', 'standard', 'std'],
            'control': ['control', 'ctrl', 'cntrl'],
            'blank': ['blank', 'buffer', 'negative'],
            'spike': ['spike', 'spikein', 'positive']
        }
    
    # Create a copy to avoid modifying the original
    updated_metadata = {k: v.copy() for k, v in sample_metadata.items()}
    changes = []
    identified_patterns = {}
    
    def classify_sample_name(sample_name: str) -> Optional[str]:
        """
        Classify a sample based on its name using hierarchical pattern matching.
        Returns specific pool names when possible, otherwise general categories.
        """
        sample_lower = sample_name.lower()
        
        # First, look for specific pool patterns (most specific)
        pool_patterns = {
            'eisaipool': 'EISAIPool',
            'gwpool': 'GWPool', 
            'hoof': 'HoofPool',  # Can be HoofPool, Hoof_Pool, etc.
        }
        
        for pattern, group_name in pool_patterns.items():
            if pattern in sample_lower:
                identified_patterns[pattern] = identified_patterns.get(pattern, 0) + 1
                return group_name
        
        # Then look for general control categories
        for control_type, keywords in control_keywords.items():
            for keyword in keywords:
                if keyword in sample_lower:
                    identified_patterns[keyword] = identified_patterns.get(keyword, 0) + 1
                    # Return capitalized version of control type
                    if control_type == 'pool':
                        return 'Pool'
                    elif control_type == 'qc':
                        return 'QC'
                    elif control_type == 'reference':
                        return 'Reference'
                    elif control_type == 'control':
                        return 'Control'
                    elif control_type == 'blank':
                        return 'Blank'
                    elif control_type == 'spike':
                        return 'Spike'
        
        return None
    
    def classify_from_metadata(sample_name: str, metadata_row: pd.Series) -> Optional[str]:
        """Check metadata columns for control indicators"""
        if metadata is None:
            return None
            
        for col_name, value in metadata_row.items():
            if pd.isna(value) or value == '':
                continue
                
            value_str = str(value).lower()
            
            # Check for specific patterns first
            if 'eisaipool' in value_str:
                return 'EISAIPool'
            elif 'gwpool' in value_str:
                return 'GWPool'
            elif 'hoof' in value_str and 'pool' in value_str:
                return 'HoofPool'
            
            # Check general patterns
            for control_type, keywords in control_keywords.items():
                for keyword in keywords:
                    if keyword in value_str:
                        if control_type == 'pool':
                            return 'Pool'
                        elif control_type == 'qc':
                            return 'QC'
                        elif control_type == 'reference':
                            return 'Reference'
                        elif control_type == 'control':
                            return 'Control'
                        elif control_type == 'blank':
                            return 'Blank'
                        elif control_type == 'spike':
                            return 'Spike'
        
        return None
    
    # Process each sample
    for sample_name, sample_info in updated_metadata.items():
        original_group = sample_info.get('Group', 'Unknown')
        
        # Skip if update_nan_only is True and sample has a valid group
        if update_nan_only and pd.notna(original_group) and str(original_group) != 'Unknown':
            continue
        
        # Try classification from sample name first
        new_group = classify_sample_name(sample_name)
        
        # If not found and metadata is available, check metadata
        if new_group is None and metadata is not None:
            # Find matching row in metadata (try different approaches)
            metadata_row = None
            for col in metadata.columns:
                matching_rows = metadata[metadata[col].astype(str) == sample_name]
                if len(matching_rows) > 0:
                    metadata_row = matching_rows.iloc[0]
                    break
            
            if metadata_row is not None:
                new_group = classify_from_metadata(sample_name, metadata_row)
        
        # Update if we found a classification
        if new_group and new_group != original_group:
            updated_metadata[sample_name]['Group'] = new_group
            changes.append({
                'sample': sample_name,
                'original': str(original_group) if pd.notna(original_group) else 'NaN',
                'new': new_group
            })
    
    # Calculate summary statistics
    group_counts = {}
    control_types = set()
    
    for sample_info in updated_metadata.values():
        group = sample_info.get('Group', 'Unknown')
        group_counts[group] = group_counts.get(group, 0) + 1
        
        # Track control types
        group_lower = str(group).lower()
        for control_type in control_keywords.keys():
            if any(keyword in group_lower for keyword in control_keywords[control_type]):
                control_types.add(control_type)
                break
        
        # Also check for specific pools
        if group in ['EISAIPool', 'GWPool', 'HoofPool', 'Pool']:
            control_types.add('pool')
    
    summary = {
        'total_updated': len(changes),
        'changes': changes,
        'group_counts': group_counts,
        'control_types': sorted(list(control_types)),
        'identified_patterns': identified_patterns,
        'total_samples': len(updated_metadata),
        'control_samples': sum(count for group, count in group_counts.items() 
                             if group not in ['Unknown'] and 
                             any(keyword in str(group).lower() 
                                 for keywords in control_keywords.values() 
                                 for keyword in keywords)),
        'study_samples': sum(count for group, count in group_counts.items() 
                           if group not in ['Unknown', 'QC', 'Pool', 'EISAIPool', 'GWPool', 
                                          'HoofPool', 'Reference', 'Control', 'Blank', 'Spike'])
    }
    
    return updated_metadata, summary
