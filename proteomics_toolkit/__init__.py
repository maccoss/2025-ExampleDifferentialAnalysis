"""
Proteomics Analysis Toolkit
===========================

A comprehensive Python library for analyzing mass spectrometry-based proteomics data,
particularly designed for Skyline quantitation outputs. This toolkit provides a complete
workflow from data import through statistical analysis and visualization.

QUICK START EXAMPLE:
-------------------
    import proteomics_toolkit as ptk
    
    # 1. Load data
    protein_data, metadata, _ = ptk.load_skyline_data('proteins.csv', 'metadata.csv')
    
    # 2. Process and normalize
    processed = ptk.parse_protein_identifiers(protein_data)
    normalized = ptk.median_normalize(processed, sample_columns=['Sample1', 'Sample2'])
    
    # 3. Statistical analysis
    config = ptk.StatisticalConfig()
    results = ptk.run_comprehensive_statistical_analysis(normalized, sample_metadata, config)
    
    # 4. Visualization and export
    ptk.plot_volcano(results)
    ptk.export_complete_analysis(normalized, sample_metadata, config, results)

MODULE OVERVIEW:
===============
    
data_import
    Purpose: Load and parse Skyline quantitation files and metadata
    Key functions: load_skyline_data(), clean_sample_names()
    Use when: Starting analysis, need to load protein/peptide data
    
preprocessing  
    Purpose: Data cleaning, quality assessment, and sample classification
    Key functions: parse_protein_identifiers(), classify_samples(), filter_proteins_by_completeness()
    Use when: Need to clean data, assess completeness, classify control vs study samples
    
normalization
    Purpose: Various normalization methods to reduce technical variation
    Key functions: median_normalize(), vsn_normalize(), quantile_normalize()
    Use when: Raw data needs normalization before statistical analysis
    
statistical_analysis
    Purpose: Differential analysis, mixed-effects models, and statistical testing
    Key functions: run_comprehensive_statistical_analysis(), StatisticalConfig()
    Use when: Performing statistical comparisons between groups
    
visualization
    Purpose: Quality control plots, results visualization, and publication-ready figures  
    Key functions: plot_volcano(), plot_box_plot(), plot_pca()
    Use when: Need QC plots, want to visualize results, create figures
    
validation
    Purpose: Data validation, error checking, and diagnostic reporting
    Key functions: validate_metadata_data_consistency(), enhanced_sample_processing()
    Use when: Need to validate metadata matches data, troubleshoot sample mismatches
    
export
    Purpose: Export results, configurations, and create reproducible analysis records
    Key functions: export_complete_analysis(), export_timestamped_config()
    Use when: Saving results, creating configuration backups, sharing analysis

TYPICAL WORKFLOW:
================
1. ptk.load_skyline_data() → Load protein data and metadata
2. ptk.validate_metadata_data_consistency() → Validate data consistency  
3. ptk.parse_protein_identifiers() → Clean and annotate protein data
4. ptk.classify_samples() → Identify control vs study samples
5. ptk.median_normalize() → Normalize data (or other normalization methods)
6. ptk.run_comprehensive_statistical_analysis() → Perform statistical analysis
7. ptk.plot_volcano() → Visualize results
8. ptk.export_complete_analysis() → Export everything for reproducibility

ERROR HANDLING:
==============
The toolkit includes comprehensive validation with clear error messages:
- SampleMatchingError: When samples in metadata don't match protein data
- ControlSampleError: When control samples are missing or misspecified
- Use ptk.validate_metadata_data_consistency() to diagnose issues early
"""

# =============================================================================
# MODULE IMPORTS - Core functionality organized by analysis stage
# =============================================================================

from . import data_import         # Data loading and parsing
from . import preprocessing       # Data cleaning and quality assessment  
from . import normalization       # Normalization methods
from . import statistical_analysis # Statistical testing and modeling
from . import visualization       # Plotting and visualization
from . import validation         # Data validation and error checking
from . import export             # Results export and configuration management

__version__ = "1.0.0"
__author__ = "Michael MacCoss Lab, University of Washington"

# =============================================================================
# CONVENIENCE IMPORTS - Most commonly used functions available at top level
# =============================================================================

# DATA LOADING - Essential functions for starting any analysis
from .data_import import (
    load_skyline_data,        # Main function: Load protein/peptide data + metadata
    clean_sample_names        # Clean up sample column names automatically
)

# DATA PREPROCESSING - Core data preparation functions  
from .preprocessing import (
    parse_protein_identifiers,    # Extract UniProt IDs, gene names, descriptions
    classify_samples,            # Classify samples into study vs control groups
    apply_systematic_color_scheme # Apply consistent colors for visualization
)

# NORMALIZATION - All normalization methods for reducing technical variation
from .normalization import (
    median_normalize,      # Most common: Simple, robust median normalization
    vsn_normalize,        # Advanced: Variance stabilizing normalization
    quantile_normalize,   # Strong: Force identical sample distributions  
    mad_normalize,        # Robust: Median absolute deviation normalization
    z_score_normalize,    # Standard: Z-score standardization
    rlr_normalize,        # Advanced: Robust linear regression
    loess_normalize,      # Advanced: LOESS intensity-dependent correction
    handle_negative_values, # Handle negative values from VSN normalization
    analyze_negative_values # Analyze negative value patterns
)

# STATISTICAL ANALYSIS - Core statistical functions and configuration
from .statistical_analysis import (
    run_comprehensive_statistical_analysis, # Main function: Complete statistical analysis
    display_analysis_summary,              # Display analysis results summary
    StatisticalConfig,                     # Configuration class for analysis parameters
    run_statistical_analysis              # Lower-level statistical analysis function
)

# DATA VALIDATION - Comprehensive data validation and error handling
from .validation import (
    validate_metadata_data_consistency,      # Main validation: Check metadata vs data consistency
    enhanced_sample_processing,             # Sample processing with built-in validation
    generate_sample_matching_diagnostic_report, # Detailed diagnostic reports
    SampleMatchingError,                    # Exception: Samples don't match between files
    ControlSampleError                      # Exception: Control sample configuration issues
)

# DATA EXPORT - Save results and create reproducible analysis records
from .export import (
    export_complete_analysis,              # Main function: Export everything (data + config)
    export_analysis_results,              # Export data files only
    export_timestamped_config,            # Export configuration with timestamp
    create_config_dict_from_notebook_vars, # Create config from notebook variables
    export_significant_proteins_summary,   # Export summary of significant results
    export_results                        # General results export function
)

# VISUALIZATION - Publication-ready plots and quality control visualizations
from .visualization import (
    plot_volcano,                        # Main results plot: Volcano plot for differential analysis
    plot_box_plot,                      # QC plot: Sample intensity distributions
    plot_normalization_comparison,       # QC plot: Before/after normalization comparison
    plot_pca,                           # QC plot: Principal component analysis
    plot_comparative_pca,               # Advanced: Compare multiple normalization methods
    plot_control_correlation_analysis,   # QC plot: Control sample correlation analysis
    plot_control_group_correlation_analysis, # QC plot: Group-wise control correlations
    plot_individual_control_pool_analysis    # Detailed: Individual control pool analysis
)

# =============================================================================
# PUBLIC API - All functions available for import
# =============================================================================

__all__ = [
    # MODULES - Full module access for advanced workflows
    "data_import",              # Complete data import functionality
    "preprocessing",            # Complete preprocessing functionality  
    "normalization",            # All normalization methods
    "statistical_analysis",     # Complete statistical analysis suite
    "visualization",            # All plotting and visualization functions
    "validation",              # Data validation and error handling
    "export",                  # Results export and configuration management
    
    # DATA LOADING - Start here for any analysis
    "load_skyline_data",        # ESSENTIAL: Load protein data + metadata
    "clean_sample_names",       # Clean and standardize sample column names
    
    # PREPROCESSING - Data preparation and quality control
    "parse_protein_identifiers", # Extract protein annotations (UniProt, Gene, etc.)
    "classify_samples",         # Separate control vs study samples
    "apply_systematic_color_scheme", # Consistent colors for groups
    
    # NORMALIZATION - Choose the best method for your data
    "median_normalize",         # Good Starting Point: Robust, preserves scale
    "vsn_normalize",           # Use for highly variable samples: Handles heteroscedasticity
    "quantile_normalize",      # Strong normalization (identical distributions)
    "mad_normalize",           # Robust to outliers
    "z_score_normalize",       # Standardize to mean=0, std=1
    "rlr_normalize",           # Robust linear regression
    "loess_normalize",         # LOESS intensity-dependent correction
    "handle_negative_values",   # Handle VSN negative values
    "analyze_negative_values",  # Analyze negative value patterns
    
    # STATISTICAL ANALYSIS - The core of differential analysis
    "run_comprehensive_statistical_analysis", # MAIN FUNCTION: Complete analysis
    "display_analysis_summary",              # Show results summary
    "StatisticalConfig",                     # Configuration for statistical analysis
    "run_statistical_analysis",             # Lower-level analysis function
    
    # VALIDATION - Catch errors early, get helpful diagnostics
    "validate_metadata_data_consistency",      # RECOMMENDED: Validate before analysis
    "enhanced_sample_processing",             # Sample processing with validation
    "generate_sample_matching_diagnostic_report", # Detailed error diagnostics
    "SampleMatchingError",                    # Custom exception for sample mismatches
    "ControlSampleError",                     # Custom exception for control sample issues
    
    # EXPORT - Save your work and create reproducible records
    "export_complete_analysis",              # MAIN FUNCTION: Export everything
    "export_analysis_results",              # Export data files only
    "export_timestamped_config",            # Export configuration backup
    "create_config_dict_from_notebook_vars", # Create config from notebook variables
    "export_significant_proteins_summary",   # Export significant results summary
    "export_results",                       # General export function
    
    # VISUALIZATION - Publication-ready plots and QC
    "plot_volcano",                        # ESSENTIAL: Volcano plot (main results plot)
    "plot_box_plot",                      # ESSENTIAL: Sample distributions (QC)
    "plot_normalization_comparison",       # Before/after normalization (QC)
    "plot_pca",                           # Principal component analysis (QC)
    "plot_comparative_pca",               # Compare normalization methods
    "plot_control_correlation_analysis",   # Control sample QC
    "plot_control_group_correlation_analysis", # Group-wise control QC
    "plot_individual_control_pool_analysis",   # Detailed control analysis
]
