"""
Proteomics Analysis Toolkit

A comprehensive library for analyzing mass spectrometry-based proteomics data,
particularly from Skyline outputs.

Main modules:
- data_import: Loading and parsing Skyline quantitation files
- preprocessing: Data cleaning and quality assessment
- normalization: Various normalization methods (median, VSN)
- statistical_analysis: All statistical testing and differential analysis
- visualization: Plotting functions for QC and results
- export: Data export and configuration management
"""

from . import data_import
from . import preprocessing
from . import normalization
from . import statistical_analysis
from . import visualization
from . import export

__version__ = "1.0.0"
__author__ = "Michael MacCoss, University of Washington"

# Main workflow functions for easy access
from .data_import import load_skyline_data, clean_sample_names
from .preprocessing import (
    parse_protein_identifiers,
    classify_samples,
    apply_systematic_color_scheme,
)
from .normalization import (
    median_normalize,
    vsn_normalize,
    quantile_normalize,
    mad_normalize,
    z_score_normalize,
    rlr_normalize,
    loess_normalize,
    handle_negative_values,
    analyze_negative_values,
)
from .statistical_analysis import (
    run_comprehensive_statistical_analysis,
    display_analysis_summary,
    StatisticalConfig,
    run_statistical_analysis,
)
from .export import (
    export_analysis_results,
    export_timestamped_config,
    export_complete_analysis,
    create_config_dict_from_notebook_vars,
    export_significant_proteins_summary,
    export_results,
)
from .visualization import (
    plot_box_plot,
    plot_volcano,
    plot_normalization_comparison,
    plot_pca,
    plot_comparative_pca,
    plot_control_correlation_analysis,
    plot_control_group_correlation_analysis,
    plot_individual_control_pool_analysis,
)

__all__ = [
    "data_import",
    "preprocessing",
    "normalization",
    "statistical_analysis",
    "visualization",
    "export",
    "load_skyline_data",
    "parse_protein_identifiers",
    "classify_samples",
    "apply_systematic_color_scheme",
    "clean_sample_names",
    "median_normalize",
    "vsn_normalize",
    "quantile_normalize",
    "mad_normalize",
    "z_score_normalize",
    "rlr_normalize",
    "loess_normalize",
    "handle_negative_values",
    "analyze_negative_values",
    "run_comprehensive_statistical_analysis",
    "display_analysis_summary",
    "export_results",
    "run_statistical_analysis",
    "StatisticalConfig",
    "export_analysis_results",
    "export_timestamped_config",
    "export_complete_analysis",
    "create_config_dict_from_notebook_vars",
    "export_significant_proteins_summary",
    "plot_box_plot",
    "plot_volcano",
    "plot_normalization_comparison",
    "plot_pca",
    "plot_comparative_pca",
    "plot_control_correlation_analysis",
    "plot_control_group_correlation_analysis",
    "plot_individual_control_pool_analysis",
]
