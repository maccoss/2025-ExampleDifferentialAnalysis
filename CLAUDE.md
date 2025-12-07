# CLAUDE.md - AI Assistant Guide for proteomics_toolkit

This document helps AI assistants (Claude, Copilot, etc.) understand how to work with the `proteomics_toolkit` library for mass spectrometry-based proteomics data analysis.

## Project Overview

The `proteomics_toolkit` is a Python library for analyzing quantitative proteomics data, particularly from Skyline outputs. It provides a complete analysis pipeline from data import through statistical analysis and visualization.

## Key Modules and Their Purposes

### 1. `data_import` - Data Loading
```python
import proteomics_toolkit as ptk

# Load protein data and metadata
protein_data, metadata, peptide_data = ptk.load_skyline_data(
    protein_file="protein_quant.csv",
    metadata_file="metadata.csv",
    peptide_file=None  # Optional
)

# Clean sample names (removes common prefixes/suffixes)
cleaned_names = ptk.clean_sample_names(sample_columns, auto_detect=True)
```

### 2. `preprocessing` - Data Preparation
```python
# Parse protein identifiers (UniProt accessions, gene names)
processed_data = ptk.preprocessing.parse_protein_identifiers(protein_data)
processed_data = ptk.preprocessing.parse_gene_and_description(processed_data)

# Filter proteins by detection rate
filtered_data = ptk.preprocessing.filter_proteins_by_completeness(
    processed_data, sample_columns, min_detection_rate=0.5
)

# Classify samples into groups
group_distribution, control_samples, study_samples, sample_metadata, group_colors = ptk.classify_samples(
    sample_metadata=sample_metadata,
    group_column="Group",
    group_labels=["Control", "Treatment"],
    control_column="Sample_Type",
    control_labels=["Pool", "QC"]
)
```

### 3. `normalization` - Data Normalization
Available methods (choose ONE):
```python
# Median normalization - RECOMMENDED for most analyses
normalized = ptk.normalization.median_normalize(data, sample_columns=sample_columns)

# VSN - Variance Stabilizing Normalization (produces arcsinh-transformed data)
normalized = ptk.normalization.vsn_normalize(data, sample_columns=sample_columns)

# Other options: quantile_normalize, mad_normalize, z_score_normalize, rlr_normalize, loess_normalize
```

### 4. `statistical_analysis` - Statistical Testing

#### Configuration Pattern
```python
config = ptk.StatisticalConfig()

# Choose analysis type based on experimental design:
# - "paired": Before/after comparison in same subjects
# - "unpaired": Independent group comparison
# - "longitudinal": F-test for ANY change over time (categorical timepoints)
# - "linear_trend": Test for linear dose-response or time trend (continuous predictor)
# - "interaction": Group × Time interaction effects

config.analysis_type = "longitudinal"  # or "linear_trend", "paired", "unpaired"
config.statistical_test_method = "mixed_effects"  # RECOMMENDED for complex designs
config.subject_column = "Subject"  # For random effects
config.time_column = "Week"  # Time/dose variable

# For longitudinal: time treated as CATEGORICAL (F-test: do any timepoints differ?)
# For linear_trend: time treated as CONTINUOUS (tests if slope ≠ 0)

# Important: Empty covariates for longitudinal to avoid confounding with random intercept
config.covariates = []  # Between-subject variables absorbed by (1|Subject)
```

#### Running Analysis
```python
results = ptk.run_comprehensive_statistical_analysis(
    normalized_data=normalized_data,
    sample_metadata=sample_metadata,
    config=config,
    protein_annotations=filtered_data  # Optional: for Gene names
)

# Display summary
ptk.display_analysis_summary(results, config, label_top_n=25)
```

### 5. `enrichment` - Gene Set Enrichment Analysis

The enrichment module is a **general-purpose** module for gene set enrichment analysis via the Enrichr API. It can be used with any gene list, not just temporal data:

```python
import proteomics_toolkit as ptk
from proteomics_toolkit.enrichment import EnrichmentConfig

# Configure enrichment
enrich_config = EnrichmentConfig(
    enrichr_libraries=['GO_Biological_Process_2023', 'KEGG_2021_Human', 'Reactome_2022'],
    pvalue_cutoff=0.05,
    top_n=20,
    min_genes=5
)

# Option 1: Enrichment on any gene list
gene_list = ['TP53', 'BRCA1', 'BRCA2', 'ATM', 'CHEK2']
enrichment_df = ptk.run_enrichment_analysis(gene_list, config=enrich_config)

# Option 2: Enrichment on up/down regulated genes from differential analysis
enrichment_by_direction = ptk.run_differential_enrichment(
    stats_results,
    logfc_threshold=1.0,  # 2-fold change
    pvalue_threshold=0.05,
    config=enrich_config
)

# Option 3: Enrichment by any group (clusters, categories, treatments)
enrichment_by_cluster = ptk.run_enrichment_by_group(
    clustered_df,
    group_column='Cluster_Name',  # or 'Treatment', 'Response_Type', etc.
    gene_column='Gene',
    config=enrich_config
)

# Visualize
fig = ptk.plot_enrichment_barplot(enrichment_df, title='Enriched Pathways')
fig = ptk.plot_enrichment_comparison(enrichment_by_direction, title='Up vs Down')
```

### 6. `temporal_clustering` - Longitudinal Pattern Analysis

For time-series data with multiple timepoints:
```python
from proteomics_toolkit.temporal_clustering import (
    TemporalClusteringConfig,
    run_temporal_analysis
)

# Configure
temporal_config = TemporalClusteringConfig(
    auto_detect_clusters=True,  # Automatic optimal k via silhouette scores
    min_clusters=2,
    max_clusters=8,
    p_value_threshold=0.05,
    # Gene set enrichment libraries (uses enrichment module internally)
    enrichr_libraries=[
        'GO_Biological_Process_2023',
        'KEGG_2021_Human',
        'Reactome_2022'
    ]
)

# Run complete pipeline
temporal_results = run_temporal_analysis(
    data_df=normalized_data,
    metadata_dict=sample_metadata,
    stats_df=statistical_results,
    treatment_name='Treatment',
    config=temporal_config,
    run_enrichment=True,  # Query Enrichr API
    output_prefix='Temporal-Analysis'
)

# Access results
sig_proteins = temporal_results['significant_df']
heatmap_fig = temporal_results['fig_heatmap']
enrichment = temporal_results['enrichment_results']
```

### 7. `visualization` - Plots and Figures

**Standard QC and Results Plots:**
```python
# Volcano plot
ptk.plot_volcano(results, fc_threshold=0.585, p_threshold=0.05, label_top_n=25)

# Box plots
ptk.plot_box_plot(data, sample_columns, sample_metadata, group_colors)

# PCA comparison
ptk.plot_comparative_pca(original_data, median_normalized, vsn_normalized, ...)
```

**Grouped Data Visualizations (for any categorical grouping):**
```python
# Heatmap organized by any group (clusters, doses, treatments)
fig = ptk.plot_grouped_heatmap(
    data_df=merged_df,
    value_columns=['Week_0', 'Week_2', 'Week_4', 'Week_8'],  # or doses, conditions
    group_column='Cluster_Name',  # or 'Treatment', 'Response_Type'
    label_column='Gene',
    title='Expression Heatmap by Group'
)

# Trajectory plots for any grouped data
fig = ptk.plot_grouped_trajectories(
    data_df=merged_df,
    value_columns=['Dose_0', 'Dose_10', 'Dose_50', 'Dose_100'],  # or weeks
    group_column='Sensitivity',  # or 'Cluster', 'Treatment'
    x_values=[0, 10, 50, 100],  # Numeric x-axis values
    x_label='Dose (mg)',
    title='Dose Response by Sensitivity Group'
)
```

### 8. `export` - Save Results
```python
ptk.export_complete_analysis(
    normalized_data=normalized_data,
    sample_metadata=sample_metadata,
    config_dict=config_dict,
    differential_results=results,
    output_prefix="MyAnalysis"
)
```

## Common Patterns

### Module Reloading (for development)
When working in notebooks and modifying toolkit code:
```python
import sys
import importlib

if 'proteomics_toolkit' in sys.modules:
    importlib.reload(sys.modules['proteomics_toolkit'])
    submodules = ['statistical_analysis', 'visualization', 'normalization', 
                  'preprocessing', 'data_import', 'export', 'validation', 
                  'temporal_clustering']
    for submodule in submodules:
        module_name = f'proteomics_toolkit.{submodule}'
        if module_name in sys.modules:
            importlib.reload(sys.modules[module_name])
```

### Sample Metadata Structure
Sample metadata is typically a dictionary mapping sample names to their properties:
```python
sample_metadata = {
    'Sample_001': {'Subject': 'S1', 'Week': 0, 'Group': 'Treatment'},
    'Sample_002': {'Subject': 'S1', 'Week': 4, 'Group': 'Treatment'},
    # ...
}
```

### Statistical Analysis Types Reference

| Analysis Type | When to Use | Model | Test |
|--------------|-------------|-------|------|
| `paired` | Before/after in same subjects | Mixed effects with interaction | Group × Time interaction |
| `unpaired` | Independent groups | T-test or Welch | Group difference |
| `longitudinal` | Multiple timepoints, any change | `Protein ~ C(Week) + (1|Subject)` | F-test on time (categorical) |
| `linear_trend` | Dose-response or monotonic time | `Protein ~ Time + (1|Subject)` | Slope ≠ 0 (continuous) |
| `interaction` | Groups respond differently over time | `Protein ~ Group * Time + (1|Subject)` | Interaction term |

## Important Notes for AI Assistants

1. **Normalization Methods**: VSN produces arcsinh-transformed data (may include negatives); Median preserves original scale
2. **Log Transformation**: Set `config.log_transform_before_stats = "auto"` to automatically detect if needed
3. **Covariates Warning**: For longitudinal/repeated measures, between-subject covariates (Age, Sex) confound with the random intercept. Keep `config.covariates = []`
4. **P-value Selection**: Use `use_adjusted_pvalue = "adjusted"` for FDR-corrected values, or `"unadjusted"` for raw p-values
5. **Fold Change Threshold**: Typically expressed as log2. `np.log2(1.5) ≈ 0.585` for 1.5-fold change

## Visualization Guidelines

**IMPORTANT: Default font sizes should be LARGE for readability.** Scientific figures are often viewed in presentations, posters, and publications where small text is difficult to read.

### Recommended Default Font Sizes
- **Title**: 18-24 pt (bold)
- **Axis labels**: 14-18 pt (bold for important labels)
- **Tick labels**: 12-14 pt
- **Legend text**: 12-14 pt
- **Legend title**: 14-16 pt (bold)
- **Annotations**: 10-12 pt

### Key Principles
1. **Start with larger fonts** - it's easier to reduce than to increase later
2. **Use bold for titles and important labels** - improves visibility
3. **Position legends to avoid overlap** - use `bbox_to_anchor` to place outside plot area
4. **Use `tight_layout(rect=[...])` or `constrained_layout=True`** to accommodate legends outside the figure
5. **Test readability** at the expected viewing size (often smaller than on-screen)

### Example Pattern
```python
# Good defaults for most plots
label_fontsize = 14
title_fontsize = 18
tick_fontsize = 12
legend_fontsize = 12

ax.set_title("Title", fontsize=title_fontsize, fontweight='bold')
ax.set_xlabel("X Label", fontsize=label_fontsize, fontweight='bold')
ax.set_ylabel("Y Label", fontsize=label_fontsize, fontweight='bold')
ax.tick_params(labelsize=tick_fontsize)

# Legend outside plot to avoid overlap
fig.legend(loc='upper left', bbox_to_anchor=(1.02, 0.98), fontsize=legend_fontsize)
plt.tight_layout(rect=[0, 0, 0.85, 1])  # Leave space for legend
```

## Example Notebooks

- **`Minimal-Proteomics-Analysis.ipynb`**: Basic workflow for dose-response or group comparisons
- **`Verapamil-Liraglutide-Separate-Analysis.ipynb`**: Advanced longitudinal analysis with temporal clustering

## Running Tests

```bash
pytest tests/ -v
```

## Dependencies

- pandas, numpy, scipy, matplotlib, seaborn
- scikit-learn (for clustering)
- statsmodels (for mixed-effects models)
- requests (for Enrichr API)
