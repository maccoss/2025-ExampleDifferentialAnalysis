# Proteomics Analysis Toolkit

A comprehensive Python toolkit for analyzing mass spectrometry-based proteomics data, with a focus on Skyline quantitation workflows.

## Features

### 🔬 Core Analysis Modules
- **data_import**: Load Skyline data, parse protein identifiers, and manage sample metadata
- **preprocessing**: Data quality assessment, filtering, and annotation parsing  
- **normalization**: Multiple normalization methods (median, VSN, quantile)
- **differential**: Statistical analysis for differential protein expression
- **visualization**: Professional plots for data interpretation, including grouped heatmaps and trajectory plots
- **temporal_clustering**: Temporal trend analysis with K-means clustering
- **enrichment**: Gene set enrichment analysis via Enrichr API (works with any gene list)

### 🎯 Key Capabilities
- **Automated Control Detection**: Intelligent identification of QC, pool, and reference samples ([guide](README_control_identification.md))
- **Flexible Data Import**: Support for Skyline protein/peptide quantitation files
- **Multi-method Normalization**: Choose the best normalization for your data
- **Statistical Analysis**: Mixed-effects models with longitudinal analysis support
- **Temporal Clustering**: K-means clustering of temporal protein trends
- **Gene Set Enrichment**: Enrichr API integration for pathway analysis on any gene list
- **Grouped Visualizations**: Generalized heatmaps and trajectory plots for any grouped data
- **Professional Visualization**: Consistent color schemes and publication-ready plots

## Quick Start

### Basic Workflow
```python
import proteomics_toolkit as ptk

# 1. Load data
protein_data, metadata, peptide_data = ptk.load_skyline_data(
    protein_file="protein_quant.csv",
    metadata_file="metadata.csv"
)

# 2. Process and clean sample names
sample_columns = ptk.data_import.identify_sample_columns(protein_data, metadata)
cleaned_names = ptk.data_import.clean_sample_names(sample_columns)
sample_metadata = ptk.data_import.match_samples_to_metadata(cleaned_names, metadata)

# 3. Identify control samples automatically
sample_metadata, summary = ptk.data_import.identify_and_classify_controls(
    sample_metadata=sample_metadata,
    metadata=metadata
)

# 4. Process and filter data
processed_data = ptk.preprocessing.parse_protein_identifiers(protein_data)
filtered_data = ptk.preprocessing.filter_proteins_by_completeness(
    processed_data, list(cleaned_names.values()), min_detection_rate=0.5
)

# 5. Normalize data
normalized_data = ptk.normalization.median_normalize(
    filtered_data[list(cleaned_names.values())]
)

# 6. Statistical analysis
diff_results = ptk.differential.run_differential_analysis(
    data=pd.concat([filtered_data[['Protein', 'Gene']], normalized_data], axis=1),
    sample_metadata=sample_metadata,
    group1='Control', group2='Treatment'
)

# 7. Visualization
group_colors, group_counts = ptk.preprocessing.calculate_group_colors(sample_metadata)
ptk.visualization.plot_box_plot(
    data=pd.concat([filtered_data[['Protein']], normalized_data], axis=1),
    sample_columns=list(cleaned_names.values()),
    sample_metadata=sample_metadata,
    group_colors=group_colors
)

# 8. Quality Control - CV Distribution for Control Samples
cv_data = ptk.visualization.plot_control_cv_distribution(
    data=normalized_data,  # Use median normalized data for CV calculation
    sample_columns=list(cleaned_names.values()),
    sample_metadata=sample_metadata,
    control_column="Sample_Type",  # Or your control column name
    control_labels=["Pool", "QC", "Reference"],  # Your control labels
    normalization_method="Median",
    cv_threshold=20.0
)
```

## Advanced Features

### Control Sample Identification
The toolkit includes intelligent control sample identification that works across different experimental designs:

```python
# Automatic identification with defaults
sample_metadata, summary = ptk.data_import.identify_and_classify_controls(
    sample_metadata=sample_metadata,
    metadata=metadata
)

# Custom control patterns for your experiment
custom_keywords = {
    'pool': ['pool', 'pooled', 'reference_pool'],
    'qc': ['qc', 'quality', 'technical_replicate'],
    'control': ['control', 'vehicle', 'mock']
}

sample_metadata, summary = ptk.data_import.identify_and_classify_controls(
    sample_metadata=sample_metadata,
    metadata=metadata,
    control_keywords=custom_keywords
)
```

See [Control Identification Guide](README_control_identification.md) for detailed usage.

### Normalization Options
```python
# Median normalization (robust, interpretable)
normalized_data = ptk.normalization.median_normalize(raw_data)

# VSN normalization (variance stabilization)
normalized_data = ptk.normalization.vsn_normalize(raw_data, optimize_params=True)

# Quantile normalization (assumes similar distributions)  
normalized_data = ptk.normalization.quantile_normalize(raw_data)
```

### Statistical Analysis
```python
# Paired analysis (e.g., before/after treatment)
diff_results = ptk.differential.run_differential_analysis(
    data=analysis_data,
    sample_metadata=sample_metadata,
    group1='Baseline', group2='Treatment',
    analysis_type='paired',
    subject_column='Subject'
)

# Unpaired analysis with effect sizes
diff_results = ptk.differential.run_differential_analysis(
    data=analysis_data,
    sample_metadata=sample_metadata, 
    group1='Control', group2='Treatment',
    analysis_type='unpaired'
)
diff_results = ptk.differential.calculate_effect_sizes(diff_results)

# Longitudinal analysis with mixed-effects models
config = ptk.StatisticalConfig()
config.analysis_type = "longitudinal"  # F-test for any change over time
config.statistical_test_method = "mixed_effects"
config.subject_column = "Subject"
config.time_column = "Week"  # Categorical time variable
config.covariates = []  # Empty to avoid confounding with random intercept

results = ptk.run_comprehensive_statistical_analysis(
    normalized_data=normalized_data,
    sample_metadata=sample_metadata,
    config=config
)

# Linear trend analysis (continuous time/dose)
config.analysis_type = "linear_trend"  # Test if slope != 0
results = ptk.run_comprehensive_statistical_analysis(
    normalized_data=normalized_data,
    sample_metadata=sample_metadata,
    config=config
)
```

## Installation

```bash
# Install in development mode
pip install -e /path/to/proteomics_toolkit
```

## Dependencies
- pandas >= 1.3.0
- numpy >= 1.21.0
- scipy >= 1.7.0
- matplotlib >= 3.4.0
- seaborn >= 0.11.0
- scikit-learn >= 1.0.0
- statsmodels >= 0.12.0
- requests >= 2.25.0 (for Enrichr API in temporal_clustering)

## Module Overview

### data_import.py
- `load_skyline_data()`: Load protein/peptide quantitation and metadata
- `identify_sample_columns()`: Automatically detect sample columns
- `clean_sample_names()`: Remove common prefixes/suffixes
- `match_samples_to_metadata()`: Link samples to experimental metadata
- `identify_and_classify_controls()`: Intelligent control sample detection

### preprocessing.py  
- `parse_protein_identifiers()`: Extract UniProt accessions and databases
- `parse_gene_and_description()`: Parse gene names and clean descriptions
- `assess_data_completeness()`: Evaluate missing data patterns
- `filter_proteins_by_completeness()`: Remove proteins with too much missing data
- `calculate_group_colors()`: Generate consistent colors for visualization

### normalization.py
- `median_normalize()`: Median-based normalization
- `vsn_normalize()`: Variance Stabilizing Normalization
- `quantile_normalize()`: Quantile normalization
- `calculate_normalization_stats()`: Evaluate normalization effectiveness

### differential.py
- `run_differential_analysis()`: Statistical testing for differential expression
- `calculate_effect_sizes()`: Cohen's d and effect size categories
- `export_results()`: Save results in standardized format

### visualization.py
- `plot_box_plot()`: Sample intensity distributions
- `plot_pca()`: Principal component analysis
- `plot_sample_correlation_heatmap()`: Sample-sample correlation matrix
- `plot_volcano()`: Volcano plots for differential results
- `plot_normalization_comparison()`: Before/after normalization comparison
- `plot_control_cv_distribution()`: CV distribution analysis for control samples
- `plot_grouped_heatmap()`: Heatmap for any grouped data (clusters, treatments, doses)
- `plot_grouped_trajectories()`: Line plots for grouped trajectories (temporal, dose-response)
- `plot_protein_profile()`: Single protein expression profile across conditions

### enrichment.py
General-purpose gene set enrichment analysis via the Enrichr API.

**Configuration:**
- `EnrichmentConfig`: Dataclass for configuring libraries, thresholds, and API settings

**Core Functions:**
- `query_enrichr()`: Query the Enrichr API directly with a gene list
- `parse_enrichr_results()`: Parse raw API results into a tidy DataFrame
- `run_enrichment_analysis()`: Complete enrichment analysis on a gene list
- `run_enrichment_by_group()`: Run enrichment for each group (clusters, categories, etc.)
- `run_differential_enrichment()`: Run enrichment on up/down-regulated genes from differential analysis

**Visualization:**
- `plot_enrichment_barplot()`: Horizontal bar plots of enrichment results
- `plot_enrichment_comparison()`: Dot plots comparing enrichment across groups

**Utilities:**
- `get_available_libraries()`: List commonly used Enrichr libraries
- `merge_enrichment_results()`: Merge multiple enrichment DataFrames

**Example Usage:**
```python
import proteomics_toolkit as ptk

# Enrichment on a gene list (any source)
gene_list = ['TP53', 'BRCA1', 'BRCA2', 'ATM', 'CHEK2']
enrichment = ptk.run_enrichment_analysis(gene_list)

# Enrichment on differential expression results
enrichment_by_direction = ptk.run_differential_enrichment(
    stats_results,
    logfc_threshold=1.0,
    pvalue_threshold=0.05
)

# Enrichment by cluster or group
enrichment_by_cluster = ptk.run_enrichment_by_group(
    clustered_df,
    group_column='Cluster_Name',
    gene_column='Gene'
)

# Visualize results
fig = ptk.plot_enrichment_barplot(enrichment, title='Top Pathways')
fig = ptk.plot_enrichment_comparison(enrichment_by_cluster, title='Cluster Comparison')
```

### temporal_clustering.py
Comprehensive module for analyzing temporal/longitudinal proteomics data.

**Configuration:**
- `TemporalClusteringConfig`: Dataclass for configuring clustering, significance, and visualization settings

**Core Analysis Functions:**
- `calculate_temporal_means()`: Calculate mean protein abundance at each timepoint across subjects
- `cluster_temporal_trends()`: K-means or hierarchical clustering of temporal trajectories
- `name_clusters_by_pattern()`: Assign descriptive names (e.g., "Sustained Increase", "Early Response") to clusters
- `classify_trend_pattern()`: Classify individual protein temporal patterns
- `merge_with_statistics()`: Combine temporal data with differential expression results
- `filter_significant_proteins()`: Filter to statistically significant proteins only

**Gene Set Enrichment (via Enrichr API):**
- `query_enrichr()`: Query the Enrichr API for gene set enrichment analysis
- `run_enrichment_by_cluster()`: Run enrichment analysis for each cluster

**Visualization:**
- `plot_cluster_heatmap()`: Heatmap of protein expression organized by cluster
- `plot_cluster_parallel_coordinates()`: Parallel coordinate plots for temporal patterns
- `plot_enrichment_barplot()`: Horizontal bar plots of enrichment results
- `plot_enrichment_comparison()`: Dot plots comparing enrichment across clusters

**Complete Pipeline:**
- `run_temporal_analysis()`: Complete temporal analysis pipeline including clustering, visualization, and enrichment

**Example Usage:**
```python
import proteomics_toolkit as ptk

# Configure temporal analysis
config = ptk.TemporalClusteringConfig(
    n_clusters=4,
    p_value_threshold=0.05,
    linewidth=2.0,
    alpha=0.4
)

# Run complete analysis
results = ptk.run_temporal_analysis(
    data_df=normalized_data,
    metadata_dict=sample_metadata_dict,
    stat_results=statistical_results,
    treatment_name='Verapamil',
    config=config,
    run_enrichment=True
)

# Access individual results
print(f"Significant proteins: {len(results['sig_df'])}")
fig_heatmap = results['fig_heatmap']
fig_parallel = results['fig_parallel']
enrichment = results.get('enrichment_results', {})
```

## Contributing

This toolkit is designed to be modular and extensible. To add new functionality:

1. Follow existing naming conventions
2. Include comprehensive docstrings
3. Add input validation and error handling
4. Provide usage examples
5. Update relevant documentation

## License

MIT License - see LICENSE file for details.
