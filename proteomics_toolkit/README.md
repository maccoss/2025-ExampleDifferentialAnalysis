# Proteomics Analysis Toolkit

A comprehensive Python toolkit for analyzing mass spectrometry-based proteomics data, with a focus on Skyline quantitation workflows.

## Features

### 🔬 Core Analysis Modules
- **data_import**: Load Skyline data, parse protein identifiers, and manage sample metadata
- **preprocessing**: Data quality assessment, filtering, and annotation parsing  
- **normalization**: Multiple normalization methods (median, VSN, quantile)
- **differential**: Statistical analysis for differential protein expression
- **visualization**: Professional plots for data interpretation

### 🎯 Key Capabilities
- **Automated Control Detection**: Intelligent identification of QC, pool, and reference samples ([guide](README_control_identification.md))
- **Flexible Data Import**: Support for Skyline protein/peptide quantitation files
- **Multi-method Normalization**: Choose the best normalization for your data
- **Statistical Analysis**: Paired/unpaired differential expression with effect sizes  
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

## Contributing

This toolkit is designed to be modular and extensible. To add new functionality:

1. Follow existing naming conventions
2. Include comprehensive docstrings
3. Add input validation and error handling
4. Provide usage examples
5. Update relevant documentation

## License

MIT License - see LICENSE file for details.
