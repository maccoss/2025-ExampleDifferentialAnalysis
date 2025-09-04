# Proteomics Analysis Toolkit and Example Workflows

This repository provides a **proteomics analysis toolkit** with example workflows for analyzing quantitative proteomics data from Skyline. The toolkit supports end-to-end analysis from data import through statistical analysis and visualization.

## Getting Started - Minimal Analysis Notebook

**`Minimal-Analysis.ipynb`** - A simple notebook that demonstrates the analysis workflow in a single python notebook.

### Key Features:
- **Complete workflow** in one notebook with minimal code
- **Comprehensive configuration** - all analysis parameters in one place
- **All visualization capabilities** - quality control, normalization comparison, PCA, correlation analysis
- **Statistical analysis** - mixed-effects models, t-tests, and non-parametric tests
- **Flexible design** - works with dose-response, time-course, or treatment comparison studies
- **Export capabilities** - saves results, configuration, and processed data
- **Automated Tests** - Use pytest for testing

### Quick Start:
1. Configure your analysis parameters in the first cell
2. Run all cells to perform complete analysis
3. Review plots and tables for quality control and results
4. Export results for downstream analysis or reporting

---

## Data Input Requirements

The toolkit uses two standard Skyline reports:
- **Replicates report** - Export metadata annotations from the document grid
- **Protein quantitation report** - Use the included `MJM Protein Total Areas.skyr` template
- **Sample metadata CSV** - Contains experimental design information

---

## Proteomics Toolkit Capabilities

The `proteomics_toolkit` package provides modular functions for comprehensive proteomics analysis:

### **Data Import & Processing**
- **Skyline output parsing** - Handles protein quantitation matrices and metadata
- **Protein identifier extraction** - Robust UniProt accession parsing with regex
- **Sample name normalization** - Automatic cleanup of sample naming conventions
- **Data structure standardization** - Creates consistent data formats for analysis

### **Quality Assessment**
- **Missing data visualization** - Comprehensive plots of protein detection patterns
- **Sample completeness metrics** - Group-wise statistics and quality reporting
- **Control sample analysis** - Automated QC pool correlation analysis
- **Data integrity validation** - Ensures proper experimental design structure

### **Normalization Methods**
- **Median normalization** - Simple, robust correction for loading differences
- **VSN (Variance Stabilizing Normalization)** - Reduces heteroscedasticity
- **Quantile normalization** - Makes sample distributions identical
- **MAD normalization** - Median Absolute Deviation (robust to outliers)
- **Z-score normalization** - Standardizes to mean=0, std=1
- **RLR (Robust Linear Regression)** - Corrects for systematic batch effects
- **LOESS normalization** - Handles intensity-dependent bias
- **Negative value handling** - Comprehensive strategies for VSN output

### **Statistical Analysis**
- **Mixed-effects models** - Handles complex experimental designs with random effects
- **T-tests** - Paired and unpaired comparisons with multiple correction
- **Non-parametric tests** - Mann-Whitney U and Wilcoxon signed-rank tests
- **Covariate adjustment** - Incorporates additional metadata variables
- **Multiple testing correction** - Benjamini-Hochberg FDR and other methods
- **Effect size calculation** - Fold changes and confidence intervals

### **Visualization Suite**
- **Box plots** - Sample distribution analysis with group coloring
- **PCA analysis** - Principal component plots with normalization comparison
- **Correlation heatmaps** - Triangular matrices with hierarchical clustering
- **Volcano plots** - Statistical significance vs. fold change visualization
- **Control sample QC** - Specialized plots for quality assessment
- **Normalization comparison** - Before/after distribution analysis

### **Export & Reporting**
- **Results export** - Statistical tables with comprehensive annotations
- **Configuration saving** - Complete analysis parameter documentation
- **Publication-ready tables** - Formatted output for manuscripts
- **Processed data export** - Normalized data matrices for downstream analysis

---

## Example Notebooks (more to come)

### **Minimal Analysis**
- **`Minimal-Analysis.ipynb`** - Complete workflow

---

## Installation and Setup

```python
# Clone the repository
git clone https://github.com/maccoss/2025-ExampleDifferentialAnalysis.git
cd 2025-ExampleDifferentialAnalysis

# Install required packages
pip install -r requirements.txt

# Start with the minimal analysis notebook
jupyter notebook EISAI-Minimal-Analysis.ipynb
```

## Analysis Workflow

### 1. Configure Analysis Parameters
Set all analysis parameters in the comprehensive configuration cell:
- Input file paths and data filtering options
- Normalization method selection (8 options available)
- Statistical test method and experimental design
- Significance thresholds and visualization settings

### 2. Data Processing Pipeline
- **Import**: Load Skyline outputs and metadata
- **Clean**: Normalize sample names and parse protein identifiers  
- **Filter**: Remove low-quality proteins based on detection rate
- **Visualize**: Quality control plots for raw data assessment

### 3. Normalization and Quality Control
- **Normalize**: Apply selected normalization method
- **Compare**: Before/after distribution analysis
- **Assess**: PCA, correlation, and control sample analysis
- **Validate**: Statistical metrics for normalization effectiveness

### 4. Statistical Analysis
- **Configure**: Set up statistical model parameters
- **Analyze**: Run comprehensive differential analysis
- **Visualize**: Create volcano plots and results tables
- **Export**: Save results and processed data

## Key Advantages

- **Streamlined Workflow**: Complete analysis in one notebook
- **Comprehensive Configuration**: All parameters documented and configurable
- **Production Ready**: Robust error handling and quality control
- **Flexible Design**: Works with various experimental designs
- **Reproducible**: Complete configuration export for method documentation
- **Publication Ready**: High-quality plots and formatted results tables

---

## Repository Structure

```
├── EISAI-Minimal-Analysis.ipynb          # RECOMMENDED: Complete streamlined workflow
├── proteomics_toolkit/                   # Core analysis functions
│   ├── data_import.py                   # Data loading and preprocessing
│   ├── normalization.py                # All normalization methods
│   ├── statistical_analysis.py         # Statistical tests and modeling
│   ├── visualization.py                # Plotting functions
│   ├── preprocessing.py                # Data filtering and parsing
│   └── export.py                       # Results export and configuration
├── example_data/                        # Sample datasets for testing
├── old_notebooks/                       # Legacy detailed step-by-step notebooks
├── tests/                               # Comprehensive test suite
└── requirements.txt                     # Python dependencies
```

## Example Data

The repository includes example CSF proteomics data from an pilot study demonstrating:
- **Dose-response analysis**: 4 dose levels (0, 20, 40, 80 mg)
- **Longitudinal design**: Paired samples at baseline (D-02) and follow-up (D-13)
- **Quality control samples**: Multiple control pools for normalization assessment
- **Mixed-effects modeling**: Complex experimental design with subject-level random effects

---

*This toolkit provides a comprehensive, production-ready solution for quantitative proteomics analysis from Skyline outputs. The minimal notebook approach ensures reproducible, well-documented analyses suitable for publication and regulatory submissions.*


