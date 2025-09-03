# =============================================================================
# PROTEOMICS ANALYSIS CONFIGURATION
# Generated: 2025-09-03 10:27:14
# Analysis: Minimal dose-response analysis with comprehensive configuration
# =============================================================================

# =============================================================================
# 1. INPUT FILES AND PATHS
# =============================================================================
toolkit_path = '.'
metadata_file = 'example_data/2025-6-CSF-Total-Pilot-MetaData.csv'
protein_file = 'example_data/2025-6-CSF-Total-Pilot-ProteinQuant.csv'
peptide_file = ''
remove_common_prefix = True

# =============================================================================
# 2. DATA FILTERING PARAMETERS
# =============================================================================
min_detection_rate = 0.5

# =============================================================================
# 3. NORMALIZATION STRATEGY
# =============================================================================
normalization_method = 'Median'
optimize_vsn = False

# =============================================================================
# 4. NEGATIVE VALUE HANDLING STRATEGY
# =============================================================================
handle_negatives = False
negative_handling_method = 'min_positive'
min_positive_replacement = None

# =============================================================================
# 5. STATISTICAL ANALYSIS STRATEGY
# =============================================================================
statistical_test_method = 'mixed_effects'
analysis_type = 'paired'

# =============================================================================
# 6. EXPERIMENTAL DESIGN CONFIGURATION
# =============================================================================
subject_column = 'Subject'
paired_column = 'Visit'
paired_label1 = 'D-02'
paired_label2 = 'D-13'
group_column = 'DrugDose'
group_labels = ['0', '20', '40', '80']

# =============================================================================
# 7. MIXED-EFFECTS MODEL CONFIGURATION
# =============================================================================
interaction_terms = ['DrugDose', 'Visit']
additional_interactions = []
covariates = []

# =============================================================================
# 8. CONTROL SAMPLE CONFIGURATION
# =============================================================================
control_column = 'Subject'
control_labels = ['HoofPool', 'GWPool', 'EISAIPool']

# =============================================================================
# 9. VISUALIZATION SETTINGS
# =============================================================================
use_systematic_colors = True
systematic_color_palette = 'Set1'
group_order = None
group_colors = None

# =============================================================================
# 10. SIGNIFICANCE THRESHOLDS
# =============================================================================
p_value_threshold = 0.05
fold_change_threshold = 1.5
q_value_max = 0.1
use_adjusted_pvalue = 'adjusted'
enable_pvalue_fallback = True

# =============================================================================
# 11. OUTPUT AND EXPORT SETTINGS
# =============================================================================
export_results = True
output_prefix = 'EISAI-Minimal-Analysis'
label_top_proteins = 25
random_seed = 42
min_samples_per_group = 3

# =============================================================================
# COMPUTED VALUES (for reference)
# =============================================================================
# group_colors: None
# Total proteins analyzed: 1584
# Total samples: 48
