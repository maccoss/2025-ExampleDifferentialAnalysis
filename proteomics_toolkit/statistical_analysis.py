"""
Statistical Analysis Module for Proteomics Data

This module provides a clean, configuration-driven approach to differential
protein expression analysis with support for various statistical methods.
"""

import pandas as pd
import numpy as np
from scipy.stats import ttest_1samp, ttest_ind, mannwhitneyu, wilcoxon
from statsmodels.stats.multitest import multipletests
import warnings
from .normalization import is_normalization_log_transformed

# Try to import statsmodels for mixed-effects models
try:
    from statsmodels.formula.api import mixedlm

    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False
    mixedlm = None
    warnings.warn("statsmodels not available. Mixed-effects models will be disabled.")


def _apply_log_transformation_if_needed(data, config):
    """
    Apply log transformation to data if needed based on configuration.
    Uses existing normalization infrastructure to determine if data is already log-transformed.

    Parameters:
    -----------
    data : pd.DataFrame
        Input data with numeric columns to potentially transform
    config : StatisticalConfig
        Configuration object containing log transformation settings

    Returns:
    --------
    pd.DataFrame
        Data with log transformation applied if needed
    """
    # Get numeric columns (sample data)
    numeric_columns = data.select_dtypes(include=[np.number]).columns.tolist()

    if config.log_transform_before_stats == "auto":
        # Use existing normalization infrastructure to determine if data is already log-transformed
        if hasattr(config, "normalization_method") and config.normalization_method:
            already_log_transformed = is_normalization_log_transformed(
                config.normalization_method
            )

            if already_log_transformed:
                apply_log_transform = False
                print(
                    f"Log transformation: AUTO-DETECTED (not needed - {config.normalization_method} already log-transforms data)"
                )
            else:
                apply_log_transform = True
                print(
                    f"Log transformation: AUTO-DETECTED (needed - {config.normalization_method} preserves original scale)"
                )
        else:
            # No normalization method info, check data range as fallback
            if numeric_columns:
                sample_data_range = data[numeric_columns]
                mean_value = sample_data_range.mean().mean()
                apply_log_transform = mean_value > 50
                status = "needed" if apply_log_transform else "not needed"
                print(
                    f"Log transformation: AUTO-DETECTED ({status} - mean value {mean_value:.1f})"
                )
            else:
                apply_log_transform = False
                print("Log transformation: AUTO-DETECTED (no numeric columns found)")

    elif str(config.log_transform_before_stats).lower() in ["true", "1", "yes", "on"]:
        apply_log_transform = True
        print("Log transformation: ENABLED (forced by configuration)")
    else:
        apply_log_transform = False
        print("Log transformation: DISABLED (by configuration)")

    if not apply_log_transform or not numeric_columns:
        print("Using data as-is for statistical analysis")
        return data

    # Apply log transformation
    print(f"Applying {config.log_base} transformation for statistical analysis...")

    # Create a copy to avoid modifying original data
    transformed_data = data.copy()
    sample_data_subset = transformed_data[numeric_columns]

    # Handle negative values if any exist
    if (sample_data_subset < 0).any().any():
        print("  -> Handling negative values...")
        min_val = sample_data_subset.min().min()
        shift_amount = abs(min_val) + 1
        transformed_data[numeric_columns] = (
            transformed_data[numeric_columns] + shift_amount
        )
        print(f"     Shifted all values by +{shift_amount:.2f}")

    # Determine pseudocount
    if config.log_pseudocount is None:
        pseudocount = (
            max(1e-6, sample_data_subset.min().min() / 100)
            if sample_data_subset.min().min() > 0
            else 0.1
        )
    else:
        pseudocount = config.log_pseudocount

    # Apply appropriate log transformation
    if config.log_base == "log2":
        transformed_data[numeric_columns] = np.log2(
            transformed_data[numeric_columns] + pseudocount
        )
    elif config.log_base == "log10":
        transformed_data[numeric_columns] = np.log10(
            transformed_data[numeric_columns] + pseudocount
        )
    elif config.log_base == "ln":
        transformed_data[numeric_columns] = np.log(
            transformed_data[numeric_columns] + pseudocount
        )
    else:
        raise ValueError(f"Unknown log base: {config.log_base}")

    print(
        f"  -> Applied {config.log_base} transformation with pseudocount {pseudocount}"
    )

    # Verify transformation
    new_mean = transformed_data[numeric_columns].mean().mean()
    print(
        f"  -> New data range: {transformed_data[numeric_columns].min().min():.2f} to {transformed_data[numeric_columns].max().max():.2f}"
    )
    print(f"  -> New mean: {new_mean:.2f}")

    return transformed_data


class StatisticalConfig:
    """Configuration class for statistical analysis parameters"""

    def __init__(self):
        # Basic analysis parameters
        self.statistical_test_method = "mixed_effects"
        self.analysis_type = "paired"
        self.p_value_threshold = 0.05
        self.fold_change_threshold = 1.5

        # Experimental design
        self.subject_column = "Subject"
        self.paired_column = "Visit"
        self.paired_label1 = "D-02"
        self.paired_label2 = "D-13"
        self.group_column = "Comparison"
        self.group_labels = ["Placebo", "Drug"]

        # Mixed-effects specific
        self.interaction_terms = ["Comparison", "Visit"]
        self.additional_interactions = []
        self.covariates = []

        # Multiple testing correction
        self.correction_method = "fdr_bh"

        # P-value selection parameters
        self.use_adjusted_pvalue = "adjusted"  # "adjusted" or "unadjusted"
        self.enable_pvalue_fallback = (
            True  # Auto-fallback to unadjusted if no adjusted significant
        )

        # Log transformation parameters
        self.log_transform_before_stats = "auto"  # "auto", True, False
        self.log_base = "log2"  # "log2", "log10", "ln"
        self.log_pseudocount = None  # None for auto, or specific value

        # Normalization method (used for auto log transformation)
        self.normalization_method = None  # Set this to the normalization method used


def prepare_metadata_dataframe(sample_metadata_dict, sample_columns, config):
    """Convert sample metadata dictionary to DataFrame suitable for analysis"""

    print(f"Preparing metadata for {len(sample_columns)} samples...")

    # Create DataFrame from metadata dictionary
    metadata_rows = []
    for sample_name in sample_columns:
        if sample_name in sample_metadata_dict:
            row = sample_metadata_dict[sample_name].copy()
            row["Sample"] = sample_name
            metadata_rows.append(row)
        else:
            print(f"Warning: No metadata found for sample {sample_name}")

    if not metadata_rows:
        print("Warning: No metadata found for any samples - returning empty DataFrame")
        # Return empty DataFrame with expected columns for graceful handling
        expected_cols = ["Sample", config.subject_column, config.paired_column, config.group_column]
        return pd.DataFrame(columns=expected_cols)

    metadata_df = pd.DataFrame(metadata_rows)

    # Skip validation if no data (graceful handling for edge cases)
    if len(metadata_df) == 0:
        print("Empty metadata - skipping validation")
        return metadata_df

    # Validate required columns exist
    required_cols = [config.subject_column, config.paired_column, config.group_column]
    missing_cols = [col for col in required_cols if col not in metadata_df.columns]
    if missing_cols:
        raise ValueError(f"Missing required metadata columns: {missing_cols}")

    # CRITICAL: Filter out samples with missing values in required columns
    print(f"  Before filtering: {len(metadata_df)} samples")
    for col in required_cols:
        before_count = len(metadata_df)
        metadata_df = metadata_df.dropna(subset=[col])
        after_count = len(metadata_df)
        if before_count != after_count:
            print(f"  Removed {before_count - after_count} samples missing {col}")

    if len(metadata_df) == 0:
        raise ValueError("No samples remain after filtering for required metadata")

    print(f"  After filtering: {len(metadata_df)} samples")
    print(f"  Subjects: {metadata_df[config.subject_column].nunique()}")
    print(f"  Groups: {metadata_df[config.group_column].value_counts().to_dict()}")
    print(f"  Timepoints: {metadata_df[config.paired_column].value_counts().to_dict()}")

    return metadata_df


def run_paired_t_test(protein_data, metadata_df, config):
    """Run paired t-test analysis"""

    print("Running paired t-test analysis...")

    results = []
    n_proteins = len(protein_data)

    for i, (protein_idx, protein_values) in enumerate(protein_data.iterrows()):
        if (i + 1) % 200 == 0:
            print(f"  Processed {i + 1}/{n_proteins} proteins...")

        # Get data for this protein
        protein_df = pd.DataFrame(
            {"Sample": protein_values.index, "Intensity": protein_values.values}
        )

        # Merge with metadata
        protein_df = protein_df.merge(metadata_df, on="Sample", how="inner")

        # Remove missing values
        protein_df = protein_df.dropna(subset=["Intensity"])

        if len(protein_df) < 4:  # Need at least some data
            results.append(_create_empty_result(protein_idx, "Insufficient data"))
            continue

        # Calculate paired differences for each subject
        baseline_data = protein_df[
            protein_df[config.paired_column] == config.paired_label1
        ]
        followup_data = protein_df[
            protein_df[config.paired_column] == config.paired_label2
        ]

        # Merge on subject to get paired data
        paired_data = baseline_data.merge(
            followup_data, on=config.subject_column, suffixes=("_baseline", "_followup")
        )

        if len(paired_data) < 3:  # Need at least 3 pairs
            results.append(
                _create_empty_result(protein_idx, "Insufficient paired data")
            )
            continue

        # Calculate differences (followup - baseline)
        differences = (
            paired_data["Intensity_followup"] - paired_data["Intensity_baseline"]
        )

        try:
            # Paired t-test (test if mean difference != 0)
            t_stat, p_value = ttest_1samp(differences, 0)

            # Calculate effect size (Cohen's d for paired data)
            mean_diff = differences.mean()
            std_diff = differences.std()
            cohens_d = mean_diff / std_diff if std_diff > 0 else 0

            # Log fold change (approximate)
            log_fc = mean_diff  # Already in log-like space if VSN normalized

            result = {
                "Protein": protein_idx,
                "logFC": log_fc,
                "AveExpr": protein_df["Intensity"].mean(),
                "t": t_stat,
                "P.Value": p_value,
                "B": np.nan,  # Not applicable for t-test
                "n_pairs": len(paired_data),
                "mean_diff": mean_diff,
                "std_diff": std_diff,
                "cohens_d": cohens_d,
                "test_method": "Paired t-test",
            }

            results.append(result)

        except Exception as e:
            results.append(_create_empty_result(protein_idx, f"Analysis failed: {e}"))

    print(f"✓ Paired t-test completed for {len(results)} proteins")
    return pd.DataFrame(results)


def run_mixed_effects_analysis(
    protein_data, metadata_df, config, protein_annotations=None
):
    """Run mixed-effects model analysis"""

    if not HAS_STATSMODELS:
        raise ImportError("statsmodels is required for mixed-effects analysis")

    print("Running mixed-effects analysis...")
    print(
        f"  Model: Protein ~ {' * '.join(config.interaction_terms)} + (1|{config.subject_column})"
    )

    results = []
    n_proteins = len(protein_data)

    for i, (protein_idx, protein_values) in enumerate(protein_data.iterrows()):
        if (i + 1) % 100 == 0:
            print(f"  Processed {i + 1}/{n_proteins} proteins...")

        # Get actual protein name if annotations are provided
        if protein_annotations is not None and "Protein" in protein_annotations.columns:
            # Use protein_idx (the actual row index) to get the correct protein name
            actual_protein_name = protein_annotations.loc[protein_idx, "Protein"]
        else:
            actual_protein_name = protein_idx  # Fallback to index

        # Prepare data for this protein
        protein_df = pd.DataFrame(
            {"Sample": protein_values.index, "Intensity": protein_values.values}
        )

        # Merge with metadata
        protein_df = protein_df.merge(metadata_df, on="Sample", how="inner")

        # Remove missing values
        protein_df = protein_df.dropna(subset=["Intensity"])

        if len(protein_df) < 8:  # Need sufficient data for mixed model
            results.append(
                _create_empty_mixed_effects_result(
                    actual_protein_name, "Insufficient data"
                )
            )
            continue

        try:
            # Build formula with all interaction terms and additional interactions
            all_interaction_terms = (
                config.interaction_terms + config.additional_interactions
            )

            if len(all_interaction_terms) < 2:
                raise ValueError("Need at least 2 terms for interaction analysis")

            # Build main interaction (first two terms)
            formula = (
                f"Intensity ~ {all_interaction_terms[0]} * {all_interaction_terms[1]}"
            )

            # Add additional interaction terms as main effects
            if len(all_interaction_terms) > 2:
                additional_terms = " + ".join(all_interaction_terms[2:])
                formula += f" + {additional_terms}"

            # Add covariates if specified
            if config.covariates:
                formula += " + " + " + ".join(config.covariates)

            # Fit mixed-effects model
            if mixedlm is None:
                raise ImportError("statsmodels required for mixed-effects analysis")

            # Suppress convergence warnings during fitting
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=UserWarning)
                warnings.filterwarnings("ignore", category=RuntimeWarning)
                warnings.filterwarnings("ignore", message=".*convergence.*")
                warnings.filterwarnings("ignore", message=".*singular.*")

                model = mixedlm(
                    formula, protein_df, groups=protein_df[config.subject_column]
                )
                fitted_model = model.fit(method="lbfgs")

            # Extract results
            params = fitted_model.params
            pvalues = fitted_model.pvalues

            # Get interaction effect (main result of interest)
            # Look for any interaction term (statsmodels creates them based on alphabetical order)
            interaction_coef = np.nan
            interaction_pvalue = np.nan

            # Find interaction parameters (contain both variable names and ":")
            interaction_candidates = [
                p
                for p in params.index
                if config.interaction_terms[0] in p
                and config.interaction_terms[1] in p
                and ":" in p
            ]

            if interaction_candidates:
                # Use the first (and typically only) interaction term
                term_name = interaction_candidates[0]
                interaction_coef = params[term_name]
                interaction_pvalue = pvalues[term_name]
            else:
                # No interaction parameter found - will use NaN values
                pass

            # Get main effects with flexible pattern matching
            group_effect = np.nan
            group_pvalue = np.nan

            # Find group effect parameters
            group_candidates = [
                p
                for p in params.index
                if config.interaction_terms[0] in p
                and ":" not in p
                and p != "Intercept"
            ]

            if group_candidates:
                term_name = group_candidates[0]
                group_effect = params[term_name]
                group_pvalue = pvalues[term_name]

            # Find time effect parameters
            time_effect = np.nan
            time_pvalue = np.nan

            time_candidates = [
                p
                for p in params.index
                if config.interaction_terms[1] in p
                and ":" not in p
                and p != "Intercept"
            ]

            if time_candidates:
                term_name = time_candidates[0]
                time_effect = params[term_name]
                time_pvalue = pvalues[term_name]

            result = {
                "Protein": actual_protein_name,  # Use actual protein name
                "logFC": interaction_coef,  # Use interaction as primary effect
                "AveExpr": protein_df["Intensity"].mean(),
                "t": np.nan,  # Not applicable for mixed model
                "P.Value": interaction_pvalue,  # Use interaction p-value as primary
                "B": np.nan,  # Not applicable
                "group_effect": group_effect,
                "group_pvalue": group_pvalue,
                "time_effect": time_effect,
                "time_pvalue": time_pvalue,
                "aic": fitted_model.aic,
                "bic": fitted_model.bic,
                "n_obs": len(protein_df),
                "test_method": "Mixed-effects model",
            }

            results.append(result)

        except Exception as e:
            error_message = f"Model failed: {e}"
            if i < 3:  # Show details for first few failures
                print(f"  Protein {i + 1} failed: {error_message}")
            results.append(
                _create_empty_mixed_effects_result(actual_protein_name, error_message)
            )

    print(f"✓ Mixed-effects analysis completed for {len(results)} proteins")
    return pd.DataFrame(results)


def run_unpaired_t_test(protein_data, metadata_df, config):
    """Run unpaired t-test analysis"""

    print("Running unpaired t-test analysis...")

    results = []
    n_proteins = len(protein_data)

    # Filter to specific timepoint if needed
    if config.paired_column and config.paired_label2:
        metadata_df = metadata_df[
            metadata_df[config.paired_column] == config.paired_label2
        ]
        print(f"  Analyzing {config.paired_label2} timepoint only")

    for i, (protein_idx, protein_values) in enumerate(protein_data.iterrows()):
        if (i + 1) % 200 == 0:
            print(f"  Processed {i + 1}/{n_proteins} proteins...")

        # Get data for this protein
        protein_df = pd.DataFrame(
            {"Sample": protein_values.index, "Intensity": protein_values.values}
        )

        # Merge with metadata
        protein_df = protein_df.merge(metadata_df, on="Sample", how="inner")
        protein_df = protein_df.dropna(subset=["Intensity"])

        if len(protein_df) < 4:
            results.append(_create_empty_result(protein_idx, "Insufficient data"))
            continue

        # Split into groups
        group1_data = protein_df[
            protein_df[config.group_column] == config.group_labels[0]
        ]["Intensity"]
        group2_data = protein_df[
            protein_df[config.group_column] == config.group_labels[1]
        ]["Intensity"]

        if len(group1_data) < 2 or len(group2_data) < 2:
            results.append(_create_empty_result(protein_idx, "Insufficient group data"))
            continue

        try:
            # Welch's t-test (unequal variances)
            t_stat, p_value = ttest_ind(group2_data, group1_data, equal_var=False)

            # Calculate effect size
            pooled_std = np.sqrt(
                (
                    (len(group1_data) - 1) * group1_data.var()
                    + (len(group2_data) - 1) * group2_data.var()
                )
                / (len(group1_data) + len(group2_data) - 2)
            )
            cohens_d = (
                (group2_data.mean() - group1_data.mean()) / pooled_std
                if pooled_std > 0
                else 0
            )

            # Log fold change
            log_fc = group2_data.mean() - group1_data.mean()

            result = {
                "Protein": protein_idx,
                "logFC": log_fc,
                "AveExpr": protein_df["Intensity"].mean(),
                "t": t_stat,
                "P.Value": p_value,
                "B": np.nan,
                "n_group1": len(group1_data),
                "n_group2": len(group2_data),
                "cohens_d": cohens_d,
                "test_method": "Welch t-test",
            }

            results.append(result)

        except Exception as e:
            results.append(_create_empty_result(protein_idx, f"Analysis failed: {e}"))

    print(f"✓ Unpaired t-test completed for {len(results)} proteins")
    return pd.DataFrame(results)


def run_wilcoxon_test(protein_data, metadata_df, config):
    """Run paired non-parametric Wilcoxon signed-rank test"""

    print("Running Wilcoxon signed-rank test analysis...")

    results = []
    n_proteins = len(protein_data)

    for i, (protein_idx, protein_values) in enumerate(protein_data.iterrows()):
        if (i + 1) % 200 == 0:
            print(f"  Processed {i + 1}/{n_proteins} proteins...")

        # Get data for this protein
        protein_df = pd.DataFrame(
            {"Sample": protein_values.index, "Intensity": protein_values.values}
        )

        # Merge with metadata
        protein_df = protein_df.merge(metadata_df, on="Sample", how="inner")

        # Remove missing values
        protein_df = protein_df.dropna(subset=["Intensity"])

        if len(protein_df) < 4:  # Need at least some data
            results.append(_create_empty_result(protein_idx, "Insufficient data"))
            continue

        # Calculate paired differences for each subject
        baseline_data = protein_df[
            protein_df[config.paired_column] == config.paired_label1
        ]
        followup_data = protein_df[
            protein_df[config.paired_column] == config.paired_label2
        ]

        # Merge on subject to get paired data
        paired_data = baseline_data.merge(
            followup_data, on=config.subject_column, suffixes=("_baseline", "_followup")
        )

        if len(paired_data) < 3:  # Need at least 3 pairs
            results.append(
                _create_empty_result(protein_idx, "Insufficient paired data")
            )
            continue

        # Calculate differences (followup - baseline)
        differences = (
            paired_data["Intensity_followup"] - paired_data["Intensity_baseline"]
        )

        # Remove zero differences for Wilcoxon test
        non_zero_diffs = differences[differences != 0]
        
        if len(non_zero_diffs) < 3:
            results.append(
                _create_empty_result(protein_idx, "Insufficient non-zero differences")
            )
            continue

        try:
            # Wilcoxon signed-rank test
            statistic, p_value = wilcoxon(non_zero_diffs, alternative='two-sided')

            # Calculate effect size (r = z / sqrt(N))
            # For Wilcoxon, we use median and IQR
            mean_diff = differences.mean()
            median_diff = differences.median()

            # Pseudo-Cohen's d using median and MAD (more robust)
            mad = np.median(np.abs(differences - median_diff)) * 1.4826  # Scale factor for normality
            effect_size = median_diff / mad if mad > 0 else 0

            result = {
                "Protein": protein_idx,
                "logFC": median_diff,  # Use median for non-parametric
                "AveExpr": protein_df["Intensity"].mean(),
                "statistic": statistic,
                "P.Value": p_value,
                "B": np.nan,  # Not applicable for non-parametric tests
                "n_pairs": len(paired_data),
                "mean_diff": mean_diff,
                "median_diff": median_diff,
                "Effect_Size": effect_size,
                "test_method": "Wilcoxon signed-rank",
            }

            results.append(result)

        except Exception as e:
            results.append(_create_empty_result(protein_idx, f"Analysis failed: {e}"))

    print(f"✓ Wilcoxon signed-rank test completed for {len(results)} proteins")
    return pd.DataFrame(results)


def run_mann_whitney_test(protein_data, metadata_df, config):
    """Run unpaired non-parametric Mann-Whitney U test"""

    print("Running Mann-Whitney U test analysis...")

    results = []
    n_proteins = len(protein_data)

    # Filter to specific timepoint if needed (only if paired_column is present in metadata)
    if (config.paired_column and 
        config.paired_label2 and 
        hasattr(config, 'paired_column') and 
        config.paired_column in metadata_df.columns):
        metadata_df = metadata_df[
            metadata_df[config.paired_column] == config.paired_label2
        ]
        print(f"  Analyzing {config.paired_label2} timepoint only")

    for i, (protein_idx, protein_values) in enumerate(protein_data.iterrows()):
        if (i + 1) % 200 == 0:
            print(f"  Processed {i + 1}/{n_proteins} proteins...")

        # Get data for this protein
        protein_df = pd.DataFrame(
            {"Sample": protein_values.index, "Intensity": protein_values.values}
        )

        # Merge with metadata
        protein_df = protein_df.merge(metadata_df, on="Sample", how="inner")
        protein_df = protein_df.dropna(subset=["Intensity"])

        if len(protein_df) < 4:
            results.append(_create_empty_result(protein_idx, "Insufficient data"))
            continue

        # Split into groups
        group1_data = protein_df[
            protein_df[config.group_column] == config.group_labels[0]
        ]["Intensity"]
        group2_data = protein_df[
            protein_df[config.group_column] == config.group_labels[1]
        ]["Intensity"]

        if len(group1_data) < 2 or len(group2_data) < 2:
            results.append(_create_empty_result(protein_idx, "Insufficient group data"))
            continue

        try:
            # Mann-Whitney U test
            statistic, p_value = mannwhitneyu(
                group2_data, group1_data, alternative='two-sided'
            )

            # Calculate effect size (r = z / sqrt(N))
            # For Mann-Whitney, we use median and IQR-based effect size
            median1 = group1_data.median()
            median2 = group2_data.median()
            
            # Pooled MAD for effect size
            mad1 = np.median(np.abs(group1_data - median1)) * 1.4826
            mad2 = np.median(np.abs(group2_data - median2)) * 1.4826
            pooled_mad = np.sqrt((mad1**2 + mad2**2) / 2)
            
            effect_size = (median2 - median1) / pooled_mad if pooled_mad > 0 else 0

            # Log fold change using medians
            log_fc = median2 - median1

            result = {
                "Protein": protein_idx,
                "logFC": log_fc,
                "AveExpr": protein_df["Intensity"].mean(),
                "statistic": statistic,
                "P.Value": p_value,
                "B": np.nan,
                "n_group1": len(group1_data),
                "n_group2": len(group2_data),
                "Effect_Size": effect_size,
                "test_method": "Mann-Whitney U",
            }

            results.append(result)

        except Exception as e:
            results.append(_create_empty_result(protein_idx, f"Analysis failed: {e}"))

    print(f"✓ Mann-Whitney U test completed for {len(results)} proteins")
    return pd.DataFrame(results)


def _create_empty_result(protein_idx, reason):
    """Create empty result for failed analysis"""
    return {
        "Protein": protein_idx,
        "logFC": np.nan,
        "AveExpr": np.nan,
        "t": np.nan,
        "P.Value": np.nan,
        "B": np.nan,
        "test_method": f"Failed: {reason}",
    }


def _create_empty_mixed_effects_result(protein_idx, reason):
    """Create empty result for failed mixed-effects analysis"""
    return {
        "Protein": protein_idx,
        "logFC": np.nan,
        "AveExpr": np.nan,
        "t": np.nan,
        "P.Value": np.nan,
        "B": np.nan,
        "group_effect": np.nan,
        "group_pvalue": np.nan,
        "time_effect": np.nan,
        "time_pvalue": np.nan,
        "aic": np.nan,
        "bic": np.nan,
        "n_obs": 0,
        "test_method": f"Mixed-effects failed: {reason}",
    }


def apply_multiple_testing_correction(results_df, config):
    """Apply multiple testing correction"""

    if "P.Value" not in results_df.columns:
        print("Warning: No P.Value column found for correction")
        return results_df

    # Get valid p-values
    valid_pvalues = results_df["P.Value"].dropna()

    if len(valid_pvalues) == 0:
        print("Warning: No valid p-values found")
        results_df["adj.P.Val"] = np.nan
        results_df["Significant"] = False
        return results_df

    # Check if correction should be applied
    correction_method = getattr(config, 'correction_method', config.use_adjusted_pvalue)
    if correction_method == "none" or config.use_adjusted_pvalue == "none":
        # No correction - adjusted p-values are same as raw p-values
        results_df["adj.P.Val"] = results_df["P.Value"]
        results_df["Significant"] = results_df["P.Value"] < config.p_value_threshold
        print("Multiple testing correction applied:")
        print("  Method: none (no correction)")
        print(f"  Significant proteins (p < {config.p_value_threshold}): {(results_df['P.Value'] < config.p_value_threshold).sum()}")
    else:
        # Apply correction
        all_pvalues = results_df["P.Value"].fillna(1.0)
        rejected, adj_pvalues, _, _ = multipletests(
            all_pvalues, method=correction_method
        )

        results_df["adj.P.Val"] = adj_pvalues
        results_df["Significant"] = rejected

        print("Multiple testing correction applied:")
        print(f"  Method: {correction_method}")
        print(
            f"  Significant proteins (FDR < 0.05): {(results_df['adj.P.Val'] < 0.05).sum()}"
        )

    # Add significance categories
    results_df["Significance"] = "Not significant"
    results_df.loc[results_df["adj.P.Val"] < 0.05, "Significance"] = (
        "Significant (FDR < 0.05)"
    )
    results_df.loc[results_df["adj.P.Val"] < 0.01, "Significance"] = (
        "Highly significant (FDR < 0.01)"
    )

    return results_df


def export_results(
    differential_df: pd.DataFrame, output_file: str, include_all: bool = True
) -> None:
    """
    Export differential analysis results to CSV file.

    Parameters:
    -----------
    differential_df : pd.DataFrame
        Differential analysis results
    output_file : str
        Output CSV filename
    include_all : bool
        Whether to include all proteins or only significant ones
    """

    if not include_all:
        export_df = differential_df[differential_df["Significant"]].copy()
        print(f"Exporting {len(export_df)} significant proteins to {output_file}")
    else:
        export_df = differential_df.copy()
        print(f"Exporting all {len(export_df)} proteins to {output_file}")

    export_df.to_csv(output_file, index=False)
    print("Results exported successfully!")


def run_comprehensive_statistical_analysis(
    normalized_data, sample_metadata, config, protein_annotations=None
):
    """
    Comprehensive statistical analysis with automatic dataset validation and subject pairing

    Parameters:
    -----------
    normalized_data : pd.DataFrame
        Protein expression data (proteins x samples)
    sample_metadata : dict
        Dictionary mapping sample names to metadata
    config : StatisticalConfig
        Configuration object with analysis parameters
    protein_annotations : pd.DataFrame, optional
        DataFrame with protein annotations including 'Protein' column

    Returns:
    --------
    pd.DataFrame
        Results of statistical analysis
    """

    print("=" * 60)
    print("COMPREHENSIVE STATISTICAL ANALYSIS")
    print("=" * 60)

    # Step 0: Handle log transformation if needed
    statistical_data = _apply_log_transformation_if_needed(normalized_data, config)

    # Step 1: Clean and validate metadata
    print("Step 1: Cleaning and validating sample metadata...")

    # Clean subject IDs to fix whitespace issues
    cleaned_sample_metadata = {}
    for sample_name, metadata in sample_metadata.items():
        cleaned_metadata = metadata.copy()
        if (
            config.subject_column in cleaned_metadata
            and cleaned_metadata[config.subject_column]
        ):
            cleaned_metadata[config.subject_column] = str(
                cleaned_metadata[config.subject_column]
            ).strip()
        cleaned_sample_metadata[sample_name] = cleaned_metadata

    sample_metadata = cleaned_sample_metadata

    # Step 2: Prepare metadata dataframe
    # IMPORTANT: With standardized data structure, sample columns start at index 5
    # First 5 columns are always: Protein, Description, Protein Gene, UniProt_Accession, UniProt_Entry_Name
    if len(normalized_data.columns) > 5:
        sample_columns = list(
            normalized_data.columns[5:]
        )  # Everything after first 5 annotation columns
        print(
            f"  Using standardized data structure: {len(sample_columns)} sample columns (columns 6+)"
        )
    else:
        # Fallback for legacy data (shouldn't happen with create_standard_data_structure)
        sample_columns = normalized_data.select_dtypes(
            include=[np.number]
        ).columns.tolist()
        print(f"  Using legacy detection: {len(sample_columns)} sample columns")

    print(
        f"  Sample columns: {sample_columns[:3]}{'...' if len(sample_columns) > 3 else ''}"
    )

    metadata_df = prepare_metadata_dataframe(sample_metadata, sample_columns, config)

    # Step 3: Analyze experimental design
    print("\nStep 2: Analyzing experimental design...")

    # Filter to experimental samples only (exclude controls)
    valid_samples = {}
    for sample_name, metadata in sample_metadata.items():
        comparison_value = metadata.get(config.group_column)
        if comparison_value in config.group_labels:
            valid_samples[sample_name] = metadata

    print(f"  Valid experimental samples: {len(valid_samples)}")

    # Analyze subject pairing structure
    pairing_data = {}
    for sample_name, metadata in valid_samples.items():
        subject = metadata.get(config.subject_column)
        visit = metadata.get(config.paired_column)
        comparison = metadata.get(config.group_column)

        if subject and visit and comparison:
            if subject not in pairing_data:
                pairing_data[subject] = {}
            pairing_data[subject][visit] = {
                "sample": sample_name,
                "comparison": comparison,
            }

    # Check for complete pairs
    complete_pairs = []
    incomplete_subjects = []

    for subject, visits in pairing_data.items():
        if config.paired_label1 in visits and config.paired_label2 in visits:
            baseline = visits[config.paired_label1]
            followup = visits[config.paired_label2]

            if baseline["comparison"] == followup["comparison"]:
                complete_pairs.append(
                    {
                        "subject": subject,
                        "group": baseline["comparison"],
                        "baseline_sample": baseline["sample"],
                        "followup_sample": followup["sample"],
                    }
                )
            else:
                incomplete_subjects.append(f"{subject} (mixed groups)")
        else:
            available_visits = list(visits.keys())
            incomplete_subjects.append(
                f"{subject} (missing visits: {available_visits})"
            )

    # Group complete pairs by treatment group
    group_pairs = {
        group: [p for p in complete_pairs if p["group"] == group]
        for group in config.group_labels
    }

    print("  Complete paired subjects by group:")
    for group, pairs in group_pairs.items():
        print(f"    {group}: {len(pairs)} subjects")

    if incomplete_subjects:
        print(f"  Incomplete subjects: {len(incomplete_subjects)}")
        if len(incomplete_subjects) <= 5:  # Show details if few
            for subject_info in incomplete_subjects:
                print(f"    {subject_info}")

    # Step 4: Run statistical analysis
    print(f"\nStep 3: Running {config.statistical_test_method} analysis...")

    # Filter protein data to samples with metadata
    available_samples = metadata_df["Sample"].tolist()
    filtered_protein_data = statistical_data[available_samples]

    print(f"  Method: {config.statistical_test_method}")
    print(f"  Proteins: {len(filtered_protein_data)}")
    print(f"  Samples: {len(available_samples)}")
    print(f"  Groups: {config.group_labels}")
    print(f"  Timepoints: {config.paired_label1} -> {config.paired_label2}")

    if config.statistical_test_method == "mixed_effects":
        print(
            f"  Interaction terms: {config.interaction_terms + config.additional_interactions}"
        )
        if config.covariates:
            print(f"  Covariates: {config.covariates}")

    # Run appropriate analysis
    if config.statistical_test_method == "mixed_effects":
        results_df = run_mixed_effects_analysis(
            filtered_protein_data, metadata_df, config, protein_annotations
        )
    elif config.statistical_test_method in ["paired_t", "paired_welch"]:
        results_df = run_paired_t_test(filtered_protein_data, metadata_df, config)
    elif config.statistical_test_method in ["welch_t", "student_t"]:
        results_df = run_unpaired_t_test(filtered_protein_data, metadata_df, config)
    elif config.statistical_test_method == "wilcoxon":
        results_df = run_wilcoxon_test(filtered_protein_data, metadata_df, config)
    elif config.statistical_test_method == "mann_whitney":
        results_df = run_mann_whitney_test(filtered_protein_data, metadata_df, config)
    else:
        raise ValueError(
            f"Unknown statistical method: {config.statistical_test_method}. "
            f"Supported methods: mixed_effects, paired_t, paired_welch, welch_t, student_t, wilcoxon, mann_whitney"
        )

    # Apply multiple testing correction
    results_df = apply_multiple_testing_correction(results_df, config)

    # Sort by p-value
    results_df = results_df.sort_values("P.Value")

    print("\n✓ Statistical analysis completed!")
    print(f"  Total proteins analyzed: {len(results_df)}")
    print(f"  Proteins with valid results: {results_df['P.Value'].notna().sum()}")
    print(
        f"  Significant proteins (FDR < 0.05): {(results_df['adj.P.Val'] < 0.05).sum()}"
    )

    return results_df


def display_analysis_summary(differential_results, config, label_top_n=10):
    """
    Display comprehensive summary of statistical analysis results

    Parameters:
    -----------
    differential_results : pd.DataFrame
        Results from statistical analysis
    config : StatisticalConfig
        Configuration object with analysis parameters
    label_top_n : int
        Number of top significant proteins to display

    Returns:
    --------
    dict
        Summary statistics for downstream use
    """

    if differential_results is None or len(differential_results) == 0:
        print("⚠️ No differential analysis results available")
        return {}

    print("=" * 60)
    print("STATISTICAL ANALYSIS SUMMARY")
    print("=" * 60)

    # Basic statistics
    total_proteins = len(differential_results)
    valid_results = differential_results["P.Value"].notna().sum()
    significant_005 = (differential_results["adj.P.Val"] < 0.05).sum()
    significant_001 = (differential_results["adj.P.Val"] < 0.01).sum()

    print("Analysis Overview:")
    print(f"  Method: {config.statistical_test_method.upper()}")
    print(f"  Total proteins analyzed: {total_proteins:,}")
    print(f"  Proteins with valid results: {valid_results:,}")
    print(f"  Significant proteins (FDR < 0.05): {significant_005:,}")
    print(f"  Highly significant (FDR < 0.01): {significant_001:,}")

    if valid_results == 0:
        print("\n❌ No valid statistical results found")
        return {}

    # Show top significant results
    successful_results = differential_results[differential_results["P.Value"].notna()]

    if len(successful_results) > 0:
        print(f"\n=== TOP {label_top_n} MOST SIGNIFICANT PROTEINS ===")

        top_results = successful_results.nsmallest(label_top_n, "P.Value")

        # Choose appropriate columns based on analysis type
        if config.statistical_test_method == "mixed_effects":
            # Mixed-effects model results - use primary columns only
            display_cols = ["Protein", "logFC", "P.Value", "adj.P.Val", "n_obs"]

            # Note: logFC and P.Value already contain interaction results
            # No need to show duplicate interaction_coef and interaction_pvalue

        else:
            # Traditional statistical test results
            display_cols = ["Protein", "logFC", "P.Value", "adj.P.Val"]
            if "Effect_Size" in top_results.columns:
                display_cols.insert(2, "Effect_Size")

        # Filter to available columns
        available_cols = [col for col in display_cols if col in top_results.columns]

        if available_cols:
            # Create a clean display dataframe with only the specified columns
            display_df = pd.DataFrame()
            for col in available_cols:
                display_df[col] = top_results[col].copy()

            # Format numerical columns for better display
            for col in display_df.columns:
                if col in ["P.Value", "adj.P.Val"]:
                    display_df[col] = display_df[col].apply(
                        lambda x: f"{x:.2e}"
                        if pd.notna(x) and x < 0.01
                        else f"{x:.6f}"
                        if pd.notna(x)
                        else "N/A"
                    )
                elif col in ["logFC", "Effect_Size"]:
                    display_df[col] = display_df[col].apply(
                        lambda x: f"{x:.4f}" if pd.notna(x) else "N/A"
                    )

            print(display_df.to_string(index=False))

        # Additional analysis-specific summary
        if config.statistical_test_method == "mixed_effects":
            # Use P.Value column since it contains the interaction p-values
            interaction_significant = (successful_results["P.Value"] < 0.05).sum()
            print("\nInteraction Effects:")
            print(
                f"  Significant {' × '.join(config.interaction_terms)} interactions: {interaction_significant}"
            )

    else:
        print("\n❌ No proteins with valid statistical results")

        # Show failure reasons if available
        failed_results = differential_results[differential_results["P.Value"].isna()]
        if "test_method" in failed_results.columns and len(failed_results) > 0:
            print("\nFailure Analysis:")
            failure_reasons = failed_results["test_method"].value_counts()
            for reason, count in failure_reasons.items():
                print(f"  {reason}: {count}")

    # Create summary dictionary for return
    summary = {
        "total_proteins": total_proteins,
        "valid_results": valid_results,
        "significant_005": significant_005,
        "significant_001": significant_001,
        "analysis_method": config.statistical_test_method,
        "success_rate": valid_results / total_proteins if total_proteins > 0 else 0,
    }

    print("\n✓ Analysis summary complete!")

    return summary


# Maintain backwards compatibility
def run_statistical_analysis(
    normalized_data, sample_metadata, config, protein_annotations=None
):
    """Backwards compatible wrapper for run_comprehensive_statistical_analysis"""
    return run_comprehensive_statistical_analysis(
        normalized_data, sample_metadata, config, protein_annotations
    )
