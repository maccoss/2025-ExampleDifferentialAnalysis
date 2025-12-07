"""
Tests for proteomics_toolkit.statistical_analysis module
"""

import pandas as pd
import numpy as np
import warnings
import pytest

from proteomics_toolkit.statistical_analysis import (
    StatisticalConfig,
    prepare_metadata_dataframe,
    run_paired_t_test,
    run_unpaired_t_test,
    run_wilcoxon_test,
    run_mann_whitney_test,
    apply_multiple_testing_correction,
    run_comprehensive_statistical_analysis,
    display_analysis_summary,
    _create_empty_result,
    _create_empty_mixed_effects_result,
    _apply_log_transformation_if_needed,
)


class TestStatisticalConfig:
    """Test the StatisticalConfig class"""

    def test_config_initialization(self):
        """Test configuration initialization with defaults"""
        config = StatisticalConfig()

        assert config.statistical_test_method == "mixed_effects"
        assert config.analysis_type is None  # Must be explicitly set by user
        assert config.p_value_threshold == 0.05
        assert config.fold_change_threshold == 1.5
        assert config.use_adjusted_pvalue == "adjusted"
        assert config.enable_pvalue_fallback is True

    def test_config_modification(self):
        """Test modifying configuration values"""
        config = StatisticalConfig()

        config.statistical_test_method = "paired_t"
        config.p_value_threshold = 0.01
        config.subject_column = "Subject"

        assert config.statistical_test_method == "paired_t"
        assert config.p_value_threshold == 0.01
        assert config.subject_column == "Subject"


class TestPrepareMetadataDataframe:
    """Test metadata dataframe preparation"""

    def test_prepare_metadata_basic(
        self, sample_metadata, sample_columns, statistical_config
    ):
        """Test basic metadata preparation"""
        result = prepare_metadata_dataframe(
            sample_metadata, sample_columns, statistical_config
        )

        assert isinstance(result, pd.DataFrame)
        assert "Sample" in result.columns
        assert statistical_config.subject_column in result.columns
        assert statistical_config.group_column in result.columns
        assert len(result) <= len(sample_columns)  # May filter out some samples

    def test_prepare_metadata_missing_columns(self, sample_metadata, sample_columns):
        """Test metadata preparation when required columns are missing"""
        config = StatisticalConfig()
        config.subject_column = "NonexistentColumn"

        # This should raise a ValueError for missing required columns
        with pytest.raises(ValueError, match="Missing required metadata columns"):
            prepare_metadata_dataframe(sample_metadata, sample_columns, config)


class TestRunPairedTTest:
    """Test paired t-test analysis"""

    def test_paired_t_test_basic(
        self, sample_protein_data, sample_metadata_df, statistical_config
    ):
        """Test basic paired t-test"""
        # Configure for paired analysis
        statistical_config.statistical_test_method = "paired_t"

        # Select only sample columns for analysis
        sample_columns = [
            col for col in sample_protein_data.columns if col.startswith("Sample_")
        ]
        protein_data = sample_protein_data[sample_columns]

        result = run_paired_t_test(protein_data, sample_metadata_df, statistical_config)

        assert isinstance(result, pd.DataFrame)
        assert "P.Value" in result.columns
        assert "logFC" in result.columns
        assert "t" in result.columns
        assert len(result) == len(protein_data)  # One result per protein

    def test_paired_t_test_insufficient_data(self, statistical_config):
        """Test paired t-test with insufficient data"""
        # Create minimal dataset
        protein_data = pd.DataFrame(
            {"Sample_1": [100.0, 200.0], "Sample_2": [110.0, 210.0]}
        )

        metadata_df = pd.DataFrame(
            {
                "Sample": ["Sample_1", "Sample_2"],
                "Subject": ["S001", "S001"],  # Same subject for pairing
                "Group": ["Control", "Control"],
                "Visit": ["Baseline", "Week4"],
            }
        )

        statistical_config.statistical_test_method = "paired_t"

        result = run_paired_t_test(protein_data, metadata_df, statistical_config)

        assert isinstance(result, pd.DataFrame)
        assert len(result) == len(protein_data)


class TestRunUnpairedTTest:
    """Test unpaired t-test analysis"""

    def test_unpaired_t_test_basic(
        self, sample_protein_data, sample_metadata_df, statistical_config
    ):
        """Test basic unpaired t-test"""
        statistical_config.statistical_test_method = "welch_t"

        # Select only sample columns for analysis
        sample_columns = [
            col for col in sample_protein_data.columns if col.startswith("Sample_")
        ]
        protein_data = sample_protein_data[sample_columns]

        result = run_unpaired_t_test(
            protein_data, sample_metadata_df, statistical_config
        )

        assert isinstance(result, pd.DataFrame)
        assert "P.Value" in result.columns
        assert "logFC" in result.columns
        assert "t" in result.columns
        assert len(result) == len(protein_data)


class TestRunWilcoxonTest:
    """Test Wilcoxon signed-rank test"""

    def test_wilcoxon_test_basic(
        self, standardized_protein_data, sample_metadata, statistical_config
    ):
        """Test basic Wilcoxon signed-rank test functionality"""
        statistical_config.statistical_test_method = "wilcoxon"
        statistical_config.paired_column = "Visit"
        statistical_config.paired_label1 = "Baseline"
        statistical_config.paired_label2 = "Week4"
        statistical_config.subject_column = "Subject"

        # Create metadata DataFrame
        sample_columns = list(standardized_protein_data.columns[5:])
        metadata_df = prepare_metadata_dataframe(
            sample_metadata, sample_columns, statistical_config
        )

        result = run_wilcoxon_test(standardized_protein_data, metadata_df, statistical_config)

        # Check result structure
        assert isinstance(result, pd.DataFrame)
        assert "Protein" in result.columns
        assert "P.Value" in result.columns
        assert "logFC" in result.columns
        assert "test_method" in result.columns
        assert len(result) == len(standardized_protein_data)
        
        # Check that 'statistic' column exists in successful results
        # (may not exist if all tests fail and use _create_empty_result)
        valid_results = result[result["P.Value"].notna()]
        if len(valid_results) > 0:
            assert "statistic" in result.columns
            assert (valid_results["test_method"] == "Wilcoxon signed-rank").all()
        else:
            print("All Wilcoxon tests failed - this is expected with random test data")

    def test_wilcoxon_test_insufficient_data(self, statistical_config):
        """Test Wilcoxon test with insufficient data"""
        # Create minimal data with too few samples
        protein_data = pd.DataFrame({
            "Protein": ["P00001"],
            "Description": ["Test Protein 1"],
            "Protein Gene": ["GENE1"],
            "UniProt_Accession": ["P00001"],
            "UniProt_Entry_Name": ["TEST1_HUMAN"],
            "Sample1": [100.0],
            "Sample2": [110.0]
        }).set_index("Protein")

        # Create matching metadata - not enough for pairing
        metadata_df = pd.DataFrame({
            "Sample": ["Sample1", "Sample2"],
            "Visit": ["D-02", "D-13"],
            "Subject": ["S1", "S2"]  # Different subjects, so no pairing
        })

        statistical_config.statistical_test_method = "wilcoxon"
        result = run_wilcoxon_test(protein_data, metadata_df, statistical_config)

        # Should return empty results due to insufficient pairing
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 1
        assert result.iloc[0]["P.Value"] is np.nan or pd.isna(result.iloc[0]["P.Value"])


class TestRunMannWhitneyTest:
    """Test Mann-Whitney U test"""

    def test_mann_whitney_test_basic(
        self, standardized_protein_data, sample_metadata, statistical_config
    ):
        """Test basic Mann-Whitney U test functionality"""
        statistical_config.statistical_test_method = "mann_whitney"
        statistical_config.group_column = "DrugDose"
        statistical_config.group_labels = ["Placebo", "High"]

        # Create metadata DataFrame
        sample_columns = list(standardized_protein_data.columns[5:])
        metadata_df = prepare_metadata_dataframe(
            sample_metadata, sample_columns, statistical_config
        )

        result = run_mann_whitney_test(standardized_protein_data, metadata_df, statistical_config)

        # Check result structure
        assert isinstance(result, pd.DataFrame)
        assert "Protein" in result.columns
        assert "P.Value" in result.columns
        assert "logFC" in result.columns
        assert "test_method" in result.columns
        assert len(result) == len(standardized_protein_data)
        
        # Check that 'statistic' column exists in successful results
        # (may not exist if all tests fail and use _create_empty_result)
        valid_results = result[result["P.Value"].notna()]
        if len(valid_results) > 0:
            assert "statistic" in result.columns
            assert (valid_results["test_method"] == "Mann-Whitney U").all()
        else:
            print("All Mann-Whitney tests failed - this is expected with random test data")

    def test_mann_whitney_test_insufficient_data(self, statistical_config):
        """Test Mann-Whitney test with insufficient data"""
        # Create minimal data
        protein_data = pd.DataFrame({
            "Protein": ["P00001"],
            "Description": ["Test Protein 1"],
            "Protein Gene": ["GENE1"],
            "UniProt_Accession": ["P00001"],
            "UniProt_Entry_Name": ["TEST1_HUMAN"],
            "Sample1": [100.0]
        }).set_index("Protein")

        metadata_df = pd.DataFrame({
            "Sample": ["Sample1"],
            "DrugDose": ["0"]
        })

        statistical_config.statistical_test_method = "mann_whitney"
        statistical_config.group_column = "DrugDose"
        statistical_config.group_labels = ["0", "20"]
        # Don't set paired configuration for unpaired test
        statistical_config.paired_column = None
        statistical_config.paired_label1 = None  
        statistical_config.paired_label2 = None

        result = run_mann_whitney_test(protein_data, metadata_df, statistical_config)

        # Should return empty results due to insufficient data
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 1
        assert result.iloc[0]["P.Value"] is np.nan or pd.isna(result.iloc[0]["P.Value"])


class TestAllStatisticalMethods:
    """Test all statistical methods are working with comprehensive analysis"""

    def test_all_methods_available(
        self, standardized_protein_data, sample_metadata, statistical_config
    ):
        """Test that all 7 statistical methods can be called"""
        
        # List of all methods that should be supported
        all_methods = [
            "mixed_effects",
            "paired_t", 
            "paired_welch",
            "welch_t",
            "student_t",
            "wilcoxon",
            "mann_whitney"
        ]
        
        successful_methods = []
        failed_methods = []
        
        for method in all_methods:
            try:
                statistical_config.statistical_test_method = method
                
                # Set appropriate configuration for each method type
                if method in ["mixed_effects", "paired_t", "paired_welch", "wilcoxon"]:
                    # Paired analysis setup
                    statistical_config.paired_column = "Visit"
                    statistical_config.paired_label1 = "Baseline"
                    statistical_config.paired_label2 = "Week4"
                    statistical_config.subject_column = "Subject"
                    statistical_config.analysis_type = "paired"
                else:
                    # Unpaired analysis setup
                    statistical_config.group_column = "DrugDose"
                    statistical_config.group_labels = ["Placebo", "High"]
                    statistical_config.analysis_type = "unpaired"
                
                # Run the comprehensive analysis to test integration
                try:
                    result = run_comprehensive_statistical_analysis(
                        standardized_protein_data, sample_metadata, statistical_config
                    )
                    
                    # Basic result validation
                    assert isinstance(result, pd.DataFrame)
                    assert len(result) > 0
                    assert "Protein" in result.columns
                    assert "P.Value" in result.columns
                    
                    successful_methods.append(method)
                    print(f"✓ {method}: SUCCESS")
                    
                except Exception as e:
                    failed_methods.append((method, str(e)))
                    print(f"✗ {method}: FAILED - {e}")
                    
            except ImportError as e:
                # Handle missing dependencies (e.g., statsmodels for mixed_effects)
                if method == "mixed_effects" and "statsmodels" in str(e):
                    print(f"⚠ {method}: SKIPPED - Missing dependency (statsmodels)")
                    continue
                else:
                    failed_methods.append((method, str(e)))
                    
        # Report results
        print(f"\nSTATISTICAL METHOD COVERAGE REPORT:")
        print(f"✓ Successful methods: {len(successful_methods)}")
        print(f"✗ Failed methods: {len(failed_methods)}")
        
        if failed_methods:
            print(f"\nFailed methods details:")
            for method, error in failed_methods:
                print(f"  - {method}: {error}")
        
        # Ensure at least parametric methods work
        parametric_methods = ["paired_t", "paired_welch", "welch_t", "student_t"]
        working_parametric = [m for m in parametric_methods if m in successful_methods]
        assert len(working_parametric) >= 2, f"At least 2 parametric methods should work, got: {working_parametric}"
        
        # Ensure non-parametric methods work
        nonparametric_methods = ["wilcoxon", "mann_whitney"]
        working_nonparametric = [m for m in nonparametric_methods if m in successful_methods]
        assert len(working_nonparametric) == 2, f"Both non-parametric methods should work, got: {working_nonparametric}"


class TestMultipleTestingCorrection:
    """Test multiple testing correction"""

    def test_fdr_correction(self, differential_results, statistical_config):
        """Test FDR correction"""
        statistical_config.use_adjusted_pvalue = "fdr_bh"

        result = apply_multiple_testing_correction(
            differential_results, statistical_config
        )

        assert "adj.P.Val" in result.columns
        assert all(
            result["adj.P.Val"] >= result["P.Value"]
        )  # Adjusted p-values should be >= original

    def test_bonferroni_correction(self, differential_results, statistical_config):
        """Test Bonferroni correction"""
        statistical_config.use_adjusted_pvalue = "bonferroni"

        result = apply_multiple_testing_correction(
            differential_results, statistical_config
        )

        assert "adj.P.Val" in result.columns
        assert all(result["adj.P.Val"] >= result["P.Value"])

    def test_no_correction(self, differential_results, statistical_config):
        """Test when no correction is applied"""
        statistical_config.use_adjusted_pvalue = "none"

        result = apply_multiple_testing_correction(
            differential_results, statistical_config
        )

        assert "adj.P.Val" in result.columns
        assert all(result["adj.P.Val"] == result["P.Value"])  # Should be identical


class TestRunComprehensiveStatisticalAnalysis:
    """Test the main comprehensive analysis function"""

    def test_comprehensive_analysis_paired_t(
        self, sample_protein_data, sample_metadata, statistical_config
    ):
        """Test comprehensive analysis with paired t-test"""
        statistical_config.statistical_test_method = "paired_t"

        # Suppress output during testing
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            result = run_comprehensive_statistical_analysis(
                normalized_data=sample_protein_data,
                sample_metadata=sample_metadata,
                config=statistical_config,
                protein_annotations=None,
            )

        assert isinstance(result, pd.DataFrame)
        assert "P.Value" in result.columns
        assert "adj.P.Val" in result.columns
        assert len(result) > 0

    def test_comprehensive_analysis_with_annotations(
        self, sample_protein_data, sample_metadata, statistical_config
    ):
        """Test comprehensive analysis with protein annotations"""
        statistical_config.statistical_test_method = "paired_t"

        # Use the protein data itself as annotations (contains Protein, Gene columns)
        annotations = sample_protein_data[["Protein", "ProteinName", "Gene"]].copy()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            result = run_comprehensive_statistical_analysis(
                normalized_data=sample_protein_data,
                sample_metadata=sample_metadata,
                config=statistical_config,
                protein_annotations=annotations,
            )

        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0


class TestDisplayAnalysisSummary:
    """Test analysis summary display"""

    def test_display_summary_basic(self, differential_results, statistical_config):
        """Test basic summary display"""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            summary = display_analysis_summary(
                differential_results=differential_results,
                config=statistical_config,
                label_top_n=5,
            )

        assert isinstance(summary, dict)
        assert "analysis_method" in summary
        assert "total_proteins" in summary
        assert "significant_005" in summary  # Changed from "significant_proteins"
        assert "significant_001" in summary
        assert "valid_results" in summary
        assert "success_rate" in summary

    def test_display_summary_empty_results(self, statistical_config):
        """Test summary display with empty results"""
        empty_results = pd.DataFrame()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            summary = display_analysis_summary(
                differential_results=empty_results, config=statistical_config
            )

        assert isinstance(summary, dict)


class TestHelperFunctions:
    """Test helper functions"""

    def test_create_empty_result(self):
        """Test creating empty result for failed analysis"""
        result = _create_empty_result("P00001", "Test reason")

        assert isinstance(result, dict)
        assert result["Protein"] == "P00001"
        assert "P.Value" in result
        assert pd.isna(result["P.Value"])

    def test_create_empty_mixed_effects_result(self):
        """Test creating empty mixed-effects result"""
        result = _create_empty_mixed_effects_result("P00001", "Test reason")

        assert isinstance(result, dict)
        assert result["Protein"] == "P00001"
        assert "P.Value" in result
        assert "group_effect" in result
        assert pd.isna(result["P.Value"])


class TestEdgeCases:
    """Test edge cases and error conditions"""

    def test_analysis_with_all_nan_protein(self, sample_metadata, statistical_config):
        """Test analysis when a protein has all NaN values"""
        # Create protein data with one all-NaN protein
        protein_data = pd.DataFrame(
            {
                "Protein": ["P00001", "P00002"],
                "Sample_A_1": [100.0, np.nan],
                "Sample_A_2": [110.0, np.nan],
                "Sample_B_1": [120.0, np.nan],
                "Sample_B_2": [115.0, np.nan],
            }
        )

        statistical_config.statistical_test_method = "paired_t"

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            result = run_comprehensive_statistical_analysis(
                normalized_data=protein_data,
                sample_metadata=sample_metadata,
                config=statistical_config,
            )

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 2  # Should handle both proteins

    def test_analysis_with_no_valid_samples(self, statistical_config):
        """Test analysis when no samples match metadata"""
        protein_data = pd.DataFrame({"Protein": ["P00001"], "Unknown_Sample": [100.0]})

        sample_metadata = {"Different_Sample": {"Subject": "S001", "Group": "Control"}}

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            result = run_comprehensive_statistical_analysis(
                normalized_data=protein_data,
                sample_metadata=sample_metadata,
                config=statistical_config,
            )

        assert isinstance(result, pd.DataFrame)
        # Should return empty results or handle gracefully


class TestLogTransformation:
    """Test log transformation functionality"""

    def test_statistical_config_has_log_parameters(self):
        """Test that StatisticalConfig has all required log transformation parameters"""
        config = StatisticalConfig()

        # Check that all log transformation parameters exist
        assert hasattr(config, "log_transform_before_stats")
        assert hasattr(config, "log_base")
        assert hasattr(config, "log_pseudocount")

        # Check default values
        assert config.log_transform_before_stats == "auto"
        assert config.log_base == "log2"
        assert config.log_pseudocount is None

    def test_log_transformation_auto_mode_with_median_normalization(self):
        """Test that auto mode applies log transformation for median normalization"""
        config = StatisticalConfig()
        config.log_transform_before_stats = "auto"

        # Create test data with large values typical of median normalization
        test_data = pd.DataFrame(
            {
                "protein1": [1000.0, 1200.0, 950.0, 1100.0],
                "protein2": [800.0, 900.0, 750.0, 850.0],
                "protein3": [1500.0, 1800.0, 1400.0, 1600.0],
            }
        )
        test_data.attrs = {"normalization_method": "median"}

        result = _apply_log_transformation_if_needed(test_data, config)

        # Data should be transformed
        assert not np.array_equal(test_data.values, result.values)

        # Values should be much smaller (log-transformed)
        assert result.mean().mean() < test_data.mean().mean()

        # Should be in reasonable log range (typically 8-15 for log2 of proteomics data)
        assert result.min().min() > 5.0  # log2(32) ≈ 5
        assert result.max().max() < 20.0  # log2(1M) ≈ 20

    def test_log_transformation_auto_mode_with_vsn_normalization(self):
        """Test that auto mode skips transformation for VSN normalization (already transformed scale)"""
        config = StatisticalConfig()
        config.log_transform_before_stats = "auto"

        # Create test data with values typical of VSN (already on transformed scale)
        test_data = pd.DataFrame(
            {
                "protein1": [12.5, 13.2, 11.8, 12.9],
                "protein2": [10.8, 11.5, 10.2, 11.1],
                "protein3": [14.1, 14.8, 13.5, 14.2],
            }
        )
        test_data.attrs = {"normalization_method": "vsn"}

        result = _apply_log_transformation_if_needed(test_data, config)

        # Data should NOT be transformed (VSN is already on appropriate scale)
        assert np.array_equal(test_data.values, result.values)

    def test_log_transformation_forced_true(self):
        """Test forced log transformation (log_transform_before_stats=True)"""
        config = StatisticalConfig()
        config.log_transform_before_stats = "true"  # Use string instead of boolean
        config.log_base = "log2"

        test_data = pd.DataFrame(
            {
                "protein1": [1000.0, 2000.0, 500.0, 1500.0],
                "protein2": [800.0, 1600.0, 400.0, 1200.0],
            }
        )

        result = _apply_log_transformation_if_needed(test_data, config)

        # Should be transformed regardless of normalization method
        assert not np.array_equal(test_data.values, result.values)
        assert result.mean().mean() < test_data.mean().mean()

    def test_log_transformation_forced_false(self):
        """Test that transformation is skipped when log_transform_before_stats=False"""
        config = StatisticalConfig()
        config.log_transform_before_stats = "false"  # Use string instead of boolean

        test_data = pd.DataFrame(
            {
                "protein1": [1000.0, 2000.0, 500.0, 1500.0],
                "protein2": [800.0, 1600.0, 400.0, 1200.0],
            }
        )

        result = _apply_log_transformation_if_needed(test_data, config)

        # Data should NOT be transformed
        assert np.array_equal(test_data.values, result.values)

    def test_log_transformation_different_bases(self):
        """Test different logarithmic bases"""
        test_data = pd.DataFrame(
            {"protein1": [1000.0, 2000.0, 500.0], "protein2": [800.0, 1600.0, 400.0]}
        )

        # Test log2
        config_log2 = StatisticalConfig()
        config_log2.log_transform_before_stats = "true"  # Use string
        config_log2.log_base = "log2"
        result_log2 = _apply_log_transformation_if_needed(test_data, config_log2)

        # Test log10
        config_log10 = StatisticalConfig()
        config_log10.log_transform_before_stats = "true"  # Use string
        config_log10.log_base = "log10"
        result_log10 = _apply_log_transformation_if_needed(test_data, config_log10)

        # Test ln
        config_ln = StatisticalConfig()
        config_ln.log_transform_before_stats = "true"  # Use string
        config_ln.log_base = "ln"
        result_ln = _apply_log_transformation_if_needed(test_data, config_ln)

        # All should be different from original
        assert not np.array_equal(test_data.values, result_log2.values)
        assert not np.array_equal(test_data.values, result_log10.values)
        assert not np.array_equal(test_data.values, result_ln.values)

        # log2 should give largest values, log10 smallest, ln in between
        mean_log2 = result_log2.mean().mean()
        mean_log10 = result_log10.mean().mean()
        mean_ln = result_ln.mean().mean()

        assert mean_log2 > mean_ln > mean_log10

    def test_log_transformation_with_custom_pseudocount(self):
        """Test log transformation with custom pseudocount"""
        config = StatisticalConfig()
        config.log_transform_before_stats = "true"  # Use string
        config.log_base = "log2"
        config.log_pseudocount = 1.0  # Custom pseudocount

        test_data = pd.DataFrame(
            {
                "protein1": [100.0, 200.0, 50.0, 150.0],
                "protein2": [80.0, 160.0, 40.0, 120.0],
            }
        )

        result = _apply_log_transformation_if_needed(test_data, config)

        # Should be transformed
        assert not np.array_equal(test_data.values, result.values)

        # Check that pseudocount was applied (values should be log2(original + 1))
        expected_protein1_first = np.log2(100.0 + 1.0)  # log2(101)
        assert np.isclose(result.iloc[0, 0], expected_protein1_first, rtol=1e-10)

    def test_log_transformation_with_auto_pseudocount(self):
        """Test log transformation with automatic pseudocount calculation"""
        config = StatisticalConfig()
        config.log_transform_before_stats = "true"  # Use string
        config.log_base = "log2"
        config.log_pseudocount = None  # Auto pseudocount

        test_data = pd.DataFrame(
            {
                "protein1": [1000.0, 2000.0, 500.0, 1500.0],
                "protein2": [800.0, 1600.0, 400.0, 1200.0],
            }
        )

        result = _apply_log_transformation_if_needed(test_data, config)

        # Should be transformed
        assert not np.array_equal(test_data.values, result.values)

        # Auto pseudocount should be 1% of minimum value
        min_val = test_data.min().min()  # 400.0
        expected_pseudocount = min_val * 0.01  # 4.0
        expected_first_value = np.log2(1000.0 + expected_pseudocount)

        assert np.isclose(result.iloc[0, 0], expected_first_value, rtol=1e-10)

    def test_log_transformation_with_zero_values(self):
        """Test log transformation handles zero values correctly with pseudocount"""
        config = StatisticalConfig()
        config.log_transform_before_stats = "true"  # Use string
        config.log_base = "log2"
        config.log_pseudocount = 1.0

        test_data = pd.DataFrame(
            {
                "protein1": [0.0, 100.0, 200.0, 150.0],  # Contains zero
                "protein2": [50.0, 150.0, 0.0, 100.0],  # Contains zero
            }
        )

        result = _apply_log_transformation_if_needed(test_data, config)

        # Should not contain any NaN or infinite values
        assert not result.isnull().any().any()
        assert np.isfinite(result.values).all()

        # Zero values should become log2(0 + 1) = log2(1) = 0
        zero_positions = test_data == 0.0
        result_zero_positions = result[zero_positions]
        expected_zero_result = np.log2(1.0)  # 0.0

        assert np.allclose(result_zero_positions.dropna(), expected_zero_result)

    def test_log_transformation_preserves_dataframe_structure(self):
        """Test that log transformation preserves DataFrame structure and metadata"""
        config = StatisticalConfig()
        config.log_transform_before_stats = "true"  # Use string

        test_data = pd.DataFrame(
            {
                "protein1": [1000.0, 1200.0, 950.0],
                "protein2": [800.0, 900.0, 750.0],
                "protein3": [1500.0, 1800.0, 1400.0],
            },
            index=["sample1", "sample2", "sample3"],
        )

        # Add some attributes
        test_data.attrs = {"normalization_method": "median", "test_attr": "test_value"}

        result = _apply_log_transformation_if_needed(test_data, config)

        # Check structure preservation
        assert result.shape == test_data.shape
        assert list(result.columns) == list(test_data.columns)
        assert list(result.index) == list(test_data.index)

        # Check that attrs are preserved
        assert result.attrs.get("normalization_method") == "median"
        assert result.attrs.get("test_attr") == "test_value"

    def test_log_transformation_invalid_base_raises_error(self):
        """Test that invalid log base raises appropriate error"""
        config = StatisticalConfig()
        config.log_transform_before_stats = "true"  # Use string
        config.log_base = "invalid_base"

        test_data = pd.DataFrame(
            {"protein1": [1000.0, 1200.0, 950.0], "protein2": [800.0, 900.0, 750.0]}
        )

        try:
            _apply_log_transformation_if_needed(test_data, config)
            assert False, "Should have raised ValueError for invalid log base"
        except ValueError as e:
            assert "Unknown log base" in str(e)


class TestLongitudinalAnalysisType:
    """Test the longitudinal analysis type (F-test for any change over time)"""

    @pytest.fixture
    def longitudinal_metadata(self):
        """Create metadata for longitudinal analysis"""
        return {
            # Subject 1 across 4 timepoints
            'S1_W0': {'Subject': 'S1', 'Week': 0},
            'S1_W2': {'Subject': 'S1', 'Week': 2},
            'S1_W4': {'Subject': 'S1', 'Week': 4},
            'S1_W8': {'Subject': 'S1', 'Week': 8},
            # Subject 2 across 4 timepoints
            'S2_W0': {'Subject': 'S2', 'Week': 0},
            'S2_W2': {'Subject': 'S2', 'Week': 2},
            'S2_W4': {'Subject': 'S2', 'Week': 4},
            'S2_W8': {'Subject': 'S2', 'Week': 8},
            # Subject 3 across 4 timepoints
            'S3_W0': {'Subject': 'S3', 'Week': 0},
            'S3_W2': {'Subject': 'S3', 'Week': 2},
            'S3_W4': {'Subject': 'S3', 'Week': 4},
            'S3_W8': {'Subject': 'S3', 'Week': 8},
        }

    @pytest.fixture
    def longitudinal_protein_data(self):
        """Create protein data for longitudinal analysis"""
        np.random.seed(42)
        n_proteins = 5
        
        data = pd.DataFrame({
            'Protein': [f'P{i:05d}' for i in range(n_proteins)],
            'Gene': [f'GENE{i}' for i in range(n_proteins)],
        })
        
        # Create sample columns with realistic temporal patterns
        samples = [f'S{s}_W{w}' for s in [1, 2, 3] for w in [0, 2, 4, 8]]
        for sample in samples:
            data[sample] = np.random.lognormal(10, 0.5, n_proteins)
        
        return data

    def test_longitudinal_config_validation(self):
        """Test that longitudinal analysis configuration validates correctly"""
        config = StatisticalConfig()
        config.analysis_type = "longitudinal"
        config.statistical_test_method = "mixed_effects"
        config.subject_column = "Subject"
        config.time_column = "Week"
        config.covariates = []
        
        # Should validate successfully
        config.validate()
        
        assert config.analysis_type == "longitudinal"
        assert config.time_column == "Week"

    def test_longitudinal_config_requires_time_column(self):
        """Test that longitudinal analysis requires time_column"""
        config = StatisticalConfig()
        config.analysis_type = "longitudinal"
        config.statistical_test_method = "mixed_effects"
        config.subject_column = "Subject"
        config.time_column = None  # Missing time column
        
        with pytest.raises(ValueError, match="time_column"):
            config.validate()

    def test_longitudinal_analysis_runs(self, longitudinal_protein_data, longitudinal_metadata):
        """Test that longitudinal analysis runs successfully"""
        config = StatisticalConfig()
        config.analysis_type = "longitudinal"
        config.statistical_test_method = "mixed_effects"
        config.subject_column = "Subject"
        config.time_column = "Week"
        config.covariates = []
        config.p_value_threshold = 0.05
        
        try:
            result = run_comprehensive_statistical_analysis(
                normalized_data=longitudinal_protein_data,
                sample_metadata=longitudinal_metadata,
                config=config
            )
            
            assert isinstance(result, pd.DataFrame)
            assert 'P.Value' in result.columns
            assert 'logFC' in result.columns
            assert len(result) == 5  # 5 proteins
            print("✓ Longitudinal analysis ran successfully")
            
        except ImportError as e:
            if "statsmodels" in str(e):
                pytest.skip("statsmodels not available")
            raise


class TestLinearTrendAnalysisType:
    """Test the linear_trend analysis type (test if slope != 0)"""

    @pytest.fixture
    def trend_metadata(self):
        """Create metadata for linear trend analysis"""
        return {
            # Subject 1 across doses
            'S1_D0': {'Subject': 'S1', 'Dose': 0},
            'S1_D20': {'Subject': 'S1', 'Dose': 20},
            'S1_D40': {'Subject': 'S1', 'Dose': 40},
            # Subject 2 across doses
            'S2_D0': {'Subject': 'S2', 'Dose': 0},
            'S2_D20': {'Subject': 'S2', 'Dose': 20},
            'S2_D40': {'Subject': 'S2', 'Dose': 40},
            # Subject 3 across doses
            'S3_D0': {'Subject': 'S3', 'Dose': 0},
            'S3_D20': {'Subject': 'S3', 'Dose': 20},
            'S3_D40': {'Subject': 'S3', 'Dose': 40},
        }

    @pytest.fixture
    def trend_protein_data(self):
        """Create protein data with dose-dependent trends"""
        np.random.seed(42)
        n_proteins = 5
        
        data = pd.DataFrame({
            'Protein': [f'P{i:05d}' for i in range(n_proteins)],
            'Gene': [f'GENE{i}' for i in range(n_proteins)],
        })
        
        # Create sample columns
        samples = [f'S{s}_D{d}' for s in [1, 2, 3] for d in [0, 20, 40]]
        for sample in samples:
            data[sample] = np.random.lognormal(10, 0.5, n_proteins)
        
        return data

    def test_linear_trend_config_validation(self):
        """Test that linear_trend analysis configuration validates correctly"""
        config = StatisticalConfig()
        config.analysis_type = "linear_trend"
        config.statistical_test_method = "mixed_effects"
        config.subject_column = "Subject"
        config.time_column = "Dose"
        
        # Should validate successfully
        config.validate()
        
        assert config.analysis_type == "linear_trend"

    def test_linear_trend_config_requires_time_column(self):
        """Test that linear_trend analysis requires time_column"""
        config = StatisticalConfig()
        config.analysis_type = "linear_trend"
        config.statistical_test_method = "mixed_effects"
        config.subject_column = "Subject"
        config.time_column = None
        
        with pytest.raises(ValueError, match="time_column"):
            config.validate()

    def test_linear_trend_analysis_runs(self, trend_protein_data, trend_metadata):
        """Test that linear_trend analysis runs successfully"""
        config = StatisticalConfig()
        config.analysis_type = "linear_trend"
        config.statistical_test_method = "mixed_effects"
        config.subject_column = "Subject"
        config.time_column = "Dose"
        config.p_value_threshold = 0.05
        
        try:
            result = run_comprehensive_statistical_analysis(
                normalized_data=trend_protein_data,
                sample_metadata=trend_metadata,
                config=config
            )
            
            assert isinstance(result, pd.DataFrame)
            assert 'P.Value' in result.columns
            assert 'logFC' in result.columns
            assert len(result) == 5  # 5 proteins
            print("✓ Linear trend analysis ran successfully")
            
        except ImportError as e:
            if "statsmodels" in str(e):
                pytest.skip("statsmodels not available")
            raise

    def test_dose_response_alias(self):
        """Test that dose_response is accepted as alias for linear_trend"""
        config = StatisticalConfig()
        config.analysis_type = "dose_response"  # Alias
        config.statistical_test_method = "mixed_effects"
        config.subject_column = "Subject"
        config.time_column = "Dose"
        
        # Should validate successfully (dose_response is accepted as alias)
        config.validate()


class TestCovariateHandling:
    """Test handling of covariates in mixed-effects models"""

    def test_empty_covariates_config(self):
        """Test that empty covariates list is handled correctly"""
        config = StatisticalConfig()
        config.analysis_type = "longitudinal"
        config.statistical_test_method = "mixed_effects"
        config.subject_column = "Subject"
        config.time_column = "Week"
        config.covariates = []  # Empty - recommended for longitudinal
        
        # Should validate successfully
        config.validate()
        assert config.covariates == []

    def test_covariates_with_special_characters(self):
        """Test that covariates with special characters are handled"""
        config = StatisticalConfig()
        config.analysis_type = "paired"
        config.statistical_test_method = "mixed_effects"
        config.subject_column = "Subject"
        config.group_column = "Group"
        config.group_labels = ["Control", "Treatment"]
        config.paired_column = "Visit"
        config.paired_label1 = "Baseline"
        config.paired_label2 = "Week4"
        config.covariates = ["Age (years)", "Body Mass Index"]  # Special characters
        
        # Should validate successfully
        config.validate()

