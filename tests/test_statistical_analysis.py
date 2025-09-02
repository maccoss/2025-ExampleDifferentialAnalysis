"""
Tests for proteomics_toolkit.statistical_analysis module
"""
import pandas as pd
import numpy as np
import warnings

from proteomics_toolkit.statistical_analysis import (
    StatisticalConfig,
    prepare_metadata_dataframe,
    run_paired_t_test,
    run_unpaired_t_test,
    apply_multiple_testing_correction,
    run_comprehensive_statistical_analysis,
    display_analysis_summary,
    _create_empty_result,
    _create_empty_mixed_effects_result,
    _apply_log_transformation_if_needed
)


class TestStatisticalConfig:
    """Test the StatisticalConfig class"""
    
    def test_config_initialization(self):
        """Test configuration initialization with defaults"""
        config = StatisticalConfig()
        
        assert config.statistical_test_method == "mixed_effects"
        assert config.analysis_type == "interaction_analysis"
        assert config.p_value_threshold == 0.05
        assert config.fold_change_threshold == 1.5
        assert config.use_adjusted_pvalue == "fdr_bh"
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
    
    def test_prepare_metadata_basic(self, sample_metadata, sample_columns, statistical_config):
        """Test basic metadata preparation"""
        result = prepare_metadata_dataframe(sample_metadata, sample_columns, statistical_config)
        
        assert isinstance(result, pd.DataFrame)
        assert 'Sample' in result.columns
        assert statistical_config.subject_column in result.columns
        assert statistical_config.group_column in result.columns
        assert len(result) <= len(sample_columns)  # May filter out some samples
        
    def test_prepare_metadata_missing_columns(self, sample_metadata, sample_columns):
        """Test metadata preparation when required columns are missing"""
        config = StatisticalConfig()
        config.subject_column = "NonexistentColumn"
        
        # This should handle missing columns gracefully
        result = prepare_metadata_dataframe(sample_metadata, sample_columns, config)
        assert isinstance(result, pd.DataFrame)


class TestRunPairedTTest:
    """Test paired t-test analysis"""
    
    def test_paired_t_test_basic(self, sample_protein_data, sample_metadata_df, statistical_config):
        """Test basic paired t-test"""
        # Configure for paired analysis
        statistical_config.statistical_test_method = "paired_t"
        
        # Select only sample columns for analysis
        sample_columns = [col for col in sample_protein_data.columns if col.startswith('Sample_')]
        protein_data = sample_protein_data[sample_columns]
        
        result = run_paired_t_test(protein_data, sample_metadata_df, statistical_config)
        
        assert isinstance(result, pd.DataFrame)
        assert 'P.Value' in result.columns
        assert 'logFC' in result.columns
        assert 't' in result.columns
        assert len(result) == len(protein_data)  # One result per protein
        
    def test_paired_t_test_insufficient_data(self, statistical_config):
        """Test paired t-test with insufficient data"""
        # Create minimal dataset
        protein_data = pd.DataFrame({
            'Sample_1': [100.0, 200.0],
            'Sample_2': [110.0, 210.0]
        })
        
        metadata_df = pd.DataFrame({
            'Sample': ['Sample_1', 'Sample_2'],
            'Subject': ['S001', 'S001'],  # Same subject for pairing
            'Group': ['Control', 'Control'],
            'Visit': ['Baseline', 'Week4']
        })
        
        statistical_config.statistical_test_method = "paired_t"
        
        result = run_paired_t_test(protein_data, metadata_df, statistical_config)
        
        assert isinstance(result, pd.DataFrame)
        assert len(result) == len(protein_data)


class TestRunUnpairedTTest:
    """Test unpaired t-test analysis"""
    
    def test_unpaired_t_test_basic(self, sample_protein_data, sample_metadata_df, statistical_config):
        """Test basic unpaired t-test"""
        statistical_config.statistical_test_method = "welch_t"
        
        # Select only sample columns for analysis
        sample_columns = [col for col in sample_protein_data.columns if col.startswith('Sample_')]
        protein_data = sample_protein_data[sample_columns]
        
        result = run_unpaired_t_test(protein_data, sample_metadata_df, statistical_config)
        
        assert isinstance(result, pd.DataFrame)
        assert 'P.Value' in result.columns
        assert 'logFC' in result.columns
        assert 't' in result.columns
        assert len(result) == len(protein_data)


class TestMultipleTestingCorrection:
    """Test multiple testing correction"""
    
    def test_fdr_correction(self, differential_results, statistical_config):
        """Test FDR correction"""
        statistical_config.use_adjusted_pvalue = "fdr_bh"
        
        result = apply_multiple_testing_correction(differential_results, statistical_config)
        
        assert 'adj.P.Val' in result.columns
        assert all(result['adj.P.Val'] >= result['P.Value'])  # Adjusted p-values should be >= original
        
    def test_bonferroni_correction(self, differential_results, statistical_config):
        """Test Bonferroni correction"""
        statistical_config.use_adjusted_pvalue = "bonferroni"
        
        result = apply_multiple_testing_correction(differential_results, statistical_config)
        
        assert 'adj.P.Val' in result.columns
        assert all(result['adj.P.Val'] >= result['P.Value'])
        
    def test_no_correction(self, differential_results, statistical_config):
        """Test when no correction is applied"""
        statistical_config.use_adjusted_pvalue = "none"
        
        result = apply_multiple_testing_correction(differential_results, statistical_config)
        
        assert 'adj.P.Val' in result.columns
        assert all(result['adj.P.Val'] == result['P.Value'])  # Should be identical


class TestRunComprehensiveStatisticalAnalysis:
    """Test the main comprehensive analysis function"""
    
    def test_comprehensive_analysis_paired_t(self, sample_protein_data, sample_metadata, statistical_config):
        """Test comprehensive analysis with paired t-test"""
        statistical_config.statistical_test_method = "paired_t"
        
        # Suppress output during testing
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            result = run_comprehensive_statistical_analysis(
                normalized_data=sample_protein_data,
                sample_metadata=sample_metadata,
                config=statistical_config,
                protein_annotations=None
            )
        
        assert isinstance(result, pd.DataFrame)
        assert 'P.Value' in result.columns
        assert 'adj.P.Val' in result.columns
        assert len(result) > 0
        
    def test_comprehensive_analysis_with_annotations(self, sample_protein_data, sample_metadata, statistical_config):
        """Test comprehensive analysis with protein annotations"""
        statistical_config.statistical_test_method = "paired_t"
        
        # Use the protein data itself as annotations (contains Protein, Gene columns)
        annotations = sample_protein_data[['Protein', 'ProteinName', 'Gene']].copy()
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            result = run_comprehensive_statistical_analysis(
                normalized_data=sample_protein_data,
                sample_metadata=sample_metadata,
                config=statistical_config,
                protein_annotations=annotations
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
                label_top_n=5
            )
        
        assert isinstance(summary, dict)
        assert 'analysis_method' in summary
        assert 'total_proteins' in summary
        assert 'significant_proteins' in summary
        
    def test_display_summary_empty_results(self, statistical_config):
        """Test summary display with empty results"""
        empty_results = pd.DataFrame()
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            summary = display_analysis_summary(
                differential_results=empty_results,
                config=statistical_config
            )
        
        assert isinstance(summary, dict)


class TestHelperFunctions:
    """Test helper functions"""
    
    def test_create_empty_result(self):
        """Test creating empty result for failed analysis"""
        result = _create_empty_result("P00001", "Test reason")
        
        assert isinstance(result, dict)
        assert result['Protein'] == "P00001"
        assert 'P.Value' in result
        assert pd.isna(result['P.Value'])
        
    def test_create_empty_mixed_effects_result(self):
        """Test creating empty mixed-effects result"""
        result = _create_empty_mixed_effects_result("P00001", "Test reason")
        
        assert isinstance(result, dict)
        assert result['Protein'] == "P00001"
        assert 'P.Value' in result
        assert 'group_effect' in result
        assert pd.isna(result['P.Value'])


class TestEdgeCases:
    """Test edge cases and error conditions"""
    
    def test_analysis_with_all_nan_protein(self, sample_metadata, statistical_config):
        """Test analysis when a protein has all NaN values"""
        # Create protein data with one all-NaN protein
        protein_data = pd.DataFrame({
            'Protein': ['P00001', 'P00002'],
            'Sample_A_1': [100.0, np.nan],
            'Sample_A_2': [110.0, np.nan],
            'Sample_B_1': [120.0, np.nan],
            'Sample_B_2': [115.0, np.nan]
        })
        
        statistical_config.statistical_test_method = "paired_t"
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            result = run_comprehensive_statistical_analysis(
                normalized_data=protein_data,
                sample_metadata=sample_metadata,
                config=statistical_config
            )
        
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 2  # Should handle both proteins
        
    def test_analysis_with_no_valid_samples(self, statistical_config):
        """Test analysis when no samples match metadata"""
        protein_data = pd.DataFrame({
            'Protein': ['P00001'],
            'Unknown_Sample': [100.0]
        })
        
        sample_metadata = {
            'Different_Sample': {'Subject': 'S001', 'Group': 'Control'}
        }
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            result = run_comprehensive_statistical_analysis(
                normalized_data=protein_data,
                sample_metadata=sample_metadata,
                config=statistical_config
            )
        
        assert isinstance(result, pd.DataFrame)
        # Should return empty results or handle gracefully


class TestLogTransformation:
    """Test log transformation functionality"""
    
    def test_statistical_config_has_log_parameters(self):
        """Test that StatisticalConfig has all required log transformation parameters"""
        config = StatisticalConfig()
        
        # Check that all log transformation parameters exist
        assert hasattr(config, 'log_transform_before_stats')
        assert hasattr(config, 'log_base') 
        assert hasattr(config, 'log_pseudocount')
        
        # Check default values
        assert config.log_transform_before_stats == "auto"
        assert config.log_base == "log2"
        assert config.log_pseudocount is None
        
    def test_log_transformation_auto_mode_with_median_normalization(self):
        """Test that auto mode applies log transformation for median normalization"""
        config = StatisticalConfig()
        config.log_transform_before_stats = "auto"
        
        # Create test data with large values typical of median normalization
        test_data = pd.DataFrame({
            'protein1': [1000.0, 1200.0, 950.0, 1100.0],
            'protein2': [800.0, 900.0, 750.0, 850.0],
            'protein3': [1500.0, 1800.0, 1400.0, 1600.0]
        })
        test_data.attrs = {'normalization_method': 'median'}
        
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
        test_data = pd.DataFrame({
            'protein1': [12.5, 13.2, 11.8, 12.9],
            'protein2': [10.8, 11.5, 10.2, 11.1],
            'protein3': [14.1, 14.8, 13.5, 14.2]
        })
        test_data.attrs = {'normalization_method': 'vsn'}
        
        result = _apply_log_transformation_if_needed(test_data, config)
        
        # Data should NOT be transformed (VSN is already on appropriate scale)
        assert np.array_equal(test_data.values, result.values)
        
    def test_log_transformation_forced_true(self):
        """Test forced log transformation (log_transform_before_stats=True)"""
        config = StatisticalConfig()
        config.log_transform_before_stats = "true"  # Use string instead of boolean
        config.log_base = "log2"
        
        test_data = pd.DataFrame({
            'protein1': [1000.0, 2000.0, 500.0, 1500.0],
            'protein2': [800.0, 1600.0, 400.0, 1200.0]
        })
        
        result = _apply_log_transformation_if_needed(test_data, config)
        
        # Should be transformed regardless of normalization method
        assert not np.array_equal(test_data.values, result.values)
        assert result.mean().mean() < test_data.mean().mean()
        
    def test_log_transformation_forced_false(self):
        """Test that transformation is skipped when log_transform_before_stats=False"""
        config = StatisticalConfig()
        config.log_transform_before_stats = "false"  # Use string instead of boolean
        
        test_data = pd.DataFrame({
            'protein1': [1000.0, 2000.0, 500.0, 1500.0],
            'protein2': [800.0, 1600.0, 400.0, 1200.0]
        })
        
        result = _apply_log_transformation_if_needed(test_data, config)
        
        # Data should NOT be transformed
        assert np.array_equal(test_data.values, result.values)
        
    def test_log_transformation_different_bases(self):
        """Test different logarithmic bases"""
        test_data = pd.DataFrame({
            'protein1': [1000.0, 2000.0, 500.0],
            'protein2': [800.0, 1600.0, 400.0]
        })
        
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
        
        test_data = pd.DataFrame({
            'protein1': [100.0, 200.0, 50.0, 150.0],
            'protein2': [80.0, 160.0, 40.0, 120.0]
        })
        
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
        
        test_data = pd.DataFrame({
            'protein1': [1000.0, 2000.0, 500.0, 1500.0],
            'protein2': [800.0, 1600.0, 400.0, 1200.0]
        })
        
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
        
        test_data = pd.DataFrame({
            'protein1': [0.0, 100.0, 200.0, 150.0],  # Contains zero
            'protein2': [50.0, 150.0, 0.0, 100.0]   # Contains zero
        })
        
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
        
        test_data = pd.DataFrame({
            'protein1': [1000.0, 1200.0, 950.0],
            'protein2': [800.0, 900.0, 750.0],
            'protein3': [1500.0, 1800.0, 1400.0]
        }, index=['sample1', 'sample2', 'sample3'])
        
        # Add some attributes
        test_data.attrs = {'normalization_method': 'median', 'test_attr': 'test_value'}
        
        result = _apply_log_transformation_if_needed(test_data, config)
        
        # Check structure preservation
        assert result.shape == test_data.shape
        assert list(result.columns) == list(test_data.columns)
        assert list(result.index) == list(test_data.index)
        
        # Check that attrs are preserved
        assert result.attrs.get('normalization_method') == 'median'
        assert result.attrs.get('test_attr') == 'test_value'
        
    def test_log_transformation_invalid_base_raises_error(self):
        """Test that invalid log base raises appropriate error"""
        config = StatisticalConfig()
        config.log_transform_before_stats = "true"  # Use string
        config.log_base = "invalid_base"
        
        test_data = pd.DataFrame({
            'protein1': [1000.0, 1200.0, 950.0],
            'protein2': [800.0, 900.0, 750.0]
        })
        
        try:
            _apply_log_transformation_if_needed(test_data, config)
            assert False, "Should have raised ValueError for invalid log base"
        except ValueError as e:
            assert "Unknown log base" in str(e)
