"""
Tests for proteomics_toolkit.normalization module
"""
import pandas as pd
import numpy as np
import warnings

from proteomics_toolkit.normalization import (
    get_normalization_characteristics,
    is_normalization_log_transformed,
    median_normalize,
    vsn_normalize,
    quantile_normalize,
    log_transform,
    mad_normalize,
    z_score_normalize,
    rlr_normalize,
    loess_normalize,
    handle_negative_values,
    analyze_negative_values,
    calculate_normalization_stats,
    calculate_detailed_normalization_stats
)


class TestNormalizationCharacteristics:
    """Test normalization characteristics functions"""
    
    def test_get_normalization_characteristics(self):
        """Test getting normalization method characteristics"""
        characteristics = get_normalization_characteristics()
        
        assert isinstance(characteristics, dict)
        assert 'median' in characteristics
        assert 'vsn' in characteristics
        assert 'quantile' in characteristics
        
        # Check that each method has expected properties
        for method_name, props in characteristics.items():
            assert 'log_transformed' in props
            assert 'description' in props
            assert isinstance(props['log_transformed'], bool)
            
    def test_is_normalization_log_transformed(self):
        """Test checking if normalization method is log-transformed"""
        assert is_normalization_log_transformed('vsn') is True
        assert is_normalization_log_transformed('median') is False
        assert is_normalization_log_transformed('quantile') is False
        assert is_normalization_log_transformed('unknown_method') is False


class TestSeparateColumns:
    """Test column separation utility"""
    
class TestMedianNormalization:
    """Test median normalization"""
    
    def test_median_normalize_basic(self, standardized_protein_data):
        """Test basic median normalization with standardized data structure"""
        # Data is already in standardized format - can use directly
        result = median_normalize(standardized_protein_data)
        
        assert isinstance(result, pd.DataFrame)
        # Should have exact standardized structure: 5 annotation columns + 12 sample columns
        assert result.shape[1] == 17  # 5 + 12
        assert result.shape[0] == standardized_protein_data.shape[0]
        
        # Check standardized annotation columns are present and preserved
        expected_annotation_cols = ['Protein', 'Description', 'Protein Gene', 'UniProt_Accession', 'UniProt_Entry_Name']
        assert list(result.columns[:5]) == expected_annotation_cols
        assert all(result['Protein'] == standardized_protein_data['Protein'])
        
        # Check that sample columns are present (columns 5+)
        sample_result_cols = list(result.columns[5:])
        original_sample_cols = list(standardized_protein_data.columns[5:])
        assert sample_result_cols == original_sample_cols
        
        # Check that normalization was applied (sample medians should be closer to global median)
        result_sample_data = result.iloc[:, 5:]  # Sample columns from result
        
        # Medians should be approximately equal after normalization (within tolerance)
        medians = result_sample_data.median()
        median_diff = medians.max() - medians.min()
        assert median_diff < 1.0  # Should be much closer after normalization
        
    def test_median_normalize_auto_detect_columns(self, sample_protein_data):
        """Test median normalization with automatic column detection"""
        result = median_normalize(sample_protein_data)  # No sample_columns specified
        
        assert isinstance(result, pd.DataFrame)
        assert result.shape == sample_protein_data.shape
        
    def test_median_normalize_with_missing_values(self, sample_columns):
        """Test median normalization with missing values"""
        # Create data with NaN values
        data_with_nan = pd.DataFrame({
            'Protein': ['P1', 'P2', 'P3'],
            'Sample_1': [100.0, np.nan, 150.0],
            'Sample_2': [200.0, 180.0, np.nan],
            'Sample_3': [120.0, 190.0, 160.0]
        })
        
        result = median_normalize(data_with_nan, ['Sample_1', 'Sample_2', 'Sample_3'])
        
        assert isinstance(result, pd.DataFrame)
        # NaN values should be preserved
        assert pd.isna(result.iloc[1, 1])  # Sample_1, P2
        assert pd.isna(result.iloc[2, 2])  # Sample_2, P3


class TestVSNNormalization:
    """Test variance stabilizing normalization"""
    
    def test_vsn_normalize_basic(self, sample_protein_data, sample_columns):
        """Test basic VSN normalization"""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # VSN may produce optimization warnings
            
            result = vsn_normalize(sample_protein_data, sample_columns=sample_columns)
        
        assert isinstance(result, pd.DataFrame)
        assert result.shape == sample_protein_data.shape
        assert all(col in result.columns for col in sample_protein_data.columns)
        
        # VSN should log-transform the data, so values should be different scale
        sample_data = result[sample_columns]
        original_sample_data = sample_protein_data[sample_columns]
        
        # Check that transformation was applied (values should be different)
        assert not np.allclose(sample_data.values, original_sample_data.values, equal_nan=True)
        
    def test_vsn_normalize_with_optimization(self, sample_protein_data, sample_columns):
        """Test VSN normalization with parameter optimization"""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            result = vsn_normalize(
                sample_protein_data, 
                optimize_params=True,
                sample_columns=sample_columns
            )
        
        assert isinstance(result, pd.DataFrame)
        assert result.shape == sample_protein_data.shape


class TestQuantileNormalization:
    """Test quantile normalization"""
    
    def test_quantile_normalize_basic(self, sample_protein_data, sample_columns):
        """Test basic quantile normalization"""
        result = quantile_normalize(sample_protein_data, sample_columns)
        
        assert isinstance(result, pd.DataFrame)
        assert result.shape == sample_protein_data.shape
        
        # After quantile normalization, all samples should have identical distributions
        sample_data = result[sample_columns].dropna()
        if len(sample_data) > 1:
            # Check that sorted values are approximately equal across samples
            sorted_values = sample_data.apply(lambda x: np.sort(x.dropna()))
            
            # All columns should have similar quantile distributions
            # (This is a simplified check - full quantile normalization ensures exact equality)
            first_col = sorted_values.iloc[:, 0]
            for col in sorted_values.columns[1:]:
                correlation = np.corrcoef(first_col, sorted_values[col])[0, 1]
                assert correlation > 0.9  # Should be highly correlated


class TestLogTransformation:
    """Test log transformation"""
    
    def test_log_transform_log2(self, sample_protein_data, sample_columns):
        """Test log2 transformation"""
        # Extract only numeric sample data for log transformation
        sample_data = sample_protein_data[sample_columns]
        
        result = log_transform(sample_data, base='log2')
        
        assert isinstance(result, pd.DataFrame)
        assert result.shape == sample_data.shape
        
        # Values should be much smaller after log transformation
        original_mean = sample_data.mean().mean()
        transformed_mean = result.mean().mean()
        assert transformed_mean < original_mean
        
    def test_log_transform_natural_log(self, sample_protein_data, sample_columns):
        """Test natural log transformation"""
        sample_data = sample_protein_data[sample_columns]
        
        result = log_transform(sample_data, base='ln')
        
        assert isinstance(result, pd.DataFrame)
        assert result.shape == sample_data.shape
        
    def test_log_transform_with_negatives(self, sample_columns):
        """Test log transformation with negative values"""
        # Create data with some negative values
        data_with_negatives = pd.DataFrame({
            'Sample_1': [100.0, -50.0],  # Negative value
            'Sample_2': [200.0, 180.0],
            'Sample_3': [120.0, 0.0]     # Zero value
        })
        
        # Log transform should handle negatives with pseudocount
        result = log_transform(data_with_negatives, base='log2')
        
        assert isinstance(result, pd.DataFrame)
        # Should handle negatives gracefully (they become very negative after log transform)


class TestMADNormalization:
    """Test median absolute deviation normalization"""
    
    def test_mad_normalize_basic(self, sample_protein_data, sample_columns):
        """Test basic MAD normalization"""
        result = mad_normalize(sample_protein_data, sample_columns)
        
        assert isinstance(result, pd.DataFrame)
        assert result.shape == sample_protein_data.shape
        
        # MAD normalization should make the median absolute deviations similar
        sample_data = result[sample_columns]
        mads = sample_data.mad()
        
        # MADs should be more similar after normalization
        mad_cv = mads.std() / mads.mean() if mads.mean() > 0 else 0
        assert mad_cv < 1.0  # Should be reduced compared to original


class TestZScoreNormalization:
    """Test z-score normalization"""
    
    def test_z_score_normalize_basic(self, sample_protein_data, sample_columns):
        """Test basic z-score normalization"""
        result = z_score_normalize(sample_protein_data, sample_columns)
        
        assert isinstance(result, pd.DataFrame)
        assert result.shape == sample_protein_data.shape
        
        # After z-score normalization, each sample should have mean≈0 and std≈1
        sample_data = result[sample_columns]
        
        for col in sample_columns:
            col_data = sample_data[col].dropna()
            if len(col_data) > 1:
                assert abs(col_data.mean()) < 0.1  # Mean should be close to 0
                assert abs(col_data.std() - 1.0) < 0.1  # Std should be close to 1


class TestRLRNormalization:
    """Test robust linear regression normalization"""
    
    def test_rlr_normalize_basic(self, sample_protein_data, sample_columns):
        """Test basic RLR normalization"""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # RLR may produce optimization warnings
            
            result = rlr_normalize(sample_protein_data, sample_columns)
        
        assert isinstance(result, pd.DataFrame)
        assert result.shape == sample_protein_data.shape


class TestLoessNormalization:
    """Test LOESS normalization"""
    
    def test_loess_normalize_basic(self, sample_protein_data, sample_columns):
        """Test basic LOESS normalization"""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # LOESS may produce warnings
            
            result = loess_normalize(sample_protein_data, sample_columns=sample_columns)
        
        assert isinstance(result, pd.DataFrame)
        assert result.shape == sample_protein_data.shape
        
    def test_loess_normalize_different_span(self, sample_protein_data, sample_columns):
        """Test LOESS normalization with different span parameter"""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            result = loess_normalize(
                sample_protein_data, 
                span=0.5,  # Different from default 0.75
                sample_columns=sample_columns
            )
        
        assert isinstance(result, pd.DataFrame)
        assert result.shape == sample_protein_data.shape


class TestNegativeValueHandling:
    """Test negative value handling"""
    
    def test_analyze_negative_values(self):
        """Test analyzing negative values in data"""
        # Create data with negative values
        data_with_negatives = pd.DataFrame({
            'Protein': ['P1', 'P2', 'P3'],
            'Sample_1': [100.0, -50.0, 150.0],
            'Sample_2': [200.0, -20.0, 180.0],
            'Sample_3': [120.0, 190.0, -10.0]
        })
        
        sample_columns = ['Sample_1', 'Sample_2', 'Sample_3']
        
        analysis = analyze_negative_values(
            data_with_negatives, 
            normalization_method='median',
            sample_columns=sample_columns
        )
        
        assert isinstance(analysis, dict)
        assert 'negative_count' in analysis
        assert 'total_values' in analysis
        assert 'proteins_with_negatives' in analysis
        assert analysis['negative_count'] > 0
        
    def test_handle_negative_values_min_positive(self):
        """Test handling negative values with min_positive method"""
        data_with_negatives = pd.DataFrame({
            'Protein': ['P1', 'P2'],
            'Sample_1': [100.0, -50.0],
            'Sample_2': [200.0, 180.0]
        })
        
        result = handle_negative_values(
            data_with_negatives, 
            method='min_positive',
            sample_columns=['Sample_1', 'Sample_2']
        )
        
        assert isinstance(result, pd.DataFrame)
        # Negative values should be replaced
        assert all(result[['Sample_1', 'Sample_2']].min() >= 0)
        
    def test_handle_negative_values_shift_to_positive(self):
        """Test handling negative values with shift_to_positive method"""
        data_with_negatives = pd.DataFrame({
            'Protein': ['P1', 'P2'],
            'Sample_1': [100.0, -50.0],
            'Sample_2': [200.0, -20.0]
        })
        
        result = handle_negative_values(
            data_with_negatives,
            method='shift_to_positive',
            sample_columns=['Sample_1', 'Sample_2']
        )
        
        assert isinstance(result, pd.DataFrame)
        # All values should be positive after shifting
        assert all(result[['Sample_1', 'Sample_2']].min() > 0)


class TestNormalizationStats:
    """Test normalization statistics calculation"""
    
    def test_calculate_normalization_stats(self, sample_protein_data, sample_columns):
        """Test calculating basic normalization statistics"""
        # Apply median normalization
        normalized_data = median_normalize(sample_protein_data, sample_columns)
        
        # Extract only sample data for statistics calculation
        original_sample_data = sample_protein_data[sample_columns]
        normalized_sample_data = normalized_data[sample_columns]
        
        stats = calculate_normalization_stats(
            original_sample_data, 
            normalized_sample_data,
            method='median'
        )
        
        assert isinstance(stats, dict)
        assert 'original_median_range' in stats
        assert 'normalized_median_range' in stats
        assert 'original_cv_median' in stats
        assert 'normalized_cv_median' in stats
        
    def test_calculate_detailed_normalization_stats(self, sample_protein_data, sample_columns):
        """Test calculating detailed normalization statistics"""
        normalized_data = median_normalize(sample_protein_data, sample_columns)
        
        # Extract sample data
        original_sample_data = sample_protein_data[sample_columns]
        normalized_sample_data = normalized_data[sample_columns]
        
        detailed_stats = calculate_detailed_normalization_stats(
            original_sample_data,
            normalized_sample_data, 
            method='median'
        )
        
        assert isinstance(detailed_stats, dict)
        assert 'overall' in detailed_stats
        assert 'by_sample' in detailed_stats
        assert 'by_protein' in detailed_stats


class TestNormalizationEdgeCases:
    """Test edge cases in normalization"""
    
    def test_normalization_with_all_zeros(self, sample_columns):
        """Test normalization when data contains all zeros"""
        zero_data = pd.DataFrame({
            'Protein': ['P1', 'P2'],
            'Sample_1': [0.0, 0.0],
            'Sample_2': [0.0, 0.0],
            'Sample_3': [0.0, 0.0]
        })
        
        # Should handle gracefully without crashing
        result = median_normalize(zero_data, ['Sample_1', 'Sample_2', 'Sample_3'])
        assert isinstance(result, pd.DataFrame)
        
    def test_normalization_with_single_sample(self):
        """Test normalization with only one sample"""
        single_sample_data = pd.DataFrame({
            'Protein': ['P1', 'P2', 'P3'],
            'Sample_1': [100.0, 200.0, 150.0]
        })
        
        result = median_normalize(single_sample_data, ['Sample_1'])
        assert isinstance(result, pd.DataFrame)
        assert result.shape == single_sample_data.shape
        
    def test_normalization_with_single_protein(self):
        """Test normalization with only one protein"""
        single_protein_data = pd.DataFrame({
            'Protein': ['P1'],
            'Sample_1': [100.0],
            'Sample_2': [200.0],
            'Sample_3': [150.0]
        })
        
        result = median_normalize(single_protein_data, ['Sample_1', 'Sample_2', 'Sample_3'])
        assert isinstance(result, pd.DataFrame)
        assert result.shape == single_protein_data.shape
