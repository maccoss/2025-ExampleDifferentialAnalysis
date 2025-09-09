"""
Tests for CV distribution analysis and control sample quality assessment
"""

import pandas as pd
import numpy as np
import pytest
import matplotlib.pyplot as plt
from unittest.mock import patch
from proteomics_toolkit.visualization import plot_control_cv_distribution
from proteomics_toolkit.normalization import median_normalize


class TestCVDistributionAnalysis:
    """Test CV distribution analysis for control samples"""

    @pytest.fixture
    def mock_cv_data_setup(self):
        """Create mock data for CV analysis testing"""
        # Create sample data with known CV characteristics
        np.random.seed(42)
        
        # Control samples with different CV characteristics
        control_samples = {
            'HoofPool': ['HoofPool_1', 'HoofPool_2', 'HoofPool_3', 'HoofPool_4'],
            'GWPool': ['GWPool_1', 'GWPool_2', 'GWPool_3', 'GWPool_4'],
            'PlatePool': ['PlatePool_1', 'PlatePool_2', 'PlatePool_3', 'PlatePool_4']
        }
        
        all_control_samples = []
        for samples in control_samples.values():
            all_control_samples.extend(samples)
        
        # Create protein data with 100 proteins
        protein_data = pd.DataFrame({
            'Protein': [f'P{i:03d}' for i in range(100)],
            'Description': [f'Protein {i} description' for i in range(100)],
            'Protein Gene': [f'GENE{i}' for i in range(100)],
            'UniProt_Accession': [f'P{i:05d}' for i in range(100)],
            'UniProt_Entry_Name': [f'PROT{i}_HUMAN' for i in range(100)],
        })
        
        # Add control sample data with different CV patterns
        for i, sample_name in enumerate(all_control_samples):
            if 'HoofPool' in sample_name:
                # Low CV controls (good reproducibility)
                base_values = np.random.normal(1000, 100, 100)
                noise = np.random.normal(0, 50, 100)  # 5% CV
            elif 'GWPool' in sample_name:
                # Medium CV controls
                base_values = np.random.normal(800, 150, 100)
                noise = np.random.normal(0, 120, 100)  # 15% CV
            else:  # PlatePool
                # Higher CV controls
                base_values = np.random.normal(1200, 200, 100)
                noise = np.random.normal(0, 240, 100)  # 20% CV
            
            protein_data[sample_name] = np.maximum(base_values + noise, 1)  # Ensure positive values
        
        # Create sample metadata
        sample_metadata = {}
        for control_type, samples in control_samples.items():
            for sample in samples:
                sample_metadata[sample] = {
                    'Subject': control_type,
                    'SampleType': 'Control',
                    'ControlType': control_type
                }
        
        return protein_data, all_control_samples, sample_metadata, control_samples

    def test_cv_distribution_basic_functionality(self, mock_cv_data_setup):
        """Test basic CV distribution plotting functionality"""
        protein_data, sample_columns, sample_metadata, control_samples = mock_cv_data_setup
        control_labels = list(control_samples.keys())
        
        with patch('matplotlib.pyplot.show'):
            cv_data = plot_control_cv_distribution(
                data=protein_data,
                sample_columns=sample_columns,
                sample_metadata=sample_metadata,
                control_column='Subject',
                control_labels=control_labels,
                normalization_method="Test",
                figsize=(15, 5),
                cv_threshold=20.0
            )
        
        # Verify return data structure
        assert isinstance(cv_data, dict)
        assert len(cv_data) == len(control_labels)
        
        # Verify each control type has CV data
        for control_type in control_labels:
            assert control_type in cv_data
            assert isinstance(cv_data[control_type], list)
            assert len(cv_data[control_type]) > 0  # Should have CV values
        
        plt.close('all')

    def test_cv_calculation_accuracy(self, mock_cv_data_setup):
        """Test that CV calculations are mathematically correct"""
        protein_data, sample_columns, sample_metadata, control_samples = mock_cv_data_setup
        
        # Test CV calculation for HoofPool samples (should have low CV)
        hoof_samples = control_samples['HoofPool']
        hoof_data = protein_data[hoof_samples]
        
        # Calculate CV manually for first protein
        first_protein_values = hoof_data.iloc[0].values
        manual_cv = (np.std(first_protein_values) / np.mean(first_protein_values)) * 100
        
        with patch('matplotlib.pyplot.show'):
            cv_data = plot_control_cv_distribution(
                data=protein_data,
                sample_columns=sample_columns,
                sample_metadata=sample_metadata,
                control_column='Subject',
                control_labels=['HoofPool'],
                normalization_method="Test"
            )
        
        # The first CV value should match our manual calculation (within reasonable tolerance)
        calculated_cv = cv_data['HoofPool'][0]
        assert abs(calculated_cv - manual_cv) < 2.0  # Allow for reasonable difference due to implementation details
        
        plt.close('all')

    def test_cv_distribution_with_insufficient_samples(self):
        """Test CV analysis behavior with insufficient samples per control type"""
        # Create data with only 1 sample per control type
        protein_data = pd.DataFrame({
            'Protein': ['P001', 'P002'],
            'Description': ['Protein 1', 'Protein 2'],
            'Protein Gene': ['GENE1', 'GENE2'],
            'UniProt_Accession': ['P00001', 'P00002'],
            'UniProt_Entry_Name': ['PROT1_HUMAN', 'PROT2_HUMAN'],
            'Control_A_1': [100, 200],
            'Control_B_1': [150, 250]
        })
        
        sample_columns = ['Control_A_1', 'Control_B_1']
        sample_metadata = {
            'Control_A_1': {'Subject': 'ControlA'},
            'Control_B_1': {'Subject': 'ControlB'}
        }
        
        with patch('matplotlib.pyplot.show'):
            cv_data = plot_control_cv_distribution(
                data=protein_data,
                sample_columns=sample_columns,
                sample_metadata=sample_metadata,
                control_column='Subject',
                control_labels=['ControlA', 'ControlB'],
                normalization_method="Test"
            )
        
        # Should return empty lists for controls with insufficient samples
        assert cv_data['ControlA'] == []
        assert cv_data['ControlB'] == []
        
        plt.close('all')

    def test_cv_distribution_with_missing_controls(self, mock_cv_data_setup):
        """Test CV analysis when some control types are missing from data"""
        protein_data, sample_columns, sample_metadata, control_samples = mock_cv_data_setup
        
        # Request a control type that doesn't exist in the data
        control_labels = ['HoofPool', 'NonexistentPool']
        
        with patch('matplotlib.pyplot.show'):
            cv_data = plot_control_cv_distribution(
                data=protein_data,
                sample_columns=sample_columns,
                sample_metadata=sample_metadata,
                control_column='Subject',
                control_labels=control_labels,
                normalization_method="Test"
            )
        
        # Should handle missing control gracefully
        assert 'HoofPool' in cv_data
        assert len(cv_data['HoofPool']) > 0
        assert 'NonexistentPool' in cv_data
        assert cv_data['NonexistentPool'] == []
        
        plt.close('all')

    def test_cv_distribution_zero_mean_handling(self):
        """Test CV calculation when protein has zero mean (division by zero)"""
        protein_data = pd.DataFrame({
            'Protein': ['P001', 'P002'],
            'Description': ['Protein 1', 'Protein 2'],
            'Protein Gene': ['GENE1', 'GENE2'],
            'UniProt_Accession': ['P00001', 'P00002'],
            'UniProt_Entry_Name': ['PROT1_HUMAN', 'PROT2_HUMAN'],
            'Control_1': [0, 100],  # First protein has zero mean
            'Control_2': [0, 200],
            'Control_3': [0, 150],
            'Control_4': [0, 180]
        })
        
        sample_columns = ['Control_1', 'Control_2', 'Control_3', 'Control_4']
        sample_metadata = {
            sample: {'Subject': 'TestControl'} for sample in sample_columns
        }
        
        with patch('matplotlib.pyplot.show'):
            cv_data = plot_control_cv_distribution(
                data=protein_data,
                sample_columns=sample_columns,
                sample_metadata=sample_metadata,
                control_column='Subject',
                control_labels=['TestControl'],
                normalization_method="Test"
            )
        
        # Should only include one protein (the one with non-zero mean)
        assert len(cv_data['TestControl']) == 1
        
        plt.close('all')

    def test_cv_distribution_different_thresholds(self, mock_cv_data_setup):
        """Test CV analysis with different CV thresholds"""
        protein_data, sample_columns, sample_metadata, control_samples = mock_cv_data_setup
        
        with patch('matplotlib.pyplot.show'):
            # Test with 10% threshold
            cv_data_10 = plot_control_cv_distribution(
                data=protein_data,
                sample_columns=sample_columns,
                sample_metadata=sample_metadata,
                control_column='Subject',
                control_labels=['HoofPool'],
                cv_threshold=10.0
            )
            
            # Test with 30% threshold
            cv_data_30 = plot_control_cv_distribution(
                data=protein_data,
                sample_columns=sample_columns,
                sample_metadata=sample_metadata,
                control_column='Subject',
                control_labels=['HoofPool'],
                cv_threshold=30.0
            )
        
        # Data should be the same regardless of threshold
        assert cv_data_10['HoofPool'] == cv_data_30['HoofPool']
        
        plt.close('all')

    def test_cv_distribution_integration_with_normalization(self):
        """Test CV analysis integration with median normalization"""
        # Create test data
        np.random.seed(42)
        protein_data = pd.DataFrame({
            'Protein': [f'P{i:03d}' for i in range(50)],
            'Description': [f'Protein {i}' for i in range(50)],
            'Protein Gene': [f'GENE{i}' for i in range(50)],
            'UniProt_Accession': [f'P{i:05d}' for i in range(50)],
            'UniProt_Entry_Name': [f'PROT{i}_HUMAN' for i in range(50)],
        })
        
        # Add control samples with different scales
        control_samples = ['Pool_1', 'Pool_2', 'Pool_3', 'Pool_4']
        for i, sample in enumerate(control_samples):
            # Different sample scales to test normalization effect
            scale_factor = 1 + i * 0.5
            protein_data[sample] = np.random.normal(1000 * scale_factor, 100, 50)
        
        sample_metadata = {
            sample: {'Subject': 'TestPool'} for sample in control_samples
        }
        
        # Apply median normalization
        normalized_data = median_normalize(protein_data, sample_columns=control_samples)
        
        with patch('matplotlib.pyplot.show'):
            cv_data = plot_control_cv_distribution(
                data=normalized_data,
                sample_columns=control_samples,
                sample_metadata=sample_metadata,
                control_column='Subject',
                control_labels=['TestPool'],
                normalization_method="Median"
            )
        
        # Should produce valid CV data
        assert len(cv_data['TestPool']) > 0
        # CV values should be reasonable (0-100%)
        assert all(0 <= cv <= 100 for cv in cv_data['TestPool'])
        
        plt.close('all')


class TestCVDistributionPlotting:
    """Test plotting aspects of CV distribution analysis"""

    def test_plot_creation_and_cleanup(self):
        """Test that plots are created and can be cleaned up properly"""
        # Simple test data
        protein_data = pd.DataFrame({
            'Protein': ['P1', 'P2'],
            'Description': ['Protein 1', 'Protein 2'],
            'Protein Gene': ['GENE1', 'GENE2'],
            'UniProt_Accession': ['P00001', 'P00002'],
            'UniProt_Entry_Name': ['PROT1_HUMAN', 'PROT2_HUMAN'],
            'Control_1': [100, 200],
            'Control_2': [110, 190],
            'Control_3': [95, 210]
        })
        
        sample_columns = ['Control_1', 'Control_2', 'Control_3']
        sample_metadata = {
            sample: {'Subject': 'TestControl'} for sample in sample_columns
        }
        
        with patch('matplotlib.pyplot.show'):
            cv_data = plot_control_cv_distribution(
                data=protein_data,
                sample_columns=sample_columns,
                sample_metadata=sample_metadata,
                control_column='Subject',
                control_labels=['TestControl'],
                figsize=(12, 4)
            )
        
        # Verify plot was created (figure should exist)
        assert plt.get_fignums()  # Should have at least one figure
        
        # Clean up
        plt.close('all')
        
        # Verify cleanup worked
        assert not plt.get_fignums()  # Should have no figures

    def test_plot_with_multiple_control_types(self):
        """Test plotting with multiple control types"""
        protein_data = pd.DataFrame({
            'Protein': ['P1'] * 1,
            'Description': ['Protein 1'],
            'Protein Gene': ['GENE1'],
            'UniProt_Accession': ['P00001'],
            'UniProt_Entry_Name': ['PROT1_HUMAN'],
            'ControlA_1': [100], 'ControlA_2': [110], 'ControlA_3': [95],
            'ControlB_1': [200], 'ControlB_2': [220], 'ControlB_3': [190],
            'ControlC_1': [300], 'ControlC_2': [330], 'ControlC_3': [285]
        })
        
        sample_columns = [f'Control{t}_{i}' for t in ['A', 'B', 'C'] for i in [1, 2, 3]]
        sample_metadata = {}
        for sample in sample_columns:
            control_type = sample.split('_')[0].replace('Control', 'Type')
            sample_metadata[sample] = {'Subject': control_type}
        
        with patch('matplotlib.pyplot.show'):
            cv_data = plot_control_cv_distribution(
                data=protein_data,
                sample_columns=sample_columns,
                sample_metadata=sample_metadata,
                control_column='Subject',
                control_labels=['TypeA', 'TypeB', 'TypeC'],
                figsize=(18, 6)
            )
        
        # Should return data for all three control types
        assert len(cv_data) == 3
        assert all(control_type in cv_data for control_type in ['TypeA', 'TypeB', 'TypeC'])
        
        plt.close('all')
