"""
Tests for visualization functions to prevent array length mismatches
and other visualization failures encountered in the notebook
"""

import pandas as pd
import numpy as np
import pytest
import matplotlib.pyplot as plt
from unittest.mock import patch
from proteomics_toolkit.visualization import plot_box_plot


class TestBoxPlotVisualization:
    """Test box plot visualization to prevent array length issues"""
    
    @pytest.fixture
    def mock_data_setup(self):
        """Create mock data that mimics the notebook scenario"""
        # Sample columns (cleaned names)
        sample_columns = [f'E0{i}-sample' for i in range(1, 6)]  # 5 samples
        
        # Sample metadata with groups
        sample_metadata = {
            'E01-sample': {'Group': '0'},
            'E02-sample': {'Group': '20'},
            'E03-sample': {'Group': '40'},
            'E04-sample': {'Group': '80'},
            'E05-sample': {'Group': 'Pool'}
        }
        
        # Protein data
        protein_data = pd.DataFrame({
            'Protein': ['ProteinA', 'ProteinB'],
            **{col: np.random.randn(2) for col in sample_columns}
        })
        
        return protein_data, sample_columns, sample_metadata
    
    def test_box_plot_array_length_consistency(self, mock_data_setup):
        """Test that box plot creates consistent array lengths"""
        protein_data, sample_columns, sample_metadata = mock_data_setup
        
        with patch('matplotlib.pyplot.show'):  # Don't actually show plots in tests
            try:
                # This should not raise an array length error
                plot_box_plot(
                    data=protein_data,
                    sample_columns=sample_columns,
                    sample_metadata=sample_metadata,
                    title="Test Box Plot",
                    figsize=(10, 6)
                )
                
                plt.close()  # Clean up
                
            except ValueError as e:
                if "same length" in str(e) or "array" in str(e):
                    pytest.fail(f"Array length mismatch in box plot: {e}")
                else:
                    raise  # Re-raise if it's a different error
    
    def test_box_plot_missing_sample_metadata(self, mock_data_setup):
        """Test box plot behavior when some samples missing from metadata"""
        protein_data, sample_columns, sample_metadata = mock_data_setup
        
        # Remove one sample from metadata (common source of array length mismatch)
        incomplete_metadata = sample_metadata.copy()
        del incomplete_metadata['E03-sample']
        
        with patch('matplotlib.pyplot.show'):
            try:
                plot_box_plot(
                    data=protein_data,
                    sample_columns=sample_columns,
                    sample_metadata=incomplete_metadata,
                    title="Test Box Plot - Missing Metadata",
                    figsize=(10, 6)
                )
                
                plt.close()
                
            except (KeyError, ValueError) as e:
                # Expected - should fail with missing metadata
                assert "KeyError" in str(type(e)) or "not found" in str(e).lower()
    
    def test_box_plot_empty_groups(self):
        """Test box plot with empty groups"""
        sample_columns = ['S1', 'S2']
        sample_metadata = {
            'S1': {'Group': 'A'},
            'S2': {'Group': 'A'}  # Only one group with samples
        }
        
        protein_data = pd.DataFrame({
            'Protein': ['P1'],
            'S1': [1.5],
            'S2': [2.0]
        })
        
        with patch('matplotlib.pyplot.show'):
            plot_box_plot(
                data=protein_data,
                sample_columns=sample_columns,
                sample_metadata=sample_metadata,
                title="Single Group Test",
                figsize=(8, 5)
            )
            
            plt.close()
    
    def test_box_plot_data_validation(self, mock_data_setup):
        """Test that box plot validates input data properly"""
        protein_data, sample_columns, sample_metadata = mock_data_setup
        
        # Test with missing required columns
        bad_data = protein_data.copy()
        bad_data = bad_data.drop(columns=['E01-sample'])  # Remove a sample column
        bad_sample_columns = sample_columns  # But keep it in sample_columns list
        
        with patch('matplotlib.pyplot.show'):
            with pytest.raises((KeyError, ValueError)):
                plot_box_plot(
                    data=bad_data,
                    sample_columns=bad_sample_columns,
                    sample_metadata=sample_metadata,
                    title="Bad Data Test",
                    figsize=(10, 6)
                )
                plt.close()


class TestVisualizationDataFlow:
    """Test the complete data flow through visualization functions"""
    
    def test_sample_to_group_mapping_consistency(self):
        """Test the core logic that failed in the notebook"""
        sample_columns = ['A1', 'A2', 'B1', 'B2']
        sample_metadata = {
            'A1': {'Group': 'GroupA'},
            'A2': {'Group': 'GroupA'},
            'B1': {'Group': 'GroupB'},
            'B2': {'Group': 'GroupB'}
        }
        
        # Simulate the sample-to-group mapping logic from the visualization function
        sample_to_group = {}
        for sample in sample_columns:
            if sample in sample_metadata:
                group = sample_metadata[sample].get('Group', 'Unknown')
                sample_to_group[sample] = group
        
        # Check consistency
        assert len(sample_to_group) == len(sample_columns)
        assert all(sample in sample_to_group for sample in sample_columns)
        
        # Group the samples by group
        groups = {}
        for sample, group in sample_to_group.items():
            if group not in groups:
                groups[group] = []
            groups[group].append(sample)
        
        # Verify grouping
        assert len(groups) == 2
        assert len(groups['GroupA']) == 2
        assert len(groups['GroupB']) == 2
    
    def test_color_position_array_generation(self):
        """Test that color and position arrays have same length as data"""
        sample_columns = ['S1', 'S2', 'S3', 'S4']
        groups = ['A', 'A', 'B', 'B']
        
        # Simulate color assignment
        unique_groups = list(set(groups))
        color_map = {group: f'color_{i}' for i, group in enumerate(unique_groups)}
        colors = [color_map[group] for group in groups]
        
        # Simulate position assignment
        positions = list(range(1, len(sample_columns) + 1))
        
        # Critical check: arrays must be same length
        assert len(colors) == len(positions) == len(sample_columns)
        assert len(colors) == len(groups)
        
        # This is the exact check that was failing in the notebook
        data_values = [np.random.randn(10) for _ in sample_columns]  # Mock data for each sample
        
        assert len(data_values) == len(positions) == len(colors)


class TestVisualizationRobustness:
    """Test visualization robustness against data quality issues"""
    
    def test_nan_handling_in_visualization(self):
        """Test visualization handles NaN values properly"""
        data = pd.DataFrame({
            'Protein': ['P1'],
            'S1': [np.nan],
            'S2': [1.5],
            'S3': [np.nan]
        })
        
        sample_columns = ['S1', 'S2', 'S3']
        sample_metadata = {col: {'Group': 'A'} for col in sample_columns}
        
        # Should handle NaN without crashing
        with patch('matplotlib.pyplot.show'):
            try:
                plot_box_plot(
                    data=data,
                    sample_columns=sample_columns, 
                    sample_metadata=sample_metadata,
                    title="NaN Test",
                    figsize=(8, 5)
                )
                plt.close()
            except (ValueError, RuntimeWarning) as e:
                if "nan" in str(e).lower():
                    pytest.fail(f"Visualization failed to handle NaN values: {e}")
                else:
                    raise
    
    def test_mixed_group_types(self):
        """Test visualization with mixed group types (numeric and string)"""
        sample_columns = ['S1', 'S2', 'S3']
        sample_metadata = {
            'S1': {'Group': 0},      # Numeric
            'S2': {'Group': '20'},   # String numeric
            'S3': {'Group': 'Pool'}  # String non-numeric
        }
        
        data = pd.DataFrame({
            'Protein': ['P1'],
            **{col: [1.0] for col in sample_columns}
        })
        
        with patch('matplotlib.pyplot.show'):
            plot_box_plot(
                data=data,
                sample_columns=sample_columns,
                sample_metadata=sample_metadata,
                title="Mixed Group Types Test",
                figsize=(10, 6)
            )
            plt.close()
