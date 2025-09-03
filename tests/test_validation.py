"""
Tests for the validation module - metadata/data consistency checking.

Tests various scenarios where metadata samples match or don't match with 
protein data columns, including control sample validation.
"""

import pytest
import pandas as pd
from unittest.mock import MagicMock
import sys
import os

# Add the parent directory to the path so we can import the module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import proteomics_toolkit as ptk


class TestMetadataDataConsistency:
    """Test validation of metadata/data file consistency."""
    
    def setup_method(self):
        """Set up test data for each test."""
        # Create mock metadata
        self.metadata = pd.DataFrame({
            'Sample_Name': ['Total-PT-001', 'Total-PT-002', 'Total-PT-PlatePool1', 'Total-PT-HoofPool1', 'Total-PT-003'],
            'DrugDose': [0, 20, 'PlatePool', 'HoofPool', 40],
            'Subject': ['Patient1', 'Patient2', 'PlatePool', 'HoofPool', 'Patient3'],
            'Visit': ['D-02', 'D-02', 'Pool', 'Pool', 'D-13']
        })
        
        # Mock protein data columns - some samples missing
        self.protein_columns = [
            'Protein', 'UniProt', 'Gene', 'Description',
            'Total-PT-001 Sum Normalized Area', 
            'Total-PT-002 Sum Normalized Area',
            'Total-PT-PlatePool1 Sum Normalized Area',
            # 'Total-PT-HoofPool1 Sum Normalized Area', # Missing!
            'Total-PT-003 Sum Normalized Area'
        ]
        
        self.metadata_sample_names = self.metadata.iloc[:, 0].tolist()
        self.control_column = 'Subject'
        self.control_labels = ['PlatePool', 'HoofPool', 'GWPool']

    def test_validate_all_samples_found(self):
        """Test validation when all samples are found in data."""
        # All samples present
        complete_protein_columns = [
            'Protein', 'UniProt', 'Gene', 'Description',
            'Total-PT-001 Sum Normalized Area', 
            'Total-PT-002 Sum Normalized Area',
            'Total-PT-PlatePool1 Sum Normalized Area',
            'Total-PT-HoofPool1 Sum Normalized Area',
            'Total-PT-003 Sum Normalized Area'
        ]
        
        results = ptk.validate_metadata_data_consistency(
            metadata=self.metadata,
            metadata_sample_names=self.metadata_sample_names,
            protein_columns=complete_protein_columns,
            control_column=self.control_column,
            control_labels=self.control_labels,
            verbose=False
        )
        
        assert results['is_valid'] is True
        assert len(results['errors']) == 0
        assert results['diagnostics']['samples_missing_from_data'] == 0
        assert results['diagnostics']['control_samples_missing_from_data'] == 0
        assert results['diagnostics']['control_samples_found_in_data'] == 2  # PlatePool, HoofPool

    def test_validate_missing_samples(self):
        """Test validation when some samples are missing from data."""
        results = ptk.validate_metadata_data_consistency(
            metadata=self.metadata,
            metadata_sample_names=self.metadata_sample_names,
            protein_columns=self.protein_columns,
            control_column=self.control_column,
            control_labels=self.control_labels,
            verbose=False
        )
        
        assert results['is_valid'] is False
        assert len(results['errors']) > 0
        assert results['diagnostics']['samples_missing_from_data'] == 1
        assert 'Total-PT-HoofPool1' in results['diagnostics']['missing_samples']
        assert results['diagnostics']['control_samples_missing_from_data'] == 1

    def test_validate_missing_control_column(self):
        """Test validation when control column doesn't exist in metadata."""
        results = ptk.validate_metadata_data_consistency(
            metadata=self.metadata,
            metadata_sample_names=self.metadata_sample_names,
            protein_columns=self.protein_columns,
            control_column='NonexistentColumn',
            control_labels=self.control_labels,
            verbose=False
        )
        
        assert results['is_valid'] is False
        assert any("Control column 'NonexistentColumn' not found" in error 
                  for error in results['errors'])

    def test_validate_no_control_samples(self):
        """Test validation when no control samples are defined."""
        results = ptk.validate_metadata_data_consistency(
            metadata=self.metadata,
            metadata_sample_names=self.metadata_sample_names,
            protein_columns=self.protein_columns,
            control_column=self.control_column,
            control_labels=['NonexistentControl'],  # No matching controls
            verbose=False
        )
        
        assert results['diagnostics']['total_control_samples_in_metadata'] == 0
        assert results['diagnostics']['control_samples_missing_from_data'] == 0

    def test_enhanced_sample_processing_strict_mode(self):
        """Test enhanced sample processing with strict validation (should raise error)."""
        protein_data = pd.DataFrame({
            'Protein': ['P1', 'P2'],
            'Total-PT-001 Sum Normalized Area': [100, 200],
            'Total-PT-002 Sum Normalized Area': [150, 250],
            'Total-PT-PlatePool1 Sum Normalized Area': [120, 220],
            # Missing HoofPool sample
            'Total-PT-003 Sum Normalized Area': [180, 280]
        })
        
        mock_toolkit = MagicMock()
        mock_toolkit.data_import.clean_sample_names.return_value = {
            'Total-PT-001 Sum Normalized Area': 'PT-001',
            'Total-PT-002 Sum Normalized Area': 'PT-002', 
            'Total-PT-PlatePool1 Sum Normalized Area': 'PT-PlatePool1',
            'Total-PT-003 Sum Normalized Area': 'PT-003'
        }
        mock_toolkit.classify_samples.return_value = ({}, [], [], {}, {})
        
        with pytest.raises(ptk.SampleMatchingError):
            ptk.enhanced_sample_processing(
                metadata=self.metadata,
                protein_data=protein_data,
                group_column='DrugDose',
                group_labels=['0', '20', '40'],
                control_column=self.control_column,
                control_labels=self.control_labels,
                toolkit_module=mock_toolkit,
                strict_validation=True
            )

    def test_enhanced_sample_processing_lenient_mode(self):
        """Test enhanced sample processing with lenient validation (should warn but continue)."""
        protein_data = pd.DataFrame({
            'Protein': ['P1', 'P2'],
            'Total-PT-001 Sum Normalized Area': [100, 200],
            'Total-PT-002 Sum Normalized Area': [150, 250],
            'Total-PT-PlatePool1 Sum Normalized Area': [120, 220],
            'Total-PT-003 Sum Normalized Area': [180, 280]
        })
        
        mock_toolkit = MagicMock()
        mock_toolkit.data_import.clean_sample_names.return_value = {
            'Total-PT-001 Sum Normalized Area': 'PT-001',
            'Total-PT-002 Sum Normalized Area': 'PT-002',
            'Total-PT-PlatePool1 Sum Normalized Area': 'PT-PlatePool1', 
            'Total-PT-003 Sum Normalized Area': 'PT-003'
        }
        mock_toolkit.classify_samples.return_value = ({'0': 1, '20': 1, '40': 1, 'PlatePool': 1}, 
                                                     ['PT-PlatePool1'], 
                                                     ['PT-001', 'PT-002', 'PT-003'], 
                                                     {}, 
                                                     {})
        
        # Should not raise an exception in lenient mode
        try:
            result = ptk.enhanced_sample_processing(
                metadata=self.metadata,
                protein_data=protein_data,
                group_column='DrugDose',
                group_labels=['0', '20', '40'],
                control_column=self.control_column,
                control_labels=self.control_labels,
                toolkit_module=mock_toolkit,
                strict_validation=False
            )
            assert result is not None
        except ptk.SampleMatchingError:
            pytest.fail("Should not raise SampleMatchingError in lenient mode")

    def test_diagnostic_report_generation(self):
        """Test generation of diagnostic reports."""
        results = ptk.validate_metadata_data_consistency(
            metadata=self.metadata,
            metadata_sample_names=self.metadata_sample_names,
            protein_columns=self.protein_columns,
            control_column=self.control_column,
            control_labels=self.control_labels,
            verbose=False
        )
        
        report = ptk.generate_sample_matching_diagnostic_report(results)
        
        assert "SAMPLE MATCHING DIAGNOSTIC REPORT" in report
        assert "SUMMARY STATISTICS" in report
        assert "CONTROL SAMPLE ANALYSIS" in report
        assert "MISSING SAMPLES DETAILS" in report
        assert "RECOMMENDATIONS" in report
        assert "Total-PT-HoofPool1" in report  # Missing sample should be listed


class TestEdgeCases:
    """Test edge cases and error conditions."""
    
    def test_empty_metadata(self):
        """Test validation with empty metadata."""
        empty_metadata = pd.DataFrame()
        
        results = ptk.validate_metadata_data_consistency(
            metadata=empty_metadata,
            metadata_sample_names=[],
            protein_columns=['Protein', 'Sample1'],
            control_column='Subject',
            control_labels=['Control'],
            verbose=False
        )
        
        assert results['diagnostics']['total_metadata_samples'] == 0
        
    def test_empty_protein_columns(self):
        """Test validation with no protein data columns."""
        metadata = pd.DataFrame({
            'Sample_Name': ['Sample1', 'Sample2'],
            'Subject': ['Patient1', 'Control']
        })
        
        results = ptk.validate_metadata_data_consistency(
            metadata=metadata,
            metadata_sample_names=['Sample1', 'Sample2'],
            protein_columns=['Protein', 'Gene'],  # No sample columns
            control_column='Subject',
            control_labels=['Control'],
            verbose=False
        )
        
        assert results['is_valid'] is False
        assert results['diagnostics']['samples_missing_from_data'] == 2

    def test_partial_sample_name_matching(self):
        """Test validation with partial sample name matching."""
        metadata = pd.DataFrame({
            'Sample_Name': ['Sample-001', 'Sample-002-Control'],
            'Type': ['Study', 'Control']
        })
        
        # Protein columns have partial matches
        protein_columns = [
            'Protein', 'Gene', 
            'Sample-001-Intensity',  # Partial match for Sample-001
            'Sample-002-Control-Value'  # Partial match for Sample-002-Control 
        ]
        
        results = ptk.validate_metadata_data_consistency(
            metadata=metadata,
            metadata_sample_names=['Sample-001', 'Sample-002-Control'],
            protein_columns=protein_columns,
            control_column='Type',
            control_labels=['Control'],
            verbose=False
        )
        
        # Both samples should be found due to substring matching
        assert results['diagnostics']['samples_found_in_data'] == 2
        assert results['diagnostics']['samples_missing_from_data'] == 0


class TestRealWorldScenarios:
    """Test scenarios based on real-world data issues."""
    
    def test_eisai_to_plate_pool_scenario(self):
        """Test the specific EISAIPool -> PlatePool labeling issue we encountered."""
        metadata = pd.DataFrame({
            'Sample_Name': ['Total-PT-001', 'Total-PT-PlatePool1'],
            'Subject': ['Patient1', 'PlatePool'],
            'DrugDose': [0, 'PlatePool']
        })
        
        # Protein data still has old EISAIPool naming
        protein_columns = [
            'Protein', 'Gene',
            'Total-PT-001 Sum Normalized Area',
            'Total-PT-EISAIPool1 Sum Normalized Area'  # Mismatched name!
        ]
        
        results = ptk.validate_metadata_data_consistency(
            metadata=metadata,
            metadata_sample_names=['Total-PT-001', 'Total-PT-PlatePool1'],
            protein_columns=protein_columns,
            control_column='Subject',
            control_labels=['PlatePool'],
            verbose=False
        )
        
        assert results['is_valid'] is False
        assert 'Total-PT-PlatePool1' in results['diagnostics']['missing_samples']
        assert results['diagnostics']['control_samples_missing_from_data'] == 1

    def test_multiple_control_types_missing(self):
        """Test scenario with multiple control types, some missing."""
        metadata = pd.DataFrame({
            'Sample_Name': ['Total-PT-001', 'Total-PT-HoofPool1', 'Total-PT-GWPool1', 'Total-PT-PlatePool1'],
            'Subject': ['Patient1', 'HoofPool', 'GWPool', 'PlatePool']
        })
        
        protein_columns = [
            'Protein', 'Gene',
            'Total-PT-001 Sum Normalized Area',
            'Total-PT-HoofPool1 Sum Normalized Area',
            # Missing GWPool and PlatePool
        ]
        
        results = ptk.validate_metadata_data_consistency(
            metadata=metadata,
            metadata_sample_names=metadata.iloc[:, 0].tolist(),
            protein_columns=protein_columns,
            control_column='Subject',
            control_labels=['HoofPool', 'GWPool', 'PlatePool'],
            verbose=False
        )
        
        assert results['is_valid'] is False
        assert results['diagnostics']['control_samples_missing_from_data'] == 2
        
        # Check the diagnostic report mentions the missing control types
        report = ptk.generate_sample_matching_diagnostic_report(results)
        assert 'GWPool' in report
        assert 'PlatePool' in report

    def test_sample_name_prefix_variations(self):
        """Test different sample name prefix patterns."""
        metadata = pd.DataFrame({
            'Sample_Name': ['PT-001', 'CSF-Total-002', 'Sample_003'],
            'Group': ['A', 'B', 'C']
        })
        
        protein_columns = [
            'Protein', 
            'PT-001 Intensity',  # Direct match
            'CSF-Total-002 Area',  # Direct match
            'Sample_003_Normalized'  # Underscore variation
        ]
        
        results = ptk.validate_metadata_data_consistency(
            metadata=metadata,
            metadata_sample_names=metadata.iloc[:, 0].tolist(),
            protein_columns=protein_columns,
            control_column='Group',
            control_labels=[],  # No controls in this test
            verbose=False
        )
        
        # All should be found due to flexible matching
        assert results['diagnostics']['samples_found_in_data'] == 3
        assert results['is_valid'] is True
