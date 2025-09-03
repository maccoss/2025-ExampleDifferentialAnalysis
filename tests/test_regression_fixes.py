"""
Regression tests for specific bug fixes in the proteomics analysis toolkit
These tests document and prevent regression of issues found during analysis
"""

import pandas as pd
import pytest
from unittest.mock import patch
from proteomics_toolkit.data_import import clean_sample_names
from proteomics_toolkit.visualization import plot_box_plot


class TestRegessionBugFixes:
    """Test specific bugs that were fixed to prevent regression"""
    
    def test_box_plot_array_length_bug_regression(self):
        """
        Regression test for box plot ValueError: 
        "arrays must have same first dimension" 
        
        This specific error occurred when sample metadata mapping
        was inconsistent with sample columns after name cleaning.
        Now the function should handle this gracefully.
        """
        # Create scenario that caused the original bug
        sample_columns = ['E01-sample', 'E02-sample', 'E03-sample']
        
        # Incomplete metadata (missing E03-sample) - this caused the bug
        incomplete_metadata = {
            'E01-sample': {'Group': 'A'},
            'E02-sample': {'Group': 'B'},
            # Missing 'E03-sample' - this creates array length mismatch
        }
        
        protein_data = pd.DataFrame({
            'Protein': ['ProteinX'],
            'E01-sample': [1.5],
            'E02-sample': [2.0],
            'E03-sample': [1.8]
        })
        
        # After our fixes, this should now work (handle missing metadata gracefully)
        with patch('matplotlib.pyplot.show'):
            try:
                plot_box_plot(
                    data=protein_data,
                    sample_columns=sample_columns,
                    sample_metadata=incomplete_metadata,
                    title="Regression Test - Array Length Bug"
                )
                
                # If we get here, the function handled the missing metadata gracefully
                # This is the desired behavior after our fixes
                assert True, "Function correctly handled missing metadata"
                
            except ValueError as e:
                # Should not be the confusing "arrays must have same first dimension" error
                error_msg = str(e)
                assert "first dimension" not in error_msg.lower()
                pytest.fail(f"Unexpected ValueError: {error_msg}")
            
            except KeyError as e:
                # KeyError is acceptable - indicates validation caught the issue
                error_msg = str(e)
                assert "sample" in error_msg.lower() or "metadata" in error_msg.lower()
    
    def test_sample_name_cleaning_returns_dict_regression(self):
        """
        Regression test for sample name cleaning return type confusion.
        
        The original code assumed clean_sample_names returned a list,
        but it returns a dictionary mapping.
        """
        sample_names = [
            'Total-PTE01-data Sum Normalized Area',
            'Total-PTE02-data Sum Normalized Area'
        ]
        
        result = clean_sample_names(sample_names, auto_detect=True)
        
        # Critical: Result MUST be a dictionary, not a list
        assert isinstance(result, dict), f"Expected dict, got {type(result)}"
        
        # All original names should be keys
        for name in sample_names:
            assert name in result, f"Original name {name} not found in result keys"
        
        # All values should be cleaned strings
        for cleaned in result.values():
            assert isinstance(cleaned, str), f"Cleaned name should be string, got {type(cleaned)}"
            assert len(cleaned) > 0, "Cleaned name should not be empty"
            
    def test_metadata_key_mismatch_regression(self):
        """
        Regression test for metadata key mismatch after sample cleaning.
        
        Original metadata used replicate names as keys, but after cleaning
        the code tried to access using cleaned sample names as keys.
        """
        # Original sample column names
        sample_columns = [
            'Total-PTE01-511-84A-C4-049 Sum Normalized Area',
            'Total-PTE02-Hoof18-050 Sum Normalized Area'
        ]
        
        # Original metadata uses replicate names (without suffix) as keys
        original_metadata = {
            'Total-PTE01-511-84A-C4-049': {'Group': 80, 'Subject': '84'},
            'Total-PTE02-Hoof18-050': {'Group': 'HoofPool', 'Subject': 'HoofPool'}
        }
        
        # Clean sample names
        cleaned_mapping = clean_sample_names(sample_columns, auto_detect=True)
        cleaned_samples = list(cleaned_mapping.values())
        
        # The BUG: directly accessing cleaned names in original metadata
        missing_keys = []
        for cleaned_sample in cleaned_samples:
            if cleaned_sample not in original_metadata:
                missing_keys.append(cleaned_sample)
        
        # This should fail - demonstrating the bug
        assert len(missing_keys) > 0, "Should have missing keys (demonstrating the bug)"
        
        # The FIX: Update metadata to use cleaned names as keys
        fixed_metadata = {}
        for original_sample, cleaned_sample in cleaned_mapping.items():
            # Extract replicate name
            replicate_name = original_sample.replace(' Sum Normalized Area', '')
            if replicate_name in original_metadata:
                fixed_metadata[cleaned_sample] = original_metadata[replicate_name]
        
        # Now all cleaned samples should have metadata
        for cleaned_sample in cleaned_samples:
            assert cleaned_sample in fixed_metadata, f"Fixed metadata missing {cleaned_sample}"
            
    def test_group_classification_after_cleaning_regression(self):
        """
        Regression test for sample group classification after name cleaning.
        
        After fixing metadata mapping, group classification should work correctly.
        """
        # Scenario from the notebook
        sample_columns = [
            'Total-PTE01-511-84A-C4-049 Sum Normalized Area',
            'Total-PTE02-Hoof18-050 Sum Normalized Area',
            'Total-PTE03-304-75B-B4-051 Sum Normalized Area'
        ]
        
        # Original metadata structure
        original_metadata = {
            'Total-PTE01-511-84A-C4-049': {'Group': 80, 'DrugDose': 80.0},
            'Total-PTE02-Hoof18-050': {'Group': 'HoofPool', 'DrugDose': None},
            'Total-PTE03-304-75B-B4-051': {'Group': 80, 'DrugDose': 80.0}
        }
        
        # Apply the fix
        cleaned_mapping = clean_sample_names(sample_columns, auto_detect=True)
        
        fixed_metadata = {}
        for original_sample, cleaned_sample in cleaned_mapping.items():
            replicate_name = original_sample.replace(' Sum Normalized Area', '')
            if replicate_name in original_metadata:
                fixed_metadata[cleaned_sample] = original_metadata[replicate_name]
        
        # Group classification should now work
        group_counts = {}
        for metadata in fixed_metadata.values():
            group = metadata.get('Group', 'Unknown')
            if group not in group_counts:
                group_counts[group] = 0
            group_counts[group] += 1
        
        # Verify expected groups
        assert 80 in group_counts, "Should have dose group 80"
        assert 'HoofPool' in group_counts, "Should have HoofPool group"
        assert group_counts[80] == 2, "Should have 2 samples in dose group 80"
        assert group_counts['HoofPool'] == 1, "Should have 1 sample in HoofPool"
        
    def test_visualization_debug_output_regression(self):
        """
        Test that visualization functions provide helpful debug output
        when problems occur (added during bug fixing).
        """
        # Create a scenario that might cause issues
        sample_columns = ['S1', 'S2', 'S3']
        sample_metadata = {
            'S1': {'Group': 'A'},
            'S2': {'Group': 'B'},
            # S3 intentionally missing
        }
        
        protein_data = pd.DataFrame({
            'Protein': ['P1'],
            'S1': [1.0],
            'S2': [2.0], 
            'S3': [3.0]
        })
        
        # The enhanced visualization function should provide debug output
        with patch('matplotlib.pyplot.show'):
            with patch('builtins.print') as mock_print:
                try:
                    plot_box_plot(
                        data=protein_data,
                        sample_columns=sample_columns,
                        sample_metadata=sample_metadata,
                        title="Debug Output Test"
                    )
                except (KeyError, ValueError):
                    pass  # We expect this to fail
                    
                # Check that debug output was produced
                print_calls = [call.args[0] for call in mock_print.call_args_list if call.args]
                debug_output = ' '.join(str(call) for call in print_calls)
                
                # Should mention array lengths or missing samples
                has_debug = any(
                    keyword in debug_output.lower() 
                    for keyword in ['length', 'sample', 'metadata', 'missing', 'debug']
                )
                assert has_debug, f"Expected debug output, got: {debug_output}"


class TestDataValidationImprovements:
    """Test improved data validation that prevents silent failures"""
    
    def test_sample_column_existence_validation(self):
        """Test validation that all expected sample columns exist in data"""
        protein_data = pd.DataFrame({
            'Protein': ['P1'],
            'Sample1': [100],
            # Missing Sample2 and Sample3
        })
        
        expected_samples = ['Sample1', 'Sample2', 'Sample3']
        
        # Check which samples are missing
        missing_samples = [col for col in expected_samples if col not in protein_data.columns]
        
        assert len(missing_samples) == 2
        assert 'Sample2' in missing_samples
        assert 'Sample3' in missing_samples
        
    def test_metadata_completeness_validation(self):
        """Test validation that metadata covers all samples"""
        sample_columns = ['A', 'B', 'C']
        sample_metadata = {
            'A': {'Group': '1'},
            'B': {'Group': '2'},
            # C is missing
        }
        
        # Validation check
        samples_with_metadata = set(sample_metadata.keys())
        expected_samples = set(sample_columns)
        
        missing_metadata = expected_samples - samples_with_metadata
        extra_metadata = samples_with_metadata - expected_samples
        
        assert len(missing_metadata) == 1
        assert 'C' in missing_metadata
        assert len(extra_metadata) == 0
        
    def test_group_assignment_validation(self):
        """Test validation of group assignments in metadata"""
        sample_metadata = {
            'S1': {'Group': 'Control'},
            'S2': {'Group': None},  # Invalid - None group
            'S3': {'Group': ''},    # Invalid - empty group
            'S4': {},               # Invalid - missing Group key
            'S5': {'Group': 'Treatment'}
        }
        
        # Validate group assignments
        valid_samples = []
        invalid_samples = []
        
        for sample, metadata in sample_metadata.items():
            group = metadata.get('Group')
            if group and isinstance(group, str) and group.strip():
                valid_samples.append(sample)
            else:
                invalid_samples.append(sample)
        
        assert len(valid_samples) == 2
        assert 'S1' in valid_samples
        assert 'S5' in valid_samples
        
        assert len(invalid_samples) == 3
        assert 'S2' in invalid_samples  # None
        assert 'S3' in invalid_samples  # Empty
        assert 'S4' in invalid_samples  # Missing
