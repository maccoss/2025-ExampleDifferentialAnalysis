"""
Tests for sample name cleaning and metadata mapping integration
These tests address the critical issues encountered in the notebook where
sample name cleaning broke the metadata mapping.
"""

import pandas as pd
import pytest
from proteomics_toolkit.data_import import clean_sample_names


class TestSampleNameMetadataIntegration:
    """Test integration between sample name cleaning and metadata mapping"""
    
    @pytest.fixture
    def sample_data_with_metadata(self):
        """Create test data that mimics the real notebook scenario"""
        # Original sample names with common prefix (like the notebook)
        original_samples = [
            "Total-PTE01-511-84A-C4-049 Sum Normalized Area",
            "Total-PTE02-Hoof18-050 Sum Normalized Area", 
            "Total-PTE03-304-75B-B4-051 Sum Normalized Area"
        ]
        
        # Metadata with original sample names as keys
        metadata_df = pd.DataFrame({
            'Replicate': ["Total-PTE01-511-84A-C4-049", "Total-PTE02-Hoof18-050", "Total-PTE03-304-75B-B4-051"],
            'Group': [80, 'HoofPool', 80],
            'DrugDose': [80.0, None, 80.0],
            'Subject': ['84', 'HoofPool', '75'],
            'Visit': ['D-02', None, 'D-13']
        })
        
        # Original sample metadata dictionary (from metadata.set_index(...).to_dict('index'))
        sample_metadata = metadata_df.set_index('Replicate').to_dict('index')
        
        return original_samples, sample_metadata, metadata_df
    
    def test_sample_name_cleaning_preserves_metadata_mapping(self, sample_data_with_metadata):
        """Test that sample name cleaning properly updates metadata mapping"""
        original_samples, original_metadata, metadata_df = sample_data_with_metadata
        
        # Clean sample names (this returns a dict mapping)
        cleaned_mapping = clean_sample_names(
            original_samples,
            auto_detect=True
        )
        
        # Extract cleaned names list
        cleaned_samples = list(cleaned_mapping.values())
        
        # The critical step: Update metadata to use cleaned names as keys
        updated_metadata = {}
        for original_name, cleaned_name in cleaned_mapping.items():
            # Extract the replicate name from the original sample column name
            replicate_name = original_name.replace(' Sum Normalized Area', '')
            if replicate_name in original_metadata:
                updated_metadata[cleaned_name] = original_metadata[replicate_name]
        
        # Verify the mapping works
        assert len(updated_metadata) == len(cleaned_samples)
        assert all(sample in updated_metadata for sample in cleaned_samples)
        
        # Verify metadata content is preserved
        for cleaned_name, metadata in updated_metadata.items():
            assert 'Group' in metadata
            assert 'Subject' in metadata
            
    def test_sample_classification_after_name_cleaning(self, sample_data_with_metadata):
        """Test that sample classification works after name cleaning"""
        original_samples, original_metadata, metadata_df = sample_data_with_metadata
        
        # Clean names and update metadata (simulating notebook fix)
        cleaned_mapping = clean_sample_names(original_samples, auto_detect=True)
        
        updated_metadata = {}
        for original_name, cleaned_name in cleaned_mapping.items():
            replicate_name = original_name.replace(' Sum Normalized Area', '')
            if replicate_name in original_metadata:
                updated_metadata[cleaned_name] = original_metadata[replicate_name]
        
        # This should work without errors now
        try:
            # Mock the classify_samples call (simplified)
            group_distribution = {}
            for sample, metadata in updated_metadata.items():
                group = metadata.get('Group', 'Unknown')
                if pd.isna(group):
                    group = metadata.get('Subject', 'Unknown')
                
                if group not in group_distribution:
                    group_distribution[group] = 0
                group_distribution[group] += 1
                    
            assert len(group_distribution) > 0
            assert 80 in group_distribution or '80' in group_distribution
            assert 'HoofPool' in group_distribution
            
        except Exception as e:
            pytest.fail(f"Sample classification failed after name cleaning: {e}")

    def test_metadata_mapping_consistency_check(self):
        """Test function to verify metadata mapping consistency"""
        # This is a utility test function that should be part of the toolkit
        sample_columns = ['E01-sample', 'E02-sample', 'E03-sample']
        metadata = {
            'E01-sample': {'Group': 'A'}, 
            'E02-sample': {'Group': 'B'}, 
            'different-name': {'Group': 'C'}  # This one doesn't match
        }
        
        matching_samples = [s for s in sample_columns if s in metadata]
        missing_samples = [s for s in sample_columns if s not in metadata]
        
        assert len(matching_samples) == 2
        assert len(missing_samples) == 1
        assert missing_samples[0] == 'E03-sample'


class TestSampleNameCleaningEdgeCases:
    """Test edge cases in sample name cleaning that caused issues"""
    
    def test_empty_cleaning_result(self):
        """Test when cleaning produces no changes"""
        samples = ['SampleA', 'SampleB', 'SampleC']
        result = clean_sample_names(samples, auto_detect=False)  # Don't auto-detect
        
        # Should return identity mapping when no cleaning parameters provided
        expected = {s: s for s in samples}
        assert result == expected
    
    def test_partial_prefix_cleaning(self):
        """Test when only some samples have the prefix"""
        samples = ['Total-PT01-data', 'Total-PT02-data', 'Different03-data']
        result = clean_sample_names(samples, auto_detect=True)
        
        # Auto-detect should find the common pattern
        assert len(result) == 3
        assert all(isinstance(k, str) and isinstance(v, str) for k, v in result.items())
    
    def test_complex_prefix_suffix_combination(self):
        """Test the complex real-world naming pattern from the notebook"""
        samples = [
            'Total-PTE01-511-84A-C4-049 Sum Normalized Area',
            'Total-PTE02-Hoof18-050 Sum Normalized Area'  
        ]
        
        result = clean_sample_names(samples, auto_detect=True)
        
        # Should remove common prefix and suffix
        for original, cleaned in result.items():
            assert not cleaned.startswith('Total-PT')
            assert not cleaned.endswith(' Sum Normalized Area')
            assert len(cleaned) > 0
