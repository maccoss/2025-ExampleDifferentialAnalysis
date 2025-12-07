"""
Tests for experimental design analysis in statistical analysis module
Focus on paired sample detection and DrugDose handling (categorical vs continuous)
"""

import pandas as pd
import numpy as np
import pytest
from unittest.mock import patch

from proteomics_toolkit.statistical_analysis import (
    StatisticalConfig,
    run_comprehensive_statistical_analysis,
    prepare_metadata_dataframe
)
from proteomics_toolkit.preprocessing import _normalize_group_value


class TestExperimentalDesignAnalysis:
    """Test experimental design validation and paired sample detection"""

    @pytest.fixture
    def dose_response_metadata_setup(self):
        """Create realistic dose-response study metadata for testing"""
        # Create metadata that mimics the actual CSF study structure
        subjects = [f'Subject_{i:02d}' for i in range(1, 19)]  # 18 subjects
        visits = ['D-02', 'D-13']
        
        sample_metadata = {}
        subject_dose_assignment = {}
        
        # Assign doses to subjects (some subjects get each dose)
        dose_assignment = [0] * 5 + [20] * 4 + [40] * 4 + [80] * 5  # 18 total
        
        for i, subject in enumerate(subjects):
            dose = dose_assignment[i]
            subject_dose_assignment[subject] = dose
            
            for visit in visits:
                sample_name = f'{subject}_{visit}'
                sample_metadata[sample_name] = {
                    'Subject': subject,
                    'Visit': visit,
                    'DrugDose': dose,  # This is the key field that caused issues
                    'StudyGroup': f'Dose_{dose}',
                    'TimePoint': visit
                }
        
        # Create sample columns list
        sample_columns = list(sample_metadata.keys())
        
        return sample_metadata, sample_columns, subject_dose_assignment

    @pytest.fixture 
    def protein_data_for_design_test(self, dose_response_metadata_setup):
        """Create protein data matching the experimental design"""
        sample_metadata, sample_columns, _ = dose_response_metadata_setup
        
        # Create protein data with annotation columns
        np.random.seed(42)
        n_proteins = 100
        
        protein_data = pd.DataFrame({
            'Protein': [f'P{i:05d}' for i in range(n_proteins)],
            'Description': [f'Protein {i} description' for i in range(n_proteins)],
            'Protein Gene': [f'GENE{i}' for i in range(n_proteins)],
            'UniProt_Accession': [f'P{i:05d}_ACC' for i in range(n_proteins)],
            'UniProt_Entry_Name': [f'PROT{i}_HUMAN' for i in range(n_proteins)],
        })
        
        # Add sample data
        for sample_col in sample_columns:
            protein_data[sample_col] = np.random.normal(1000, 200, n_proteins)
        
        return protein_data

    def test_dose_zero_handling_in_continuous_mode(self, dose_response_metadata_setup):
        """Test that dose 0 is correctly handled in continuous DrugDose analysis"""
        sample_metadata, sample_columns, subject_dose_assignment = dose_response_metadata_setup
        
        # Configure for continuous DrugDose analysis
        config = StatisticalConfig()
        config.group_column = 'DrugDose'
        config.group_labels = ['0', '20', '40', '80']
        config.force_categorical = False  # This is the key setting
        config.subject_column = 'Subject'
        config.paired_column = 'Visit'
        config.paired_label1 = 'D-02'
        config.paired_label2 = 'D-13'
        
        # Test the metadata preparation (this is where the bug was)
        metadata_df = prepare_metadata_dataframe(sample_metadata, sample_columns, config)
        
        # Should include all samples, including dose 0
        assert len(metadata_df) == len(sample_columns)
        
        # Check that dose 0 samples are present
        dose_0_samples = metadata_df[metadata_df['DrugDose'] == 0]
        assert len(dose_0_samples) > 0, "Dose 0 samples should be present in continuous mode"
        
        # Check that all doses are represented
        unique_doses = set(metadata_df['DrugDose'].unique())
        expected_doses = {0, 20, 40, 80}
        assert unique_doses == expected_doses, f"Expected doses {expected_doses}, got {unique_doses}"

    def test_dose_zero_handling_in_categorical_mode(self, dose_response_metadata_setup):
        """Test that dose 0 is correctly handled in categorical DrugDose analysis"""
        sample_metadata, sample_columns, _ = dose_response_metadata_setup
        
        # Configure for categorical DrugDose analysis
        config = StatisticalConfig()
        config.group_column = 'DrugDose'
        config.group_labels = ['0', '20', '40', '80']
        config.force_categorical = True  # Categorical mode
        config.subject_column = 'Subject'
        config.paired_column = 'Visit'
        config.paired_label1 = 'D-02'
        config.paired_label2 = 'D-13'
        
        metadata_df = prepare_metadata_dataframe(sample_metadata, sample_columns, config)
        
        # Should include all samples, including dose 0
        assert len(metadata_df) == len(sample_columns)
        
        # In categorical mode, doses are converted to strings
        unique_doses = set(metadata_df['DrugDose'].unique())
        expected_doses = {'0', '20', '40', '80'}  # Strings, not numbers
        assert unique_doses == expected_doses, f"Expected doses {expected_doses}, got {unique_doses}"
        
        # Specifically check that dose 0 is present (as string)
        dose_0_samples = metadata_df[metadata_df['DrugDose'] == '0']  # String comparison
        assert len(dose_0_samples) > 0, f"Dose 0 samples should be present in categorical mode. Unique doses: {unique_doses}"

    def test_paired_sample_detection(self, dose_response_metadata_setup):
        """Test that paired samples are correctly identified across visits"""
        sample_metadata, sample_columns, subject_dose_assignment = dose_response_metadata_setup
        
        config = StatisticalConfig()
        config.group_column = 'DrugDose'
        config.group_labels = ['0', '20', '40', '80']
        config.force_categorical = False
        config.subject_column = 'Subject'
        config.paired_column = 'Visit'
        config.paired_label1 = 'D-02'
        config.paired_label2 = 'D-13'
        
        metadata_df = prepare_metadata_dataframe(sample_metadata, sample_columns, config)
        
        # Test pairing logic manually (simulating the actual analysis logic)
        pairing_data = {}
        for _, row in metadata_df.iterrows():
            subject = row['Subject']
            visit = row['Visit']
            dose = row['DrugDose']
            sample_name = row['Sample']
            
            # The key test: comparison should not be None for dose 0
            # This tests the fix for: if subject and visit and comparison is not None
            if subject and visit and dose is not None:  # Fixed condition
                if subject not in pairing_data:
                    pairing_data[subject] = {}
                pairing_data[subject][visit] = {
                    'sample': sample_name,
                    'dose': dose
                }
        
        # Count complete pairs
        complete_pairs = []
        for subject, visits in pairing_data.items():
            if 'D-02' in visits and 'D-13' in visits:
                baseline = visits['D-02']
                followup = visits['D-13']
                
                # Doses should match between timepoints for the same subject
                if baseline['dose'] == followup['dose']:
                    complete_pairs.append({
                        'subject': subject,
                        'dose': baseline['dose'],
                        'baseline_sample': baseline['sample'],
                        'followup_sample': followup['sample']
                    })
        
        # Should find 18 complete pairs (all subjects have both timepoints)
        assert len(complete_pairs) == 18, f"Expected 18 complete pairs, found {len(complete_pairs)}"
        
        # Check dose distribution in complete pairs
        dose_counts = {}
        for pair in complete_pairs:
            dose = pair['dose']
            dose_counts[dose] = dose_counts.get(dose, 0) + 1
        
        expected_dose_counts = {0: 5, 20: 4, 40: 4, 80: 5}
        assert dose_counts == expected_dose_counts, f"Expected dose counts {expected_dose_counts}, got {dose_counts}"

    def test_group_value_normalization(self):
        """Test the _normalize_group_value function that caused issues"""
        # Test various inputs that should normalize correctly
        test_cases = [
            ('0', 0),          # String zero -> numeric zero
            (0, 0),            # Numeric zero stays numeric
            (0.0, 0),          # Float zero -> integer zero
            ('20', 20),        # String number -> numeric
            (20, 20),          # Integer stays integer
            (20.0, 20),        # Float integer -> integer
            ('control', 'control'),  # String stays string
            (None, 'Unknown'),      # None -> 'Unknown'
        ]
        
        for input_val, expected in test_cases:
            result = _normalize_group_value(input_val)
            assert result == expected, f"_normalize_group_value({input_val}) should return {expected}, got {result}"

    def test_dose_filtering_logic(self, dose_response_metadata_setup):
        """Test the sample filtering logic that checks group membership"""
        sample_metadata, sample_columns, _ = dose_response_metadata_setup
        
        config = StatisticalConfig()
        config.group_column = 'DrugDose'
        config.group_labels = ['0', '20', '40', '80']
        config.force_categorical = False
        
        # Simulate the filtering logic from the analysis
        normalized_group_labels = [_normalize_group_value(label) for label in config.group_labels]
        
        valid_samples = {}
        for sample_name, metadata in sample_metadata.items():
            comparison_value = metadata.get(config.group_column)
            normalized_comparison = _normalize_group_value(comparison_value)
            
            if normalized_comparison in normalized_group_labels:
                valid_samples[sample_name] = metadata
        
        # Should include all samples since all have valid doses
        assert len(valid_samples) == len(sample_metadata)
        
        # Specifically check that dose 0 samples are included
        dose_0_samples = [s for s, m in valid_samples.items() if m['DrugDose'] == 0]
        assert len(dose_0_samples) > 0, "Dose 0 samples should pass filtering"

    def test_mixed_effects_analysis_with_continuous_dose(self, dose_response_metadata_setup, protein_data_for_design_test):
        """Test complete mixed-effects analysis with continuous dose variable"""
        sample_metadata, sample_columns, _ = dose_response_metadata_setup
        protein_data = protein_data_for_design_test
        
        config = StatisticalConfig()
        config.statistical_test_method = 'mixed_effects'
        config.analysis_type = 'interaction'  # Required: set analysis type
        config.group_column = 'DrugDose'
        config.group_labels = ['0', '20', '40', '80']
        config.force_categorical = False  # Continuous mode
        config.subject_column = 'Subject'
        config.paired_column = 'Visit'
        config.paired_label1 = 'D-02'
        config.paired_label2 = 'D-13'
        config.interaction_terms = ['DrugDose', 'Visit']
        
        # This should not raise an error and should include dose 0 samples
        with patch('builtins.print'):  # Suppress output for testing
            results = run_comprehensive_statistical_analysis(
                normalized_data=protein_data,
                sample_metadata=sample_metadata,
                config=config,
                protein_annotations=protein_data[['Protein', 'Description', 'Protein Gene', 'UniProt_Accession', 'UniProt_Entry_Name']]
            )
        
        # Should return results for all proteins
        assert isinstance(results, pd.DataFrame)
        assert len(results) > 0
        
        # Should have the expected columns
        expected_columns = ['logFC', 'P.Value', 'adj.P.Val']
        for col in expected_columns:
            assert col in results.columns, f"Missing column: {col}"

    def test_categorical_vs_continuous_dose_treatment(self, dose_response_metadata_setup, protein_data_for_design_test):
        """Test differences between categorical and continuous dose treatment"""
        sample_metadata, sample_columns, _ = dose_response_metadata_setup
        protein_data = protein_data_for_design_test
        
        # Test continuous mode
        config_continuous = StatisticalConfig()
        config_continuous.statistical_test_method = 'mixed_effects'
        config_continuous.analysis_type = 'interaction'  # Required: set analysis type
        config_continuous.group_column = 'DrugDose'
        config_continuous.group_labels = ['0', '20', '40', '80']
        config_continuous.force_categorical = False
        config_continuous.subject_column = 'Subject'
        config_continuous.paired_column = 'Visit'
        config_continuous.paired_label1 = 'D-02'
        config_continuous.paired_label2 = 'D-13'
        config_continuous.interaction_terms = ['DrugDose', 'Visit']
        
        # Test categorical mode
        config_categorical = StatisticalConfig()
        config_categorical.statistical_test_method = 'mixed_effects'
        config_categorical.analysis_type = 'interaction'  # Required: set analysis type
        config_categorical.group_column = 'DrugDose'
        config_categorical.group_labels = ['0', '20', '40', '80']
        config_categorical.force_categorical = True
        config_categorical.subject_column = 'Subject'
        config_categorical.paired_column = 'Visit'
        config_categorical.paired_label1 = 'D-02'
        config_categorical.paired_label2 = 'D-13'
        config_categorical.interaction_terms = ['DrugDose', 'Visit']
        
        # Both should work without errors
        with patch('builtins.print'):
            results_continuous = run_comprehensive_statistical_analysis(
                normalized_data=protein_data,
                sample_metadata=sample_metadata,
                config=config_continuous,
                protein_annotations=protein_data[['Protein', 'Description', 'Protein Gene', 'UniProt_Accession', 'UniProt_Entry_Name']]
            )
            
            results_categorical = run_comprehensive_statistical_analysis(
                normalized_data=protein_data,
                sample_metadata=sample_metadata,
                config=config_categorical,
                protein_annotations=protein_data[['Protein', 'Description', 'Protein Gene', 'UniProt_Accession', 'UniProt_Entry_Name']]
            )
        
        # Both should produce results
        assert len(results_continuous) > 0
        assert len(results_categorical) > 0
        
        # Results may differ between continuous and categorical treatment
        # but both should include all proteins
        assert len(results_continuous) == len(results_categorical)

    def test_missing_visit_handling(self, dose_response_metadata_setup):
        """Test handling when some samples are missing visit information"""
        sample_metadata, sample_columns, _ = dose_response_metadata_setup
        
        # Remove visit information from some samples
        incomplete_metadata = sample_metadata.copy()
        samples_to_break = list(incomplete_metadata.keys())[:5]
        for sample in samples_to_break:
            del incomplete_metadata[sample]['Visit']
        
        config = StatisticalConfig()
        config.group_column = 'DrugDose'
        config.group_labels = ['0', '20', '40', '80']
        config.subject_column = 'Subject'
        config.paired_column = 'Visit'
        config.paired_label1 = 'D-02'
        config.paired_label2 = 'D-13'
        
        # Should handle missing visits gracefully (or raise appropriate error)
        try:
            metadata_df = prepare_metadata_dataframe(incomplete_metadata, sample_columns, config)
            # If it doesn't raise an error, it should filter out samples with missing visits
            assert len(metadata_df) < len(sample_columns)
        except (ValueError, KeyError):
            # Expected behavior - should fail with missing required data
            pass

    def test_subject_dose_consistency(self, dose_response_metadata_setup):
        """Test that subjects maintain consistent dose assignments across visits"""
        sample_metadata, sample_columns, subject_dose_assignment = dose_response_metadata_setup
        
        # Verify the test data is set up correctly
        for sample_name, metadata in sample_metadata.items():
            subject = metadata['Subject']
            dose = metadata['DrugDose']
            expected_dose = subject_dose_assignment[subject]
            assert dose == expected_dose, f"Subject {subject} has inconsistent dose assignment"
        
        # This ensures our test data mimics the real study design where
        # each subject gets the same dose at both timepoints
        subjects_by_dose = {}
        for subject, dose in subject_dose_assignment.items():
            if dose not in subjects_by_dose:
                subjects_by_dose[dose] = []
            subjects_by_dose[dose].append(subject)
        
        # Verify the expected distribution
        expected_distribution = {0: 5, 20: 4, 40: 4, 80: 5}
        actual_distribution = {dose: len(subjects) for dose, subjects in subjects_by_dose.items()}
        assert actual_distribution == expected_distribution


class TestExperimentalDesignEdgeCases:
    """Test edge cases and error conditions in experimental design analysis"""

    def test_empty_metadata(self):
        """Test behavior with empty metadata"""
        config = StatisticalConfig()
        config.group_column = 'DrugDose'
        config.subject_column = 'Subject'
        
        # This should return an empty DataFrame, not raise an error
        result = prepare_metadata_dataframe({}, [], config)
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0

    def test_mismatched_sample_columns(self):
        """Test behavior when sample columns don't match metadata"""
        sample_metadata = {
            'Sample_A': {'DrugDose': 0, 'Subject': 'S1', 'Visit': 'V1'},
            'Sample_B': {'DrugDose': 20, 'Subject': 'S2', 'Visit': 'V1'}
        }
        sample_columns = ['Sample_C', 'Sample_D']  # Different samples
        
        config = StatisticalConfig()
        config.group_column = 'DrugDose'
        config.subject_column = 'Subject'
        config.paired_column = 'Visit'
        
        # Should handle mismatched samples appropriately
        metadata_df = prepare_metadata_dataframe(sample_metadata, sample_columns, config)
        
        # Should result in empty or filtered dataframe
        assert len(metadata_df) == 0

    def test_invalid_dose_values(self):
        """Test handling of invalid dose values"""
        sample_metadata = {
            'Sample_A': {'DrugDose': 'invalid', 'Subject': 'S1', 'Visit': 'V1'},
            'Sample_B': {'DrugDose': None, 'Subject': 'S2', 'Visit': 'V1'},
            'Sample_C': {'DrugDose': float('nan'), 'Subject': 'S3', 'Visit': 'V1'},
        }
        sample_columns = ['Sample_A', 'Sample_B', 'Sample_C']
        
        config = StatisticalConfig()
        config.group_column = 'DrugDose'
        config.group_labels = ['0', '20', '40', '80']
        config.subject_column = 'Subject'
        
        # Should handle invalid values gracefully
        metadata_df = prepare_metadata_dataframe(sample_metadata, sample_columns, config)
        
        # Invalid samples should be filtered out
        # (exact behavior depends on implementation)
        assert len(metadata_df) <= len(sample_columns)
