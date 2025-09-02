"""
Tests for proteomics_toolkit.preprocessing module
"""
import pandas as pd
import numpy as np
from unittest.mock import patch

from proteomics_toolkit.preprocessing import (
    parse_protein_identifiers,
    parse_gene_and_description,
    assess_data_completeness,
    filter_proteins_by_completeness,
    calculate_group_colors,
    identify_annotation_columns,
    classify_samples,
    apply_systematic_color_scheme,
    _normalize_group_value
)


class TestNormalizeGroupValue:
    """Test group value normalization helper function"""
    
    def test_normalize_string_values(self):
        """Test normalizing string group values"""
        assert _normalize_group_value("Control") == "control"
        assert _normalize_group_value("TREATMENT") == "treatment"
        assert _normalize_group_value("  Mixed Case  ") == "mixed case"
        
    def test_normalize_numeric_values(self):
        """Test normalizing numeric group values"""
        assert _normalize_group_value(1) == 1
        assert _normalize_group_value(2.5) == 2.5
        assert _normalize_group_value(0) == 0
        
    def test_normalize_none_values(self):
        """Test normalizing None/NaN values"""
        result = _normalize_group_value(None)
        assert pd.isna(result)
        
        result = _normalize_group_value(np.nan)
        assert pd.isna(result)


class TestParseProteinIdentifiers:
    """Test protein identifier parsing"""
    
    def test_parse_protein_identifiers_basic(self, protein_identifiers):
        """Test basic protein identifier parsing"""
        # Create DataFrame with protein identifiers
        data = pd.DataFrame({
            'Protein': protein_identifiers,
            'Value1': [100, 200, 300, 400, 500],
            'Value2': [150, 250, 350, 450, 550]
        })
        
        result = parse_protein_identifiers(data, protein_col='Protein')
        
        assert isinstance(result, pd.DataFrame)
        assert result.shape[0] == data.shape[0]  # Same number of rows
        
        # Should have additional columns for parsed information
        expected_new_cols = ['database', 'accession', 'name', 'organism', 'entry_type']
        for col in expected_new_cols:
            assert col in result.columns
            
        # Check that valid UniProt IDs were parsed correctly
        sp_rows = result[result['database'] == 'sp']
        assert len(sp_rows) >= 3  # At least the 3 sp| entries
        
        tr_rows = result[result['database'] == 'tr']
        assert len(tr_rows) >= 1  # At least the 1 tr| entry
        
    def test_parse_protein_identifiers_custom_column(self):
        """Test parsing with custom protein column name"""
        data = pd.DataFrame({
            'ProteinID': ['sp|P12345|PROT1_HUMAN', 'tr|Q67890|Q67890_HUMAN'],
            'Value': [100, 200]
        })
        
        result = parse_protein_identifiers(data, protein_col='ProteinID')
        
        assert isinstance(result, pd.DataFrame)
        assert 'database' in result.columns
        assert 'accession' in result.columns


class TestParseGeneAndDescription:
    """Test gene and description parsing"""
    
    def test_parse_gene_and_description_basic(self):
        """Test parsing gene names and descriptions from protein descriptions"""
        data = pd.DataFrame({
            'Protein': ['P1', 'P2', 'P3'],
            'ProteinName': [
                'Protein kinase B GN=AKT1 PE=1 SV=1',
                'Cytochrome c GN=CYCS PE=1 SV=2', 
                'Simple protein description'
            ]
        })
        
        result = parse_gene_and_description(data)
        
        assert isinstance(result, pd.DataFrame)
        assert 'Gene' in result.columns
        assert 'CleanDescription' in result.columns
        
        # Check gene extraction
        assert result.iloc[0]['Gene'] == 'AKT1'
        assert result.iloc[1]['Gene'] == 'CYCS'
        assert result.iloc[2]['Gene'] == ''  # No gene info
        
        # Check description cleaning
        assert 'PE=1' not in result.iloc[0]['CleanDescription']
        assert 'SV=1' not in result.iloc[0]['CleanDescription']


class TestAssessDataCompleteness:
    """Test data completeness assessment"""
    
    def test_assess_data_completeness_basic(self, sample_protein_data, sample_columns):
        """Test basic data completeness assessment"""
        completeness = assess_data_completeness(
            sample_protein_data, 
            sample_columns,
            min_detection_rate=0.5
        )
        
        assert isinstance(completeness, dict)
        assert 'total_proteins' in completeness
        assert 'proteins_above_threshold' in completeness
        assert 'detection_rates' in completeness
        
        assert completeness['total_proteins'] == len(sample_protein_data)
        assert isinstance(completeness['detection_rates'], pd.Series)
        
    def test_assess_data_completeness_high_threshold(self, sample_protein_data, sample_columns):
        """Test completeness assessment with high threshold"""
        completeness = assess_data_completeness(
            sample_protein_data,
            sample_columns, 
            min_detection_rate=0.9  # High threshold
        )
        
        # With high threshold, fewer proteins should pass
        assert completeness['proteins_above_threshold'] <= completeness['total_proteins']


class TestFilterProteinsByCompleteness:
    """Test protein filtering by completeness"""
    
    def test_filter_proteins_basic(self, sample_protein_data, sample_columns):
        """Test basic protein filtering by completeness"""
        filtered_data = filter_proteins_by_completeness(
            sample_protein_data,
            sample_columns,
            min_detection_rate=0.5
        )
        
        assert isinstance(filtered_data, pd.DataFrame)
        assert len(filtered_data) <= len(sample_protein_data)  # Same or fewer proteins
        assert filtered_data.columns.tolist() == sample_protein_data.columns.tolist()
        
    def test_filter_proteins_strict_threshold(self, sample_protein_data, sample_columns):
        """Test filtering with strict threshold"""
        # Create data with some missing values
        data_with_missing = sample_protein_data.copy()
        # Add NaN values to first protein
        data_with_missing.loc[0, sample_columns[:3]] = np.nan
        
        filtered_data = filter_proteins_by_completeness(
            data_with_missing,
            sample_columns,
            min_detection_rate=0.9  # Strict threshold
        )
        
        # First protein should be filtered out due to missing values
        assert len(filtered_data) < len(data_with_missing)


class TestCalculateGroupColors:
    """Test group color calculation"""
    
    def test_calculate_group_colors_basic(self, sample_metadata):
        """Test basic group color calculation"""
        group_colors, group_counts = calculate_group_colors(sample_metadata)
        
        assert isinstance(group_colors, dict)
        assert isinstance(group_counts, pd.Series)
        
        # Check that all groups are assigned colors
        groups_in_metadata = set()
        for metadata in sample_metadata.values():
            if 'Group' in metadata:
                groups_in_metadata.add(metadata['Group'])
        
        for group in groups_in_metadata:
            assert group in group_colors
            assert group in group_counts.index
            
    def test_calculate_group_colors_empty_metadata(self):
        """Test group color calculation with empty metadata"""
        group_colors, group_counts = calculate_group_colors({})
        
        assert isinstance(group_colors, dict)
        assert isinstance(group_counts, pd.Series)
        assert len(group_colors) == 0
        assert len(group_counts) == 0


class TestIdentifyAnnotationColumns:
    """Test annotation column identification"""
    
    def test_identify_annotation_columns_basic(self, sample_protein_data):
        """Test basic annotation column identification"""
        annotation_cols = identify_annotation_columns(sample_protein_data)
        
        assert isinstance(annotation_cols, list)
        
        # Should identify non-numeric columns as annotations
        expected_annotation_cols = ['Protein', 'ProteinName', 'Gene']
        assert set(annotation_cols) == set(expected_annotation_cols)
        
    def test_identify_annotation_columns_numeric_only(self):
        """Test with data containing only numeric columns"""
        numeric_data = pd.DataFrame({
            'Sample_1': [100, 200, 300],
            'Sample_2': [150, 250, 350],
            'Sample_3': [120, 220, 320]
        })
        
        annotation_cols = identify_annotation_columns(numeric_data)
        
        assert isinstance(annotation_cols, list)
        assert len(annotation_cols) == 0  # No annotation columns


class TestClassifySamples:
    """Test sample classification"""
    
    def test_classify_samples_basic(self, sample_metadata):
        """Test basic sample classification"""
        result = classify_samples(
            sample_metadata,
            group_column='Group',
            group_labels=['Control', 'Treatment'],
            control_column='Group',
            control_labels=['Control']
        )
        
        assert isinstance(result, dict)
        assert 'study_samples' in result
        assert 'control_samples' in result
        assert 'group_distribution' in result
        
        # Check that samples were classified correctly
        study_samples = result['study_samples']
        control_samples = result['control_samples']
        
        assert isinstance(study_samples, list)
        assert isinstance(control_samples, list)
        
        # All samples should be classified as either study or control
        total_classified = len(study_samples) + len(control_samples)
        assert total_classified <= len(sample_metadata)
        
    def test_classify_samples_no_controls(self, sample_metadata):
        """Test classification when no control labels match"""
        result = classify_samples(
            sample_metadata,
            group_column='Group', 
            group_labels=['Control', 'Treatment'],
            control_column='Group',
            control_labels=['NonexistentControl']  # No matching controls
        )
        
        assert isinstance(result, dict)
        assert len(result['control_samples']) == 0  # No controls found


class TestApplySystematicColorScheme:
    """Test systematic color scheme application"""
    
    def test_apply_systematic_color_scheme_basic(self, sample_metadata):
        """Test applying systematic color scheme"""
        result = apply_systematic_color_scheme(
            sample_metadata,
            group_labels=['Control', 'Treatment'],
            control_labels=['Control'],
            systematic_color_palette='husl'
        )
        
        assert isinstance(result, tuple)
        assert len(result) == 2  # Should return (updated_metadata, group_colors)
        
        updated_metadata, group_colors = result
        
        assert isinstance(updated_metadata, dict)
        assert isinstance(group_colors, dict)
        
        # Check that colors were assigned
        assert 'Control' in group_colors
        assert 'Treatment' in group_colors
        
        # Colors should be valid hex colors or named colors
        for color in group_colors.values():
            assert isinstance(color, str)
            assert len(color) > 0
            
    def test_apply_systematic_color_scheme_different_palette(self, sample_metadata):
        """Test applying different color palette"""
        result = apply_systematic_color_scheme(
            sample_metadata,
            group_labels=['Control', 'Treatment'],
            control_labels=['Control'], 
            systematic_color_palette='Set1'
        )
        
        updated_metadata, group_colors = result
        
        assert isinstance(group_colors, dict)
        assert len(group_colors) >= 2  # At least Control and Treatment


class TestPreprocessingIntegration:
    """Test preprocessing module integration"""
    
    def test_full_preprocessing_pipeline(self, sample_protein_data, sample_columns, sample_metadata):
        """Test a complete preprocessing pipeline"""
        # 1. Parse protein identifiers
        data_with_ids = parse_protein_identifiers(sample_protein_data)
        
        # 2. Parse gene and description information
        data_with_genes = parse_gene_and_description(data_with_ids)
        
        # 3. Assess data completeness
        completeness = assess_data_completeness(data_with_genes, sample_columns)
        
        # 4. Filter proteins by completeness
        filtered_data = filter_proteins_by_completeness(
            data_with_genes, 
            sample_columns,
            min_detection_rate=0.3  # Lenient threshold
        )
        
        # 5. Classify samples
        sample_classification = classify_samples(
            sample_metadata,
            group_column='Group',
            group_labels=['Control', 'Treatment'],
            control_column='Group',
            control_labels=['Control']
        )
        
        # 6. Apply systematic color scheme
        updated_metadata, group_colors = apply_systematic_color_scheme(
            sample_metadata,
            group_labels=['Control', 'Treatment'],
            control_labels=['Control'],
            systematic_color_palette='husl'
        )
        
        # Verify all steps completed successfully
        assert isinstance(data_with_ids, pd.DataFrame)
        assert isinstance(data_with_genes, pd.DataFrame)
        assert isinstance(completeness, dict)
        assert isinstance(filtered_data, pd.DataFrame)
        assert isinstance(sample_classification, dict)
        assert isinstance(updated_metadata, dict)
        assert isinstance(group_colors, dict)
        
        # Check that data structure is preserved through pipeline
        assert 'Protein' in filtered_data.columns
        assert len(filtered_data) <= len(sample_protein_data)  # Same or fewer proteins


class TestPreprocessingEdgeCases:
    """Test edge cases in preprocessing"""
    
    def test_preprocessing_with_empty_data(self):
        """Test preprocessing functions with empty data"""
        empty_data = pd.DataFrame()
        
        # Most functions should handle empty data gracefully
        annotation_cols = identify_annotation_columns(empty_data)
        assert isinstance(annotation_cols, list)
        assert len(annotation_cols) == 0
        
    def test_preprocessing_with_single_sample(self):
        """Test preprocessing with only one sample"""
        single_sample_data = pd.DataFrame({
            'Protein': ['P1', 'P2'],
            'Sample_1': [100, 200]
        })
        
        completeness = assess_data_completeness(single_sample_data, ['Sample_1'])
        assert isinstance(completeness, dict)
        
        filtered_data = filter_proteins_by_completeness(
            single_sample_data, 
            ['Sample_1'],
            min_detection_rate=0.5
        )
        assert isinstance(filtered_data, pd.DataFrame)
        
    def test_preprocessing_with_all_missing_values(self):
        """Test preprocessing when sample has all missing values"""
        data_with_all_missing = pd.DataFrame({
            'Protein': ['P1', 'P2'],
            'Sample_1': [np.nan, np.nan],
            'Sample_2': [100, 200]
        })
        
        completeness = assess_data_completeness(
            data_with_all_missing, 
            ['Sample_1', 'Sample_2']
        )
        
        assert isinstance(completeness, dict)
        # Should handle the case where one sample has all missing values
