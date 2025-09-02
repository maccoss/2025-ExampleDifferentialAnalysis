"""
Tests for create_standard_data_structure function in preprocessing module.
"""

import pytest
import pandas as pd
import numpy as np
from proteomics_toolkit.preprocessing import create_standard_data_structure, parse_protein_identifiers, parse_gene_and_description


class TestCreateStandardDataStructure:
    """Test the create_standard_data_structure function."""
    
    @pytest.fixture
    def sample_protein_data(self):
        """Create sample protein data that mimics the CSV input format."""
        return pd.DataFrame({
            'Protein': [
                'sp|P02787|TRFE_HUMAN',
                'sp|O43240|KLK10_HUMAN',
                'sp|P06858|LIPL_HUMAN'
            ],
            'Protein Description': [
                'Serotransferrin OS=Homo sapiens OX=9606 GN=TF PE=1 SV=5',
                'Kallikrein-10 OS=Homo sapiens OX=9606 GN=KLK10 PE=1 SV=3', 
                'Lipoprotein lipase OS=Homo sapiens OX=9606 GN=LPL PE=1 SV=1'
            ],
            'Protein Gene': ['TF', 'KLK10', 'LPL'],
            'Total-PTE01-511-84A-C4-049 Sum Normalized Area': [1.2e8, 3.7e3, 1.7e3],
            'Total-PTE02-Hoof18-050 Sum Normalized Area': [1.9e8, 6.0e2, 3.3e3],
            'Total-PTF01-041-40B-B1-061 Sum Normalized Area': [2.1e8, 1.3e3, 2.1e3]
        })
    
    @pytest.fixture 
    def processed_data(self, sample_protein_data):
        """Create processed data with parsed identifiers and descriptions."""
        data = parse_protein_identifiers(sample_protein_data)
        data = parse_gene_and_description(data)
        return data
    
    @pytest.fixture
    def cleaned_sample_names(self):
        """Sample cleaned names mapping."""
        return {
            'Total-PTE01-511-84A-C4-049 Sum Normalized Area': 'E01-511-84A-C4-049',
            'Total-PTE02-Hoof18-050 Sum Normalized Area': 'E02-Hoof18-050',  
            'Total-PTF01-041-40B-B1-061 Sum Normalized Area': 'F01-041-40B-B1-061'
        }
    
    def test_basic_structure_creation(self, processed_data):
        """Test basic data structure creation without sample name cleaning."""
        result = create_standard_data_structure(processed_data)
        
        # Check that we have exactly the right columns in the right order
        expected_annotation_cols = [
            'Protein', 'Description', 'Protein Gene', 'UniProt_Accession', 'UniProt_Entry_Name'
        ]
        
        actual_annotation_cols = list(result.columns[:5])
        assert actual_annotation_cols == expected_annotation_cols
        
        # Check that we have sample columns after annotation columns
        assert len(result.columns) > 5
        sample_columns = list(result.columns[5:])
        assert len(sample_columns) == 3
        
        # Verify sample data is preserved
        assert result.shape[0] == processed_data.shape[0]  # Same number of proteins
        
    def test_structure_with_cleaned_sample_names(self, processed_data, cleaned_sample_names):
        """Test structure creation with cleaned sample names."""
        result = create_standard_data_structure(processed_data, cleaned_sample_names)
        
        # Check annotation columns are correct
        expected_annotation_cols = [
            'Protein', 'Description', 'Protein Gene', 'UniProt_Accession', 'UniProt_Entry_Name'
        ]
        actual_annotation_cols = list(result.columns[:5])
        assert actual_annotation_cols == expected_annotation_cols
        
        # Check that sample columns have cleaned names
        sample_columns = list(result.columns[5:])
        expected_sample_cols = [
            'E01-511-84A-C4-049',
            'E02-Hoof18-050',
            'F01-041-40B-B1-061'
        ]
        assert sample_columns == expected_sample_cols
        
    def test_annotation_column_content(self, processed_data):
        """Test that annotation columns have correct content."""
        result = create_standard_data_structure(processed_data)
        
        # Check Protein column (unchanged)
        assert list(result['Protein']) == [
            'sp|P02787|TRFE_HUMAN',
            'sp|O43240|KLK10_HUMAN', 
            'sp|P06858|LIPL_HUMAN'
        ]
        
        # Check Description column (cleaned)
        assert 'Serotransferrin' in result['Description'].iloc[0]
        assert 'OS=' not in result['Description'].iloc[0]  # Should be cleaned
        
        # Check UniProt_Accession column (parsed)
        assert list(result['UniProt_Accession']) == ['P02787', 'O43240', 'P06858']
        
        # Check UniProt_Entry_Name column (parsed)
        assert list(result['UniProt_Entry_Name']) == ['TRFE_HUMAN', 'KLK10_HUMAN', 'LIPL_HUMAN']
        
        # Check Protein Gene column
        assert list(result['Protein Gene']) == ['TF', 'KLK10', 'LPL']
        
    def test_sample_data_preservation(self, processed_data, cleaned_sample_names):
        """Test that sample data values are preserved correctly."""
        result = create_standard_data_structure(processed_data, cleaned_sample_names)
        
        # Check that sample values are preserved (first protein, first sample)
        original_value = processed_data['Total-PTE01-511-84A-C4-049 Sum Normalized Area'].iloc[0]
        result_value = result['E01-511-84A-C4-049'].iloc[0]
        assert result_value == original_value
        
        # Check another value (second protein, second sample)
        original_value = processed_data['Total-PTE02-Hoof18-050 Sum Normalized Area'].iloc[1]
        result_value = result['E02-Hoof18-050'].iloc[1]
        assert result_value == original_value
        
    def test_missing_annotation_columns_error(self, sample_protein_data):
        """Test that function raises error when required columns are missing."""
        # Try to create structure without parsing (missing UniProt columns)
        with pytest.raises(ValueError, match="Missing required annotation columns"):
            create_standard_data_structure(sample_protein_data)
            
    def test_empty_dataframe_handling(self):
        """Test handling of empty dataframe."""
        # Create processed empty dataframe with required columns
        empty_data = pd.DataFrame(columns=[
            'Protein', 'Description', 'Protein Gene', 'UniProt_Accession', 'UniProt_Entry_Name'
        ])
        
        result = create_standard_data_structure(empty_data)
        assert len(result) == 0
        assert list(result.columns) == [
            'Protein', 'Description', 'Protein Gene', 'UniProt_Accession', 'UniProt_Entry_Name'
        ]
        
    def test_large_dataset_structure(self, processed_data):
        """Test structure with larger number of samples."""
        # Add more sample columns to simulate real dataset
        additional_samples = {}
        for i in range(10, 20):
            col_name = f'Total-PTG{i:02d}-sample-{i:03d} Sum Normalized Area'
            additional_samples[col_name] = np.random.rand(len(processed_data))
            
        large_data = processed_data.copy()
        for col, values in additional_samples.items():
            large_data[col] = values
            
        result = create_standard_data_structure(large_data)
        
        # Should still have exactly 5 annotation columns
        assert list(result.columns[:5]) == [
            'Protein', 'Description', 'Protein Gene', 'UniProt_Accession', 'UniProt_Entry_Name'
        ]
        
        # Should have original 3 + 10 new = 13 sample columns
        assert len(result.columns) == 5 + 3 + 10  # 5 annotation + 13 sample columns
        
    def test_column_order_is_preserved(self, processed_data):
        """Test that the exact column order is maintained."""
        result = create_standard_data_structure(processed_data)
        
        columns = list(result.columns)
        
        # First 5 must be annotation columns in exact order
        assert columns[:5] == [
            'Protein', 'Description', 'Protein Gene', 'UniProt_Accession', 'UniProt_Entry_Name'
        ]
        
        # Remaining should be sample columns
        sample_cols = columns[5:]
        assert len(sample_cols) > 0
        
        # Should not have any annotation columns mixed in with sample columns
        annotation_col_names = [
            'Protein', 'Protein Description', 'Protein Gene', 'Gene', 'Description',
            'UniProt_Accession', 'UniProt_Entry_Name', 'UniProt_Database'
        ]
        
        for col in sample_cols:
            assert col not in annotation_col_names, f"Sample column section contains annotation column: {col}"
            
    def test_supports_peptide_data_structure(self):
        """Test that function can handle peptide data structure (future support)."""
        # Create sample peptide-like data
        peptide_data = pd.DataFrame({
            'Protein': ['sp|P02787|TRFE_HUMAN', 'sp|O43240|KLK10_HUMAN'],
            'Protein Description': [
                'Serotransferrin OS=Homo sapiens OX=9606 GN=TF PE=1 SV=5',
                'Kallikrein-10 OS=Homo sapiens OX=9606 GN=KLK10 PE=1 SV=3'
            ],
            'Protein Gene': ['TF', 'KLK10'],
            'Peptide': ['PEPTIDEONE', 'PEPTIDETWO'],  # Additional peptide-specific column
            'E01-sample-001': [1000, 2000],
            'E02-sample-002': [1500, 2500]
        })
        
        # Process the data
        processed_peptide = parse_protein_identifiers(peptide_data)
        processed_peptide = parse_gene_and_description(processed_peptide)
        
        # Create standard structure - should work and preserve peptide column as sample data
        result = create_standard_data_structure(processed_peptide)
        
        # Should have 5 annotation columns + 3 sample columns (including Peptide)
        assert len(result.columns) == 5 + 3
        assert list(result.columns[:5]) == [
            'Protein', 'Description', 'Protein Gene', 'UniProt_Accession', 'UniProt_Entry_Name'
        ]
        
        # Peptide column should be treated as sample data
        assert 'Peptide' in result.columns[5:]
