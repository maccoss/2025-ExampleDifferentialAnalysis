"""
Tests for export module
"""

import pandas as pd
import os
import tempfile
from proteomics_toolkit.export import export_analysis_results


class TestExportModule:
    """Test the export module functionality."""

    def test_sample_metadata_export_has_column_header(self):
        """Test that sample metadata export includes a proper column header for sample IDs."""
        
        # Create test data
        normalized_data = pd.DataFrame({
            'Protein': ['P1', 'P2', 'P3'],
            'Sample_A': [1.0, 2.0, 3.0],
            'Sample_B': [1.5, 2.5, 3.5],
        })
        
        sample_metadata = {
            'Sample_A': {
                'Sample Type': 'Unknown',
                'Group': 'Control',
                'Subject': '001'
            },
            'Sample_B': {
                'Sample Type': 'Unknown', 
                'Group': 'Treatment',
                'Subject': '002'
            }
        }
        
        # Create a temporary directory for output
        with tempfile.TemporaryDirectory() as temp_dir:
            output_prefix = os.path.join(temp_dir, "test_export")
            
            # Run the export
            exported_files = export_analysis_results(
                normalized_data=normalized_data,
                sample_metadata=sample_metadata,
                output_prefix=output_prefix
            )
            
            # Check that metadata file was created
            assert "sample_metadata" in exported_files
            metadata_file = exported_files["sample_metadata"]
            assert os.path.exists(metadata_file)
            
            # Read the exported metadata file
            exported_metadata = pd.read_csv(metadata_file, index_col=0)
            
            # Check that the index has a proper name/header
            assert exported_metadata.index.name == "Sample_ID"
            
            # Check that the data is correct
            assert len(exported_metadata) == 2
            assert 'Sample_A' in exported_metadata.index
            assert 'Sample_B' in exported_metadata.index
            assert 'Sample Type' in exported_metadata.columns
            assert 'Group' in exported_metadata.columns
            assert 'Subject' in exported_metadata.columns
            
            # Check values
            assert exported_metadata.loc['Sample_A', 'Group'] == 'Control'
            assert exported_metadata.loc['Sample_B', 'Group'] == 'Treatment'

    def test_sample_metadata_csv_format(self):
        """Test that the CSV format of sample metadata has the correct header structure."""
        
        sample_metadata = {
            'E01-511-84A-C4-049': {
                'Sample Type': 'Unknown',
                'Group': '80',
                'Comparison': 'Drug',
                'Subject': '84'
            },
            'E02-Hoof18-050': {
                'Sample Type': 'Quality Control',
                'Group': 'HoofPool',
                'Comparison': '',
                'Subject': 'HoofPool'
            }
        }
        
        normalized_data = pd.DataFrame({
            'Protein': ['P1', 'P2'],
            'E01-511-84A-C4-049': [1.0, 2.0],
            'E02-Hoof18-050': [1.5, 2.5],
        })
        
        with tempfile.TemporaryDirectory() as temp_dir:
            output_prefix = os.path.join(temp_dir, "test_format")
            
            exported_files = export_analysis_results(
                normalized_data=normalized_data,
                sample_metadata=sample_metadata,
                output_prefix=output_prefix
            )
            
            # Read the raw CSV file to check the header line
            metadata_file = exported_files["sample_metadata"]
            with open(metadata_file, 'r') as f:
                header_line = f.readline().strip()
            
            # The header should start with "Sample_ID" followed by the metadata columns
            expected_start = "Sample_ID,"
            assert header_line.startswith(expected_start), f"Header should start with '{expected_start}', but got: {header_line}"
            
            # Should contain all the expected column names
            assert "Sample Type" in header_line
            assert "Group" in header_line
            assert "Comparison" in header_line
            assert "Subject" in header_line
