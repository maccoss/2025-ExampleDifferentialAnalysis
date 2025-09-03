"""
Tests for proteomics_toolkit.data_import module
"""

import pandas as pd

from proteomics_toolkit.data_import import (
    load_skyline_data,
    parse_uniprot_identifier,
    parse_gene_from_description,
    clean_description,
    identify_sample_columns,
    clean_sample_names,
    match_samples_to_metadata,
    identify_and_classify_controls,
)


class TestLoadSkylineData:
    """Test the main data loading function"""

    def test_load_skyline_data_success(self, temp_csv_files):
        """Test successful loading of protein data and metadata"""
        protein_file, metadata_file = temp_csv_files

        result = load_skyline_data(
            protein_file=protein_file, metadata_file=metadata_file, peptide_file=None
        )

        assert len(result) == 3  # protein_data, metadata, peptide_data
        protein_data, metadata, peptide_data = result

        assert isinstance(protein_data, pd.DataFrame)
        assert isinstance(metadata, pd.DataFrame)
        assert peptide_data is None  # Since we didn't provide peptide file

        # Check protein data structure
        assert "Protein" in protein_data.columns
        assert "Sample_1" in protein_data.columns
        assert len(protein_data) == 3  # 3 proteins

        # Check metadata structure
        assert "Sample" in metadata.columns
        assert "Subject" in metadata.columns
        assert len(metadata) == 3  # 3 samples


class TestParseUniprotIdentifier:
    """Test UniProt identifier parsing"""

    def test_parse_sp_identifier(self):
        """Test parsing SwissProt identifier"""
        result = parse_uniprot_identifier("sp|P12345|PROT1_HUMAN")

        assert result["database"] == "SwissProt"
        assert result["accession"] == "P12345"
        assert result["entry_name"] == "PROT1_HUMAN"

    def test_parse_tr_identifier(self):
        """Test parsing TrEMBL identifier"""
        result = parse_uniprot_identifier("tr|Q67890|Q67890_HUMAN")

        assert result["database"] == "TrEMBL"
        assert result["accession"] == "Q67890"
        assert result["entry_name"] == "Q67890_HUMAN"

    def test_parse_invalid_identifier(self):
        """Test parsing invalid identifier"""
        result = parse_uniprot_identifier("INVALID_ID")

        # Function tries to extract accession pattern and finds "INVALID"
        assert result["database"] == ""
        assert result["accession"] == "INVALID"  # Function extracts this part
        assert result["entry_name"] == ""

    def test_parse_completely_invalid_identifier(self):
        """Test parsing completely invalid identifier"""
        result = parse_uniprot_identifier("xyz123")

        assert result["database"] == ""
        assert result["accession"] == ""
        assert result["entry_name"] == ""


class TestParseGeneFromDescription:
    """Test gene name extraction from protein descriptions"""

    def test_parse_gene_with_gn(self):
        """Test extracting gene name with GN= format"""
        desc = "Protein kinase B GN=AKT1 PE=1 SV=1"
        result = parse_gene_from_description(desc)
        assert result == "AKT1"

    def test_parse_gene_with_gene(self):
        """Test extracting gene name with Gene= format"""
        desc = "Protein description Gene=MYC PE=1 SV=2"
        result = parse_gene_from_description(desc)
        assert result == "MYC"

    def test_parse_gene_no_gene_info(self):
        """Test when no gene information is present"""
        desc = "Protein description without gene info PE=1 SV=1"
        result = parse_gene_from_description(desc)
        assert result == ""

    def test_parse_gene_empty_description(self):
        """Test with empty description"""
        result = parse_gene_from_description("")
        assert result == ""


class TestCleanDescription:
    """Test protein description cleaning"""

    def test_clean_description_with_flags(self):
        """Test cleaning description with PE and SV flags"""
        desc = "Protein kinase B GN=AKT1 PE=1 SV=1"
        result = clean_description(desc)
        assert result == "Protein kinase B"

    def test_clean_description_with_os(self):
        """Test cleaning description with organism info"""
        desc = "Protein kinase B OS=Homo sapiens GN=AKT1"
        result = clean_description(desc)
        assert result == "Protein kinase B"

    def test_clean_description_no_flags(self):
        """Test cleaning description without flags"""
        desc = "Simple protein description"
        result = clean_description(desc)
        assert result == "Simple protein description"


class TestIdentifySampleColumns:
    """Test sample column identification"""

    def test_identify_sample_columns_basic(self):
        """Test basic sample column identification"""
        protein_data = pd.DataFrame(
            {
                "Protein": ["P1", "P2"],
                "Gene": ["G1", "G2"],
                "Sample_1": [100, 200],
                "Sample_2": [150, 250],
                "Sample_3": [120, 180],
            }
        )

        metadata = pd.DataFrame({"Sample": ["Sample_1", "Sample_2", "Sample_3"]})

        result = identify_sample_columns(protein_data, metadata)
        expected = ["Sample_1", "Sample_2", "Sample_3"]
        assert result == expected

    def test_identify_sample_columns_partial_match(self):
        """Test when some samples don't match metadata"""
        protein_data = pd.DataFrame(
            {
                "Protein": ["P1", "P2"],
                "Sample_1": [100, 200],
                "Sample_2": [150, 250],
                "Sample_Unknown": [120, 180],  # This won't match metadata
            }
        )

        metadata = pd.DataFrame({"Sample": ["Sample_1", "Sample_2"]})

        result = identify_sample_columns(protein_data, metadata)
        expected = ["Sample_1", "Sample_2"]
        assert result == expected


class TestCleanSampleNames:
    """Test sample name cleaning"""

    def test_clean_sample_names_with_prefix(self):
        """Test removing common prefix from sample names"""
        sample_columns = ["Prefix_Sample_1", "Prefix_Sample_2", "Prefix_Sample_3"]

        result = clean_sample_names(sample_columns, common_prefix="Prefix_")

        expected = {
            "Prefix_Sample_1": "Sample_1",
            "Prefix_Sample_2": "Sample_2",
            "Prefix_Sample_3": "Sample_3",
        }
        assert result == expected

    def test_clean_sample_names_with_suffix(self):
        """Test removing common suffix from sample names"""
        sample_columns = ["Sample_1_Suffix", "Sample_2_Suffix", "Sample_3_Suffix"]

        result = clean_sample_names(sample_columns, common_suffix="_Suffix")

        expected = {
            "Sample_1_Suffix": "Sample_1",
            "Sample_2_Suffix": "Sample_2",
            "Sample_3_Suffix": "Sample_3",
        }
        assert result == expected

    def test_clean_sample_names_no_changes(self):
        """Test when no cleaning is needed"""
        sample_columns = ["Sample_1", "Sample_2", "Sample_3"]

        result = clean_sample_names(sample_columns)

        expected = {name: name for name in sample_columns}
        assert result == expected

    def test_sample_metadata_mapping_consistency(self):
        """
        Test that sample cleaning doesn't break metadata mapping.
        This addresses the specific issue found in the notebook.
        """
        # Simulate the real-world scenario from the notebook
        original_samples = [
            "Total-PTE01-511-84A-C4-049 Sum Normalized Area",
            "Total-PTE02-Hoof18-050 Sum Normalized Area"
        ]
        
        # Clean the sample names
        result = clean_sample_names(original_samples, auto_detect=True)
        
        # Verify that the result is a dictionary mapping original to cleaned names
        assert isinstance(result, dict)
        assert len(result) == len(original_samples)
        
        # Verify that all original samples are keys
        for sample in original_samples:
            assert sample in result
            
        # Verify that cleaned names are different (shorter)
        for original, cleaned in result.items():
            assert len(cleaned) < len(original)
            assert isinstance(cleaned, str)
            assert len(cleaned.strip()) > 0  # Not empty
            
    def test_metadata_synchronization_workflow(self):
        """
        Test the complete workflow of cleaning sample names and updating metadata.
        This replicates the fix implemented in the notebook.
        """
        # Original data structure from notebook
        original_samples = [
            "Total-PTE01-511-84A-C4-049 Sum Normalized Area",
            "Total-PTE02-Hoof18-050 Sum Normalized Area"
        ]
        
        # Original metadata (keyed by replicate name, not sample column name)
        original_metadata = {
            "Total-PTE01-511-84A-C4-049": {"Group": 80, "Subject": "84"},
            "Total-PTE02-Hoof18-050": {"Group": "HoofPool", "Subject": "HoofPool"}
        }
        
        # Step 1: Clean sample names
        cleaned_mapping = clean_sample_names(original_samples, auto_detect=True)
        
        # Step 2: Update metadata to use cleaned names as keys (the fix)
        updated_metadata = {}
        for original_sample, cleaned_sample in cleaned_mapping.items():
            # Extract replicate name by removing the suffix
            replicate_name = original_sample.replace(' Sum Normalized Area', '')
            if replicate_name in original_metadata:
                updated_metadata[cleaned_sample] = original_metadata[replicate_name]
        
        # Verify the fix worked
        cleaned_samples = list(cleaned_mapping.values())
        
        # All cleaned samples should have metadata
        assert len(updated_metadata) == len(cleaned_samples)
        for sample in cleaned_samples:
            assert sample in updated_metadata
            assert 'Group' in updated_metadata[sample]
            assert 'Subject' in updated_metadata[sample]


class TestMatchSamplesToMetadata:
    """Test sample-to-metadata matching"""

    def test_match_samples_exact(self):
        """Test exact matching of samples to metadata"""
        cleaned_sample_names = {"Sample_1": "Sample_1", "Sample_2": "Sample_2"}

        metadata = pd.DataFrame(
            {"Sample": ["Sample_1", "Sample_2"], "Group": ["Control", "Treatment"]}
        )

        result = match_samples_to_metadata(cleaned_sample_names, metadata)

        assert len(result) == 2
        assert "Sample_1" in result
        assert "Sample_2" in result
        assert result["Sample_1"]["Group"] == "Control"
        assert result["Sample_2"]["Group"] == "Treatment"

    def test_match_samples_partial(self):
        """Test when some samples don't match metadata"""
        cleaned_sample_names = {
            "Sample_1": "Sample_1",
            "Sample_Unknown": "Sample_Unknown",
        }

        metadata = pd.DataFrame({"Sample": ["Sample_1"], "Group": ["Control"]})

        result = match_samples_to_metadata(cleaned_sample_names, metadata, include_unmatched=False)

        assert len(result) == 1  # Only matched sample included
        assert "Sample_1" in result
        assert "Sample_Unknown" not in result


class TestIdentifyAndClassifyControls:
    """Test control sample identification"""

    def test_identify_controls_from_names(self, sample_metadata):
        """Test identifying controls from sample names"""
        # Add some control samples to metadata
        control_metadata = sample_metadata.copy()
        control_metadata["Pool_1"] = {"Sample_Type": "Pool"}
        control_metadata["QC_1"] = {"Sample_Type": "QC"}

        updated_metadata, summary = identify_and_classify_controls(control_metadata)

        assert isinstance(updated_metadata, dict)
        assert isinstance(summary, dict)
        assert "total_updated" in summary
        assert "control_types" in summary

    def test_identify_controls_empty_metadata(self):
        """Test control identification with empty metadata"""
        updated_metadata, summary = identify_and_classify_controls({})

        assert isinstance(updated_metadata, dict)
        assert isinstance(summary, dict)
        assert summary["total_updated"] == 0
