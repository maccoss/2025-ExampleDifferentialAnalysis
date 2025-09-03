#!/usr/bin/env python
"""Test the flexible configuration export functionality."""

import os
import tempfile
from proteomics_toolkit.export import create_config_dict_from_notebook_vars, export_timestamped_config


class TestFlexibleConfigExport:
    """Test flexible configuration dictionary creation and export."""
    
    def test_minimal_config_only_includes_provided_params(self):
        """Test that minimal config only includes provided parameters."""
        # Create minimal config with only essentials
        config = create_config_dict_from_notebook_vars(
            metadata_file="data/metadata.csv",
            protein_file="data/proteins.csv",
            normalization_method="Median",
            min_detection_rate=0.5
        )
        
        # Check that only provided params (plus minimal defaults) are included
        assert config["toolkit_path"] == "."  # Default
        assert config["metadata_file"] == "data/metadata.csv"
        assert config["protein_file"] == "data/proteins.csv"
        assert config["normalization_method"] == "Median"
        assert config["min_detection_rate"] == 0.5
        
        # Check that unprovided params are NOT included
        assert "subject_column" not in config
        assert "paired_column" not in config
        assert "group_labels" not in config
        assert "interaction_terms" not in config
        assert "control_labels" not in config
    
    def test_paired_analysis_config(self):
        """Test configuration for paired analysis."""
        config = create_config_dict_from_notebook_vars(
            metadata_file="data/metadata.csv",
            protein_file="data/proteins.csv",
            normalization_method="VSN",
            statistical_test_method="paired_t",
            paired_column="Visit",
            paired_label1="Baseline",
            paired_label2="Treatment",
            p_value_threshold=0.05
        )
        
        # Check paired analysis params are included
        assert config["statistical_test_method"] == "paired_t"
        assert config["paired_column"] == "Visit"
        assert config["paired_label1"] == "Baseline"
        assert config["paired_label2"] == "Treatment"
        assert config["p_value_threshold"] == 0.05
        
        # Check unprovided params are still not included
        assert "subject_column" not in config
        assert "group_column" not in config
        assert "covariates" not in config
    
    def test_unpaired_analysis_config(self):
        """Test configuration for simple unpaired analysis."""
        config = create_config_dict_from_notebook_vars(
            metadata_file="data/metadata.csv",
            protein_file="data/proteins.csv",
            normalization_method="Quantile",
            statistical_test_method="welch_t",
            group_column="Treatment",
            group_labels=["Control", "Drug"]
        )
        
        # Check group analysis params
        assert config["group_column"] == "Treatment"
        assert config["group_labels"] == ["Control", "Drug"]
        
        # Paired params should not be included
        assert "paired_column" not in config
        assert "paired_label1" not in config
        assert "paired_label2" not in config
    
    def test_visualization_only_config(self):
        """Test configuration with only visualization settings."""
        config = create_config_dict_from_notebook_vars(
            use_systematic_colors=True,
            systematic_color_palette="Dark2",
            label_top_proteins=10
        )
        
        # Check viz params
        assert config["use_systematic_colors"] is True
        assert config["systematic_color_palette"] == "Dark2"
        assert config["label_top_proteins"] == 10
        
        # Should still have minimal defaults
        assert config["toolkit_path"] == "."
        
        # But no file or analysis params
        assert "metadata_file" not in config
        assert "protein_file" not in config
        assert "normalization_method" not in config
    
    def test_mixed_effects_config(self):
        """Test configuration for mixed effects analysis."""
        config = create_config_dict_from_notebook_vars(
            metadata_file="data/metadata.csv",
            protein_file="data/proteins.csv",
            statistical_test_method="mixed_effects",
            subject_column="Subject",
            group_column="DrugDose",
            group_labels=["0", "20", "40", "80"],
            interaction_terms=["DrugDose", "Visit"],
            covariates=["Age", "Sex"]
        )
        
        # Check mixed effects params
        assert config["statistical_test_method"] == "mixed_effects"
        assert config["subject_column"] == "Subject"
        assert config["interaction_terms"] == ["DrugDose", "Visit"]
        assert config["covariates"] == ["Age", "Sex"]
        
        # Other analysis types' params should not be included
        assert "paired_label1" not in config
        assert "paired_label2" not in config
    
    def test_empty_config(self):
        """Test that empty config only has minimal defaults."""
        config = create_config_dict_from_notebook_vars()
        
        # Should only have the absolute minimum
        assert config == {"toolkit_path": "."}
    
    def test_overriding_defaults(self):
        """Test that provided values override defaults."""
        config = create_config_dict_from_notebook_vars(
            toolkit_path="/custom/path/to/toolkit"
        )
        
        assert config["toolkit_path"] == "/custom/path/to/toolkit"
    
    def test_timestamped_config_export_minimal(self):
        """Test that timestamped config export only includes provided sections."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create minimal config
            config = create_config_dict_from_notebook_vars(
                metadata_file="data/metadata.csv",
                protein_file="data/proteins.csv",
                min_detection_rate=0.5,
                p_value_threshold=0.05
            )
            
            # Export it
            config_file = export_timestamped_config(
                config_dict=config,
                output_prefix=os.path.join(tmpdir, "test"),
                analysis_description="Test minimal config"
            )
            
            # Read the file
            with open(config_file, 'r') as f:
                content = f.read()
            
            # Check that only relevant sections are included
            assert "INPUT FILES AND PATHS" in content
            assert "DATA FILTERING PARAMETERS" in content
            assert "SIGNIFICANCE THRESHOLDS" in content
            
            # Check that irrelevant sections are NOT included
            assert "EXPERIMENTAL DESIGN CONFIGURATION" not in content
            assert "MIXED-EFFECTS MODEL CONFIGURATION" not in content
            assert "CONTROL SAMPLE CONFIGURATION" not in content
            
            # Check that only provided parameters are written
            assert "metadata_file = 'data/metadata.csv'" in content
            assert "min_detection_rate = 0.5" in content
            assert "p_value_threshold = 0.05" in content
            
            # Check that unprovided parameters are NOT written
            assert "subject_column" not in content
            assert "paired_column" not in content
            assert "group_labels" not in content
    
    def test_timestamped_config_export_comprehensive(self):
        """Test timestamped config export with many parameters."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create comprehensive config
            config = create_config_dict_from_notebook_vars(
                # Files
                metadata_file="data/metadata.csv",
                protein_file="data/proteins.csv",
                peptide_file="data/peptides.csv",
                
                # Normalization
                normalization_method="VSN",
                optimize_vsn=True,
                
                # Stats
                statistical_test_method="mixed_effects",
                subject_column="Subject",
                interaction_terms=["DrugDose", "Visit"],
                
                # Viz
                use_systematic_colors=True,
                systematic_color_palette="Set1"
            )
            
            # Export it
            config_file = export_timestamped_config(
                config_dict=config,
                output_prefix=os.path.join(tmpdir, "test"),
                analysis_description="Test comprehensive config"
            )
            
            # Read the file
            with open(config_file, 'r') as f:
                content = f.read()
            
            # Check multiple sections are included
            assert "INPUT FILES AND PATHS" in content
            assert "NORMALIZATION STRATEGY" in content
            assert "STATISTICAL ANALYSIS STRATEGY" in content
            assert "EXPERIMENTAL DESIGN CONFIGURATION" in content
            assert "MIXED-EFFECTS MODEL CONFIGURATION" in content
            assert "VISUALIZATION SETTINGS" in content
            
            # Check parameters are written
            assert "peptide_file = 'data/peptides.csv'" in content
            assert "optimize_vsn = True" in content
            assert "subject_column = 'Subject'" in content
            assert "interaction_terms = ['DrugDose', 'Visit']" in content
    
    def test_config_section_numbering_is_dynamic(self):
        """Test that section numbering adjusts based on included sections."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Config with only a few sections
            config = create_config_dict_from_notebook_vars(
                normalization_method="Median",  # Section: NORMALIZATION
                p_value_threshold=0.05,         # Section: SIGNIFICANCE
                output_prefix="test"            # Section: OUTPUT
            )
            
            config_file = export_timestamped_config(
                config_dict=config,
                output_prefix=os.path.join(tmpdir, "test"),
                analysis_description="Test section numbering"
            )
            
            with open(config_file, 'r') as f:
                content = f.read()
            
            # Check that sections are numbered sequentially
            assert "# 1. INPUT FILES AND PATHS" in content  # toolkit_path is always included
            assert "# 2. NORMALIZATION STRATEGY" in content
            assert "# 3. SIGNIFICANCE THRESHOLDS" in content
            assert "# 4. OUTPUT AND EXPORT SETTINGS" in content
            
            # There should be no gaps in numbering
            assert "# 5." not in content  # No 5th section
    
    def test_computed_values_as_comments(self):
        """Test that computed values are included as comments."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config = create_config_dict_from_notebook_vars(
                metadata_file="data/metadata.csv",
                group_colors={"Control": "#1f77b4", "Treatment": "#ff7f0e"}
            )
            
            computed_values = {
                "Total proteins analyzed": 1500,
                "Total samples": 24,
                "group_colors": {"Control": "#1f77b4", "Treatment": "#ff7f0e"}
            }
            
            config_file = export_timestamped_config(
                config_dict=config,
                output_prefix=os.path.join(tmpdir, "test"),
                analysis_description="Test computed values",
                computed_values=computed_values
            )
            
            with open(config_file, 'r') as f:
                content = f.read()
            
            # Check computed values section
            assert "# COMPUTED VALUES (for reference)" in content
            assert "# Total proteins analyzed: 1500" in content
            assert "# Total samples: 24" in content
            assert "# Group colors assigned:" in content
            assert "#   Control: #1f77b4" in content
            assert "#   Treatment: #ff7f0e" in content
