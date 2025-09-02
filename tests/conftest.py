"""
Pytest configuration and fixtures for proteomics_toolkit tests
"""

import pytest
import pandas as pd
import numpy as np
import tempfile
import os
from proteomics_toolkit.statistical_analysis import StatisticalConfig


@pytest.fixture
def sample_protein_data():
    """Create sample protein quantification data for testing"""
    # Create 20 proteins x 12 samples
    np.random.seed(42)

    protein_names = [f"P{i:05d}" for i in range(20)]
    sample_names = [
        "Sample_A_1",
        "Sample_A_2",
        "Sample_A_3",
        "Sample_B_1",
        "Sample_B_2",
        "Sample_B_3",
        "Sample_C_1",
        "Sample_C_2",
        "Sample_C_3",
        "Sample_D_1",
        "Sample_D_2",
        "Sample_D_3",
    ]

    # Generate realistic protein abundance data (log-scale)
    data_matrix = np.random.normal(20, 2, (len(protein_names), len(sample_names)))
    # Add some missing values
    mask = np.random.random((len(protein_names), len(sample_names))) < 0.05
    data_matrix[mask] = np.nan

    df = pd.DataFrame(data_matrix, index=protein_names, columns=sample_names)

    # Add some annotation columns
    df.insert(0, "Protein", protein_names)
    df.insert(1, "ProteinName", [f"Protein_{i}" for i in range(20)])
    df.insert(2, "Gene", [f"GENE{i}" for i in range(20)])

    return df


@pytest.fixture
def standardized_protein_data():
    """Create sample protein data with standardized structure for testing normalization functions"""
    # Create 20 proteins x 12 samples
    np.random.seed(42)

    protein_names = [f"P{i:05d}" for i in range(20)]
    sample_names = [
        "Sample_A_1",
        "Sample_A_2",
        "Sample_A_3",
        "Sample_B_1",
        "Sample_B_2",
        "Sample_B_3",
        "Sample_C_1",
        "Sample_C_2",
        "Sample_C_3",
        "Sample_D_1",
        "Sample_D_2",
        "Sample_D_3",
    ]

    # Generate realistic protein abundance data (log-scale)
    data_matrix = np.random.normal(20, 2, (len(protein_names), len(sample_names)))
    # Add some missing values
    mask = np.random.random((len(protein_names), len(sample_names))) < 0.05
    data_matrix[mask] = np.nan

    # Create DataFrame with EXACT standardized structure
    df = pd.DataFrame(
        {
            # EXACTLY the 5 required annotation columns in order:
            "Protein": protein_names,
            "Description": [f"Protein_{i}_description" for i in range(20)],
            "Protein Gene": [f"GENE{i}" for i in range(20)],
            "UniProt_Accession": [f"P{i:05d}_ACC" for i in range(20)],
            "UniProt_Entry_Name": [f"PROT{i}_HUMAN" for i in range(20)],
            # Sample columns:
            **{
                sample_name: data_matrix[:, i]
                for i, sample_name in enumerate(sample_names)
            },
        }
    )

    return df


@pytest.fixture
def sample_metadata():
    """Create sample metadata for testing"""
    metadata_dict = {
        "Sample_A_1": {
            "Subject": "S001",
            "Group": "Control",
            "Visit": "Baseline",
            "DrugDose": "Placebo",
        },
        "Sample_A_2": {
            "Subject": "S002",
            "Group": "Control",
            "Visit": "Baseline",
            "DrugDose": "Placebo",
        },
        "Sample_A_3": {
            "Subject": "S003",
            "Group": "Control",
            "Visit": "Baseline",
            "DrugDose": "Placebo",
        },
        "Sample_B_1": {
            "Subject": "S001",
            "Group": "Control",
            "Visit": "Week4",
            "DrugDose": "Placebo",
        },
        "Sample_B_2": {
            "Subject": "S002",
            "Group": "Control",
            "Visit": "Week4",
            "DrugDose": "Placebo",
        },
        "Sample_B_3": {
            "Subject": "S003",
            "Group": "Control",
            "Visit": "Week4",
            "DrugDose": "Placebo",
        },
        "Sample_C_1": {
            "Subject": "S004",
            "Group": "Treatment",
            "Visit": "Baseline",
            "DrugDose": "High",
        },
        "Sample_C_2": {
            "Subject": "S005",
            "Group": "Treatment",
            "Visit": "Baseline",
            "DrugDose": "High",
        },
        "Sample_C_3": {
            "Subject": "S006",
            "Group": "Treatment",
            "Visit": "Baseline",
            "DrugDose": "High",
        },
        "Sample_D_1": {
            "Subject": "S004",
            "Group": "Treatment",
            "Visit": "Week4",
            "DrugDose": "High",
        },
        "Sample_D_2": {
            "Subject": "S005",
            "Group": "Treatment",
            "Visit": "Week4",
            "DrugDose": "High",
        },
        "Sample_D_3": {
            "Subject": "S006",
            "Group": "Treatment",
            "Visit": "Week4",
            "DrugDose": "High",
        },
    }
    return metadata_dict


@pytest.fixture
def sample_metadata_df():
    """Create sample metadata DataFrame for testing"""
    data = {
        "Sample": [
            "Sample_A_1",
            "Sample_A_2",
            "Sample_A_3",
            "Sample_B_1",
            "Sample_B_2",
            "Sample_B_3",
            "Sample_C_1",
            "Sample_C_2",
        ],
        "Subject": ["S001", "S002", "S003", "S001", "S002", "S003", "S004", "S005"],
        "Group": [
            "Control",
            "Control",
            "Control",
            "Control",
            "Control",
            "Control",
            "Treatment",
            "Treatment",
        ],
        "Visit": [
            "Baseline",
            "Baseline",
            "Baseline",
            "Week4",
            "Week4",
            "Week4",
            "Baseline",
            "Baseline",
        ],
        "DrugDose": [
            "Placebo",
            "Placebo",
            "Placebo",
            "Placebo",
            "Placebo",
            "Placebo",
            "High",
            "High",
        ],
    }
    return pd.DataFrame(data)


@pytest.fixture
def statistical_config():
    """Create sample statistical configuration"""
    config = StatisticalConfig()
    config.statistical_test_method = "mixed_effects"
    config.analysis_type = "interaction_analysis"
    config.p_value_threshold = 0.05
    config.fold_change_threshold = 1.5
    config.subject_column = "Subject"
    config.paired_column = "Visit"
    config.paired_label1 = "Baseline"
    config.paired_label2 = "Week4"
    config.group_column = "Group"
    config.group_labels = ["Control", "Treatment"]
    config.interaction_terms = ["DrugDose", "Visit"]
    config.additional_interactions = []
    config.covariates = []
    config.use_adjusted_pvalue = "fdr_bh"
    config.enable_pvalue_fallback = True
    return config


@pytest.fixture
def sample_columns():
    """Sample column names for testing"""
    return [
        "Sample_A_1",
        "Sample_A_2",
        "Sample_A_3",
        "Sample_B_1",
        "Sample_B_2",
        "Sample_B_3",
        "Sample_C_1",
        "Sample_C_2",
        "Sample_C_3",
        "Sample_D_1",
        "Sample_D_2",
        "Sample_D_3",
    ]


@pytest.fixture
def temp_csv_files():
    """Create temporary CSV files for testing file I/O operations"""
    # Create temporary directory
    temp_dir = tempfile.mkdtemp()

    # Create sample protein data file
    protein_data = {
        "Protein": ["P00001", "P00002", "P00003"],
        "ProteinName": ["Protein_1", "Protein_2", "Protein_3"],
        "Gene": ["GENE1", "GENE2", "GENE3"],
        "Sample_1": [100.0, 200.0, 150.0],
        "Sample_2": [110.0, 180.0, 140.0],
        "Sample_3": [95.0, 220.0, 160.0],
    }
    protein_file = os.path.join(temp_dir, "protein_data.csv")
    pd.DataFrame(protein_data).to_csv(protein_file, index=False)

    # Create sample metadata file
    metadata_data = {
        "Sample": ["Sample_1", "Sample_2", "Sample_3"],
        "Subject": ["S001", "S002", "S003"],
        "Group": ["Control", "Treatment", "Control"],
        "Visit": ["Baseline", "Baseline", "Baseline"],
    }
    metadata_file = os.path.join(temp_dir, "metadata.csv")
    pd.DataFrame(metadata_data).to_csv(metadata_file, index=False)

    yield protein_file, metadata_file

    # Cleanup
    import shutil

    shutil.rmtree(temp_dir)


@pytest.fixture
def differential_results():
    """Create sample differential analysis results"""
    np.random.seed(42)
    n_proteins = 10

    data = {
        "Protein": [f"P{i:05d}" for i in range(n_proteins)],
        "logFC": np.random.normal(0, 1, n_proteins),
        "AveExpr": np.random.normal(20, 2, n_proteins),
        "t": np.random.normal(0, 2, n_proteins),
        "P.Value": np.random.uniform(0.001, 0.8, n_proteins),
        "adj.P.Val": np.random.uniform(0.001, 0.9, n_proteins),
        "B": np.random.normal(0, 1, n_proteins),
        "test_method": ["Mixed-effects model"] * n_proteins,
    }

    df = pd.DataFrame(data)
    df = df.sort_values("P.Value")  # Sort by p-value like real results
    return df


@pytest.fixture
def protein_identifiers():
    """Sample protein identifiers for testing parsing"""
    return [
        "sp|P12345|PROT1_HUMAN",
        "tr|Q67890|Q67890_HUMAN",
        "sp|P11111|PROT2_HUMAN",
        "INVALID_ID",
        "sp|P22222|PROT3_HUMAN",
    ]
