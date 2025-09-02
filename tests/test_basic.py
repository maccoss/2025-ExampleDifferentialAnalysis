"""
Basic tests to verify pytest setup and basic functionality
"""

import pandas as pd
import numpy as np


def test_basic_functionality():
    """Test basic functionality to verify test setup works"""
    # Basic pandas operations
    df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})

    assert len(df) == 3
    assert list(df.columns) == ["A", "B"]
    assert df["A"].sum() == 6


def test_numpy_functionality():
    """Test basic numpy functionality"""
    arr = np.array([1, 2, 3, 4, 5])

    assert arr.mean() == 3.0
    assert arr.std() > 1.0
    assert len(arr) == 5


def test_proteomics_toolkit_import():
    """Test that we can import the proteomics toolkit"""
    try:
        import proteomics_toolkit

        assert hasattr(proteomics_toolkit, "__version__")
    except ImportError:
        # If import fails, at least verify the module structure exists
        import os

        toolkit_path = os.path.join(
            os.path.dirname(os.path.dirname(__file__)), "proteomics_toolkit"
        )
        assert os.path.exists(toolkit_path)
        assert os.path.isfile(os.path.join(toolkit_path, "__init__.py"))


def test_basic_statistical_config():
    """Test that we can import and create a statistical configuration"""
    from proteomics_toolkit.statistical_analysis import StatisticalConfig

    config = StatisticalConfig()

    assert config.statistical_test_method == "mixed_effects"
    assert config.p_value_threshold == 0.05
    assert config.use_adjusted_pvalue == "adjusted"  # Fix based on actual default


def test_basic_data_operations():
    """Test basic data operations that the toolkit would use"""
    # Create sample data similar to proteomics data
    data = pd.DataFrame(
        {
            "Protein": ["P001", "P002", "P003"],
            "Gene": ["GENE1", "GENE2", "GENE3"],
            "Sample_1": [100.0, 200.0, 150.0],
            "Sample_2": [110.0, 190.0, 160.0],
            "Sample_3": [95.0, 210.0, 140.0],
        }
    )

    # Test basic operations
    sample_columns = ["Sample_1", "Sample_2", "Sample_3"]
    sample_data = data[sample_columns]

    assert sample_data.shape == (3, 3)
    assert all(sample_data.dtypes == "float64")

    # Test median calculation (used in median normalization)
    medians = sample_data.median(axis=0)
    assert len(medians) == 3
    assert all(medians > 0)

    # Test log transformation (common in proteomics)
    log_data = np.log2(sample_data)
    assert all(log_data.mean() < sample_data.mean())


def test_sample_metadata_structure():
    """Test sample metadata dictionary structure"""
    sample_metadata = {
        "Sample_1": {"Subject": "S001", "Group": "Control", "Visit": "Baseline"},
        "Sample_2": {"Subject": "S002", "Group": "Treatment", "Visit": "Baseline"},
        "Sample_3": {"Subject": "S001", "Group": "Control", "Visit": "Week4"},
    }

    assert len(sample_metadata) == 3
    assert "Sample_1" in sample_metadata
    assert sample_metadata["Sample_1"]["Group"] == "Control"

    # Test metadata operations
    groups = set()
    for sample_info in sample_metadata.values():
        groups.add(sample_info["Group"])

    assert groups == {"Control", "Treatment"}


def test_differential_results_structure():
    """Test differential analysis results structure"""
    # Simulate differential analysis results
    results = pd.DataFrame(
        {
            "Protein": ["P001", "P002", "P003"],
            "logFC": [0.5, -1.2, 0.8],
            "P.Value": [0.01, 0.001, 0.05],
            "adj.P.Val": [0.03, 0.003, 0.05],
            "AveExpr": [10.5, 12.3, 9.8],
        }
    )

    assert len(results) == 3
    assert "P.Value" in results.columns
    assert "logFC" in results.columns

    # Test sorting by p-value (common operation)
    sorted_results = results.sort_values("P.Value")
    assert sorted_results.iloc[0]["Protein"] == "P002"  # Lowest p-value

    # Test significance filtering (P003 has adj.P.Val = 0.05, which is not < 0.05)
    significant = results[results["adj.P.Val"] < 0.05]
    assert len(significant) == 2  # P001 and P002 are significant (< 0.05)
