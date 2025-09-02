# Proteomics Toolkit Test Suite

This directory contains comprehensive tests for the proteomics_toolkit package using pytest.

## Test Structure

```
tests/
├── conftest.py                    # Pytest configuration and shared fixtures
├── test_basic.py                  # Basic functionality tests
├── test_data_import.py           # Tests for data_import module
├── test_statistical_analysis.py  # Tests for statistical_analysis module
├── test_normalization.py        # Tests for normalization module
├── test_preprocessing.py         # Tests for preprocessing module (partial)
└── README.md                     # This file
```

## Running Tests

### Install Test Dependencies

```bash
pip install pytest pytest-cov
```

### Run All Tests

```bash
# Run all tests
pytest

# Run with verbose output
pytest -v

# Run with coverage report
pytest --cov=proteomics_toolkit

# Run specific test file
pytest tests/test_basic.py

# Run specific test function
pytest tests/test_basic.py::test_basic_functionality
```

### Run Tests by Category

```bash
# Run only unit tests (when marked)
pytest -m unit

# Skip slow tests
pytest -m "not slow"
```

## Test Coverage

The test suite covers:

### ✅ Completed Test Modules

1. **Basic Functionality** (`test_basic.py`)
   - Import verification
   - Basic data operations
   - Configuration classes
   - Sample data structures

2. **Data Import** (`test_data_import.py`)
   - `load_skyline_data()` - Loading protein and metadata files
   - `parse_uniprot_identifier()` - UniProt ID parsing
   - `parse_gene_from_description()` - Gene name extraction
   - `clean_description()` - Description cleaning
   - `identify_sample_columns()` - Sample column identification
   - `clean_sample_names()` - Sample name cleaning
   - `match_samples_to_metadata()` - Sample-metadata matching
   - `identify_and_classify_controls()` - Control sample identification

3. **Statistical Analysis** (`test_statistical_analysis.py`)
   - `StatisticalConfig` - Configuration class
   - `prepare_metadata_dataframe()` - Metadata preparation
   - `run_paired_t_test()` - Paired t-test analysis
   - `run_unpaired_t_test()` - Unpaired t-test analysis
   - `apply_multiple_testing_correction()` - P-value correction
   - `run_comprehensive_statistical_analysis()` - Main analysis function
   - `display_analysis_summary()` - Results summary
   - Helper functions and edge cases

4. **Normalization** (`test_normalization.py`)
   - `get_normalization_characteristics()` - Method characteristics
   - `median_normalize()` - Median normalization
   - `vsn_normalize()` - Variance stabilizing normalization
   - `quantile_normalize()` - Quantile normalization
   - `log_transform()` - Log transformation
   - `mad_normalize()` - MAD normalization
   - `z_score_normalize()` - Z-score normalization
   - `rlr_normalize()` - Robust linear regression
   - `loess_normalize()` - LOESS normalization
   - `handle_negative_values()` - Negative value handling
   - `analyze_negative_values()` - Negative value analysis
   - Normalization statistics functions

### Partial Test Modules

5. **Preprocessing** (`test_preprocessing.py`)
   - **Note**: Some tests need function signature fixes
   - `parse_protein_identifiers()` - Protein ID parsing
   - `parse_gene_and_description()` - Gene/description parsing
   - `identify_annotation_columns()` - Column identification
   - Basic preprocessing pipeline tests

### Planned Test Modules

6. **Visualization** (not yet created)
   - Plot functions
   - Data visualization utilities
   - Chart generation and formatting

7. **Export** (not yet created) 
   - Data export functions
   - Configuration export
   - Result formatting and export

## Test Features

### Fixtures

The `conftest.py` provides shared fixtures:

- `sample_protein_data` - Mock protein quantification data
- `sample_metadata` - Mock sample metadata dictionary
- `sample_metadata_df` - Mock metadata DataFrame
- `statistical_config` - Configured StatisticalConfig object
- `sample_columns` - List of sample column names
- `temp_csv_files` - Temporary CSV files for I/O testing
- `differential_results` - Mock differential analysis results
- `protein_identifiers` - Sample protein IDs for parsing tests

### Test Categories

- **Unit Tests**: Test individual functions in isolation
- **Integration Tests**: Test workflows and function interactions
- **Edge Cases**: Test error conditions and boundary cases
- **Mock Data**: Use synthetic data that mimics real proteomics datasets

## Current Status

### Working Tests 
- Basic functionality tests pass
- Most core functions have test coverage
- Fixtures provide realistic test data
- Pytest configuration is properly set up

### Known Issues
- Some test functions need signature fixes to match actual module functions
- Mixed-effects model tests may need statsmodels dependency
- Some preprocessing tests need function signature updates

### Running Status
```bash
# Last test run results
pytest tests/test_basic.py -v
# Result: 7 passed in 0.18s ✅
```

## Development Guidelines

### Adding New Tests

1. **Create test file**: `test_module_name.py`
2. **Import modules**: Import functions to test
3. **Use fixtures**: Leverage existing fixtures from `conftest.py`
4. **Test structure**: Use classes to group related tests
5. **Descriptive names**: Use clear, descriptive test function names
6. **Documentation**: Add docstrings to test functions

### Test Naming Convention

```python
class TestFunctionName:
    """Test the function_name function"""
    
    def test_function_name_basic(self):
        """Test basic functionality"""
        pass
        
    def test_function_name_edge_case(self):
        """Test edge case handling"""
        pass
        
    def test_function_name_error_conditions(self):
        """Test error conditions"""
        pass
```

### Example Test

```python
def test_median_normalize_basic(self, sample_protein_data, sample_columns):
    """Test basic median normalization"""
    result = median_normalize(sample_protein_data, sample_columns)
    
    assert isinstance(result, pd.DataFrame)
    assert result.shape == sample_protein_data.shape
    
    # Check that medians are more similar after normalization
    sample_data = result[sample_columns]
    medians = sample_data.median()
    median_diff = medians.max() - medians.min()
    assert median_diff < 1.0
```

## Contributing

When adding new tests:

1. Check that function signatures match the actual module functions
2. Use appropriate fixtures from `conftest.py`
3. Test both success and failure cases
4. Add integration tests for workflows
5. Run tests before committing: `pytest -v`

## Debugging Tests

```bash
# Run with detailed output
pytest -v -s

# Stop on first failure
pytest -x

# Run specific test with debugging
pytest tests/test_module.py::TestClass::test_method -v -s

# Show local variables on failure
pytest --tb=long
```
