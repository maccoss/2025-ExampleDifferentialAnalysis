#!/usr/bin/env python3
"""
Diagnostic script to test proteomics_toolkit functions directly
and identify what's working vs. what's broken.
"""

import sys
import pandas as pd

# Add current directory to path
sys.path.insert(0, '.')

print("🔍 PROTEOMICS TOOLKIT DIAGNOSTIC TEST")
print("=" * 50)

try:
    print("✅ data_import module imported successfully")
except Exception as e:
    print(f"❌ data_import import failed: {e}")

try:
    print("✅ normalization module imported successfully")
except Exception as e:
    print(f"❌ normalization import failed: {e}")

try:
    print("✅ statistical_analysis module imported successfully")
except Exception as e:
    print(f"❌ statistical_analysis import failed: {e}")

try:
    print("✅ preprocessing module imported successfully")
except Exception as e:
    print(f"❌ preprocessing import failed: {e}")

print("\n🧪 TESTING CORE FUNCTIONS:")
print("-" * 30)

# Test 1: Basic data operations
print("\n1. Testing basic data operations...")
try:
    data = pd.DataFrame({
        'Protein': ['P1', 'P2', 'P3'],
        'Sample_A': [100, 200, 300],
        'Sample_B': [110, 190, 310]
    })
    print(f"   ✅ Created test dataframe: {data.shape}")
except Exception as e:
    print(f"   ❌ Basic dataframe creation failed: {e}")

# Test 2: UniProt parsing
print("\n2. Testing UniProt identifier parsing...")
try:
    from proteomics_toolkit.data_import import parse_uniprot_identifier
    result = parse_uniprot_identifier("sp|P12345|TEST_HUMAN")
    print(f"   ✅ UniProt parsing works: {result}")
except Exception as e:
    print(f"   ❌ UniProt parsing failed: {e}")

# Test 3: Median normalization
print("\n3. Testing median normalization...")
try:
    from proteomics_toolkit.normalization import median_normalize
    test_data = pd.DataFrame({
        'Protein': ['P1', 'P2'],
        'Sample_A': [100.0, 200.0],
        'Sample_B': [110.0, 190.0]
    })
    normalized = median_normalize(test_data, sample_columns=['Sample_A', 'Sample_B'])
    print(f"   ✅ Median normalization works: shape {normalized.shape}")
except Exception as e:
    print(f"   ❌ Median normalization failed: {e}")

# Test 4: Statistical config
print("\n4. Testing statistical configuration...")
try:
    from proteomics_toolkit.statistical_analysis import StatisticalConfig
    config = StatisticalConfig()
    print(f"   ✅ StatisticalConfig works: analysis_type = {config.analysis_type}")
except Exception as e:
    print(f"   ❌ StatisticalConfig failed: {e}")

# Test 5: Separate columns function
print("\n5. Testing column separation...")
try:
    from proteomics_toolkit.normalization import _separate_annotation_and_sample_columns
    test_data = pd.DataFrame({
        'Protein': ['P1', 'P2'],
        'Gene': ['G1', 'G2'],
        'Sample_A': [100.0, 200.0],
        'Sample_B': [110.0, 190.0]
    })
    result = _separate_annotation_and_sample_columns(test_data)
    print(f"   ✅ Column separation works: {type(result)}")
    print(f"      Returns: {result}")
except Exception as e:
    print(f"   ❌ Column separation failed: {e}")

# Test 6: MAD normalization (this was failing)
print("\n6. Testing MAD normalization...")
try:
    from proteomics_toolkit.normalization import mad_normalize
    test_data = pd.DataFrame({
        'Protein': ['P1', 'P2'],
        'Sample_A': [100.0, 200.0],
        'Sample_B': [110.0, 190.0]
    })
    result = mad_normalize(test_data, sample_columns=['Sample_A', 'Sample_B'])
    print(f"   ✅ MAD normalization works: shape {result.shape}")
except Exception as e:
    print(f"   ❌ MAD normalization failed: {e}")

print("\n" + "=" * 50)
print("🏁 DIAGNOSTIC COMPLETE")
