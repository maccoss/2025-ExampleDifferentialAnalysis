#!/usr/bin/env python3
"""
Simple test script for create_standard_data_structure function.
This avoids pytest and statsmodels import issues.
"""

import pandas as pd
import sys

# Add current directory to path
sys.path.append('.')

from proteomics_toolkit.preprocessing import create_standard_data_structure, parse_protein_identifiers, parse_gene_and_description

def test_basic_functionality():
    """Test basic functionality of create_standard_data_structure."""
    print("=== Testing create_standard_data_structure ===\n")
    
    # Create sample data that mimics the CSV input format
    sample_data = pd.DataFrame({
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
    
    print("1. Testing preprocessing steps...")
    
    # Process the data first (like in the notebook)
    processed_data = parse_protein_identifiers(sample_data)
    processed_data = parse_gene_and_description(processed_data)
    
    print(f"   Processed data shape: {processed_data.shape}")
    print(f"   Processed columns: {list(processed_data.columns)}")
    
    print("\n2. Testing standard data structure creation...")
    
    # Define cleaned sample names  
    cleaned_sample_names = {
        'Total-PTE01-511-84A-C4-049 Sum Normalized Area': 'E01-511-84A-C4-049',
        'Total-PTE02-Hoof18-050 Sum Normalized Area': 'E02-Hoof18-050',  
        'Total-PTF01-041-40B-B1-061 Sum Normalized Area': 'F01-041-40B-B1-061'
    }
    
    # Test the function
    try:
        result = create_standard_data_structure(processed_data, cleaned_sample_names)
        print(f"   ✅ SUCCESS! Created standard structure with shape: {result.shape}")
        
        # Verify structure
        expected_annotation_cols = [
            'Protein', 'Description', 'Protein Gene', 'UniProt_Accession', 'UniProt_Entry_Name'
        ]
        
        actual_annotation_cols = list(result.columns[:5])
        if actual_annotation_cols == expected_annotation_cols:
            print(f"   ✅ Annotation columns correct: {actual_annotation_cols}")
        else:
            print(f"   ❌ Annotation columns wrong. Expected: {expected_annotation_cols}, Got: {actual_annotation_cols}")
            return False
            
        # Check sample columns
        sample_columns = list(result.columns[5:])
        expected_sample_cols = ['E01-511-84A-C4-049', 'E02-Hoof18-050', 'F01-041-40B-B1-061']
        
        if sample_columns == expected_sample_cols:
            print(f"   ✅ Sample columns correct: {sample_columns}")
        else:
            print(f"   ❌ Sample columns wrong. Expected: {expected_sample_cols}, Got: {sample_columns}")
            return False
            
        print("\n3. Testing annotation column content...")
        
        # Check UniProt_Accession parsing
        accessions = list(result['UniProt_Accession'])
        expected_accessions = ['P02787', 'O43240', 'P06858']
        if accessions == expected_accessions:
            print(f"   ✅ UniProt_Accession parsed correctly: {accessions}")
        else:
            print(f"   ❌ UniProt_Accession parsing failed. Expected: {expected_accessions}, Got: {accessions}")
            return False
            
        # Check Description cleaning
        first_description = result['Description'].iloc[0]
        if 'Serotransferrin' in first_description and 'OS=' not in first_description:
            print(f"   ✅ Description cleaned correctly: '{first_description}'")
        else:
            print(f"   ❌ Description cleaning failed: '{first_description}'")
            return False
            
        print("\n4. Testing data preservation...")
        
        # Check that sample values are preserved
        original_value = sample_data['Total-PTE01-511-84A-C4-049 Sum Normalized Area'].iloc[0]
        result_value = result['E01-511-84A-C4-049'].iloc[0]
        
        if result_value == original_value:
            print(f"   ✅ Sample data preserved: {result_value}")
        else:
            print(f"   ❌ Sample data not preserved. Expected: {original_value}, Got: {result_value}")
            return False
            
        print("\n🎉 ALL TESTS PASSED! The function works correctly.")
        return True
        
    except Exception as e:
        print(f"   ❌ ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_error_handling():
    """Test error handling for missing columns."""
    print("\n=== Testing Error Handling ===\n")
    
    # Create data without proper preprocessing
    raw_data = pd.DataFrame({
        'Protein': ['sp|P02787|TRFE_HUMAN'],
        'Protein Description': ['Serotransferrin OS=Homo sapiens OX=9606 GN=TF PE=1 SV=5'],
        'Protein Gene': ['TF'],
        'Sample1': [1000]
    })
    
    try:
        create_standard_data_structure(raw_data)
        print("   ❌ ERROR: Should have raised ValueError for missing columns")
        return False
    except ValueError as e:
        if "Missing required annotation columns" in str(e):
            print(f"   ✅ Correctly caught missing columns error: {e}")
            return True
        else:
            print(f"   ❌ Wrong error message: {e}")
            return False
    except Exception as e:
        print(f"   ❌ Unexpected error: {e}")
        return False

if __name__ == "__main__":
    print("Running standard data structure tests...\n")
    
    success1 = test_basic_functionality()
    success2 = test_error_handling()
    
    if success1 and success2:
        print("\n🎉 ALL TESTS PASSED!")
        sys.exit(0)
    else:
        print("\n❌ SOME TESTS FAILED!")
        sys.exit(1)
