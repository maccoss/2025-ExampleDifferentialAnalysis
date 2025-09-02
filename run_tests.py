#!/usr/bin/env python3
"""
Test runner script for proteomics_toolkit

This script runs the test suite and provides a summary of results.
"""

import subprocess
import sys
import os


def run_command(cmd, description):
    """Run a command and return success status"""
    print(f"\n{'='*60}")
    print(f"🧪 {description}")
    print('='*60)
    
    try:
        result = subprocess.run(cmd, shell=True, check=False, cwd=os.path.dirname(__file__))
        if result.returncode == 0:
            print(f"✅ {description} - PASSED")
            return True
        else:
            print(f"❌ {description} - FAILED (exit code: {result.returncode})")
            return False
    except Exception as e:
        print(f"💥 {description} - ERROR: {e}")
        return False


def main():
    """Run the full test suite"""
    
    print("Proteomics Toolkit Test Suite")
    print("="*60)
    
    # Change to the project root directory
    project_root = os.path.dirname(os.path.abspath(__file__))
    os.chdir(project_root)
    
    # List of test commands to run
    test_commands = [
        ("python -m pytest tests/test_basic.py -v", "Basic Functionality Tests"),
        ("python -m pytest tests/test_data_import.py -v --tb=short", "Data Import Tests"),
        ("python -m pytest tests/test_statistical_analysis.py -v --tb=short", "Statistical Analysis Tests"), 
        ("python -m pytest tests/test_normalization.py -v --tb=short", "Normalization Tests"),
        ("python -m pytest tests/ --tb=short -q", "Complete Test Suite (Quick)"),
    ]
    
    results = []
    
    # Run each test command
    for cmd, description in test_commands:
        success = run_command(cmd, description)
        results.append((description, success))
    
    # Summary
    print(f"\n{'='*60}")
    print("📊 TEST SUMMARY")
    print('='*60)
    
    passed_count = 0
    total_count = len(results)
    
    for description, success in results:
        status = "✅ PASSED" if success else "❌ FAILED"
        print(f"{status:12} - {description}")
        if success:
            passed_count += 1
    
    print(f"\n Overall: {passed_count}/{total_count} test suites passed")
    
    if passed_count == total_count:
        print("All test suites completed successfully!")
        return 0
    else:
        print(f" {total_count - passed_count} test suite(s) had failures")
        return 1


if __name__ == "__main__":
    sys.exit(main())
