#!/usr/bin/env python3
"""
Regression test for VaxLab-report

This script runs evaluate_only.py and report_only.py with a known input
and checks that the output matches expected results.
"""

import os
import sys
import json
import tempfile
import shutil
import subprocess
from pathlib import Path

# Use test file included in repository
TEST_INPUT_FILE = "test_data.fasta"

def run_command(cmd, cwd=None):
    """Run a command and return the result"""
    try:
        # Check if we're in a conda environment or if conda is available
        conda_env = os.environ.get('CONDA_DEFAULT_ENV')
        if conda_env:
            # Already in a conda environment, run command directly
            full_cmd = cmd
        else:
            # Try to detect conda and activate environment
            conda_exe = shutil.which('conda')
            if conda_exe:
                # Use detected conda
                full_cmd = f'eval "$(conda shell.bash hook)" && conda activate base && {cmd}'
            else:
                # No conda found, run command directly (assume dependencies are available)
                full_cmd = cmd
        
        result = subprocess.run(
            full_cmd, 
            shell=True, 
            capture_output=True, 
            text=True, 
            cwd=cwd,
            timeout=120,  # 2 minute timeout
            executable='/bin/bash'  # Use bash explicitly
        )
        return result.returncode, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        return -1, "", "Command timed out"

def test_evaluate_only():
    """Test evaluate_only.py"""
    print("Testing evaluate_only.py...")
    
    # Check if test input file exists
    if not os.path.exists(TEST_INPUT_FILE):
        print(f"❌ Test input file not found: {TEST_INPUT_FILE}")
        return False
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Copy test input file
        input_file = os.path.join(tmpdir, "test_input.fasta")
        shutil.copy2(TEST_INPUT_FILE, input_file)
        
        # Run evaluate_only.py
        cmd = f"python vaxlab_report/evaluate_only.py -i {input_file} -o {tmpdir} --preset parameters.json"
        returncode, _, stderr = run_command(cmd)
        
        if returncode != 0:
            print(f"❌ evaluate_only.py failed with return code {returncode}")
            print(f"stderr: {stderr}")
            return False
        
        # Check if required output files exist
        required_files = [
            "evaluation_result_cds.json",
            "evaluation_result_mrna.json", 
            "checkpoints_cds.tsv",
            "checkpoints_mrna.tsv",
            "evaluate_only_log.txt"
        ]
        
        for filename in required_files:
            filepath = os.path.join(tmpdir, filename)
            if not os.path.exists(filepath):
                print(f"❌ Missing output file: {filename}")
                return False
        
        # Validate JSON structure
        for json_file in ["evaluation_result_cds.json", "evaluation_result_mrna.json"]:
            try:
                with open(os.path.join(tmpdir, json_file), 'r') as f:
                    data = json.load(f)
                    
                # Check required keys
                required_keys = ["input_sequence_id", "evaluation_status", "global_metrics", "local_metrics"]
                for key in required_keys:
                    if key not in data:
                        print(f"❌ Missing key '{key}' in {json_file}")
                        return False
                        
                # Check that we have some metrics
                if not data["global_metrics"]:
                    print(f"❌ No global metrics in {json_file}")
                    return False
                    
            except json.JSONDecodeError as e:
                print(f"❌ Invalid JSON in {json_file}: {e}")
                return False
        
        print("✅ evaluate_only.py test passed")
        return True

def test_report_only():
    """Test report_only.py"""
    print("Testing report_only.py...")
    
    # Check if test input file exists
    if not os.path.exists(TEST_INPUT_FILE):
        print(f"❌ Test input file not found: {TEST_INPUT_FILE}")
        return False
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Copy test input file
        input_file = os.path.join(tmpdir, "test_input.fasta")
        shutil.copy2(TEST_INPUT_FILE, input_file)
        
        # First run evaluate_only to create required JSON files
        cmd = f"python vaxlab_report/evaluate_only.py -i {input_file} -o {tmpdir} --preset parameters.json"
        returncode, _, stderr = run_command(cmd)
        
        if returncode != 0:
            print(f"❌ Setup failed: evaluate_only.py returned {returncode}")
            print(f"stderr: {stderr}")
            return False
        
        # Run report_only.py
        cmd = f"python vaxlab_report/report_only.py -i {input_file} -o {tmpdir}"
        returncode, _, stderr = run_command(cmd)
        
        if returncode != 0:
            print(f"❌ report_only.py failed with return code {returncode}")
            print(f"stderr: {stderr}")
            return False
        
        # Check if report.html was created
        report_file = os.path.join(tmpdir, "report.html")
        if not os.path.exists(report_file):
            print("❌ Missing output file: report.html")
            return False
        
        # Basic HTML validation
        with open(report_file, 'r') as f:
            html_content = f.read()
            
        required_elements = [
            "<html>",
            "</html>",
            "VaxLab mRNA Optimization Result",
            "Global Sequence Metrics",
            "Local Sequence Metrics"
        ]
        
        for element in required_elements:
            if element not in html_content:
                print(f"❌ Missing HTML element: {element}")
                return False
        
        print("✅ report_only.py test passed")
        return True

def test_with_forna():
    """Test report_only.py with --forna option"""
    print("Testing report_only.py with --forna option...")
    
    # Check if test input file exists
    if not os.path.exists(TEST_INPUT_FILE):
        print(f"❌ Test input file not found: {TEST_INPUT_FILE}")
        return False
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Copy test input file
        input_file = os.path.join(tmpdir, "test_input.fasta")
        shutil.copy2(TEST_INPUT_FILE, input_file)
        
        # First run evaluate_only
        cmd = f"python vaxlab_report/evaluate_only.py -i {input_file} -o {tmpdir} --preset parameters.json"
        returncode, _, stderr = run_command(cmd)
        
        if returncode != 0:
            print(f"❌ Setup failed: evaluate_only.py returned {returncode}")
            return False
        
        # Test both forna options
        for forna_option in ['qbio', 'tbi']:
            cmd = f"python vaxlab_report/report_only.py -i {input_file} -o {tmpdir} --forna {forna_option}"
            returncode, _, stderr = run_command(cmd)
            
            if returncode != 0:
                print(f"❌ report_only.py with --forna {forna_option} failed with return code {returncode}")
                print(f"stderr: {stderr}")
                return False
            
            # Check HTML contains forna URL
            report_file = os.path.join(tmpdir, "report.html")
            with open(report_file, 'r') as f:
                html_content = f.read()
            
            if forna_option == 'qbio':
                if "pub-forna.qbio.io" not in html_content:
                    print(f"❌ qbio Forna URL not found in HTML")
                    return False
            elif forna_option == 'tbi':
                if "nibiru.tbi.univie.ac.at" not in html_content:
                    print(f"❌ TBI Forna URL not found in HTML")
                    return False
        
        print("✅ Forna options test passed")
        return True

def main():
    """Run all regression tests"""
    print("🧪 Running VaxLab-report regression tests...\n")
    
    # Change to repository root
    repo_root = Path(__file__).parent
    os.chdir(repo_root)
    
    # Check if required files exist
    required_files = [
        "vaxlab_report/evaluate_only.py",
        "vaxlab_report/report_only.py", 
        "parameters.json"
    ]
    
    for filepath in required_files:
        if not os.path.exists(filepath):
            print(f"❌ Required file not found: {filepath}")
            sys.exit(1)
    
    # Run tests
    tests = [
        test_evaluate_only,
        test_report_only,
        test_with_forna
    ]
    
    passed = 0
    total = len(tests)
    
    for test_func in tests:
        try:
            if test_func():
                passed += 1
            else:
                print(f"❌ Test {test_func.__name__} failed")
        except Exception as e:
            print(f"❌ Test {test_func.__name__} crashed: {e}")
        print()  # Empty line between tests
    
    # Results
    print(f"📊 Test Results: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All regression tests passed!")
        print()
        print("🔍 Test Summary:")
        print("  ✅ evaluate_only.py - CDS and mRNA evaluation")
        print("  ✅ report_only.py - HTML report generation")
        print("  ✅ Forna options - Structure visualization (qbio/tbi)")
        print("  ✅ File outputs - All expected files validated")
        print("  ✅ JSON structure - Evaluation results properly formatted")
        print()
        print("📦 VaxLab-report is ready for production use!")
        sys.exit(0)
    else:
        print("💥 Some tests failed!")
        print("Please check the error messages above and fix any issues before deployment.")
        sys.exit(1)

if __name__ == "__main__":
    main()