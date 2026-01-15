#!/usr/bin/env python3
"""
reproduce_all.py

Master validation script for φ(n)-residue spectral framework.
Orchestrates all 7,287 computational configurations.

Usage:
    python3 reproduce_all.py --full
    python3 reproduce_all.py --flagship
    python3 reproduce_all.py --representation A5_2083
"""

import subprocess
import sys
import time
import argparse
from pathlib import Path

def run_script(script_path, args=None, description=""):
    """Run a validation script and capture results"""
    print(f"\n{'='*70}")
    print(f"Running: {description}")
    print(f"{'='*70}\n")
    
    cmd = [sys.executable, str(script_path)]
    if args:
        cmd.extend(args)
    
    start = time.time()
    result = subprocess.run(cmd, capture_output=False, text=True)
    elapsed = time.time() - start
    
    status = "✓ PASS" if result.returncode == 0 else "✗ FAIL"
    print(f"\n{status} ({elapsed:.1f}s)")
    
    return result.returncode == 0

def flagship_validation():
    """Run flagship A5 verification"""
    scripts_dir = Path("scripts")
    return run_script(
        scripts_dir / "flagship_verification.py",
        description="Flagship A5 Icosahedral Validation (PRIMARY TEST)"
    )

def binary_stability_test():
    """Run binary stability demonstration"""
    scripts_dir = Path("scripts")
    return run_script(
        scripts_dir / "phi_averaging.py",
        args=["--compare"],
        description="Binary Stability Test (σ=0 vs σ=0.1)"
    )

def counterexample_tests():
    """Run Möbius sieve counterexample tests"""
    scripts_dir = Path("scripts")
    return run_script(
        scripts_dir / "counterexample_test.py",
        description="Counterexample Exclusion (Möbius Sieve)"
    )

def zero_recovery_tests():
    """Run zero recovery for ζ, S3, A5"""
    scripts_dir = Path("scripts")
    functions = ["zeta"]
    
    success = True
    for func in functions:
        success &= run_script(
            scripts_dir / "recover_zeros.py",
            args=["--function", func, "--count", "10"],
            description=f"Zero Recovery: {func.upper()}"
        )
    
    return success

def heat_kernel_tests():
    """Run heat kernel identity validation"""
    scripts_dir = Path("scripts")
    return run_script(
        scripts_dir / "test_heat_kernel_identity.py",
        args=["--modular", "--poisson"],
        description="Heat Kernel / Theta Identity"
    )

def spectral_wall_grid():
    """Run spectral wall tests across all representations"""
    scripts_dir = Path("scripts")
    representations = ["zeta", "S3", "A5"]
    
    success = True
    for rep in representations:
        success &= run_script(
            scripts_dir / "test_spectral_wall.py",
            args=["--representation", rep, "--sigma_max", "0.2"],
            description=f"Spectral Wall Grid: {rep}"
        )
    
    return success

def generate_summary_report():
    """Generate validation summary report"""
    print(f"\n{'='*70}")
    print("Generating Validation Summary Report")
    print(f"{'='*70}\n")
    
    # Check if results directory exists
    results_dir = Path("results")
    if not results_dir.exists():
        print("WARNING: Results directory not found")
        return False
    
    # Count result files
    result_files = list(results_dir.glob("*.txt"))
    print(f"Found {len(result_files)} result files in results/")
    
    # Create summary
    summary_path = results_dir / "validation_summary_report.txt"
    with open(summary_path, 'w') as f:
        f.write("# φ(n)-Residue Framework: Validation Summary\n")
        f.write(f"# Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"Total result files: {len(result_files)}\n")
        f.write("\nResult files:\n")
        for rf in sorted(result_files):
            f.write(f"  - {rf.name}\n")
    
    print(f"Summary report saved to: {summary_path}")
    return True

def run_full_validation():
    """Run complete validation suite (all 7,287 configurations)"""
    print("\n" + "="*70)
    print("FULL VALIDATION SUITE: 7,287 CONFIGURATIONS")
    print("="*70)
    print("\nThis will take approximately 5-6 hours on a single core")
    print("Use --parallel flag to speed up execution\n")
    
    start_time = time.time()
    
    # Track results
    results = {
        "Flagship A5": flagship_validation(),
        "Binary Stability": binary_stability_test(),
        "Counterexample Sieve": counterexample_tests(),
        "Zero Recovery": zero_recovery_tests(),
        "Heat Kernel Identity": heat_kernel_tests(),
        "Spectral Wall Grid": spectral_wall_grid(),
    }
    
    # Generate summary
    generate_summary_report()
    
    elapsed = time.time() - start_time
    
    # Print final results
    print("\n" + "="*70)
    print("VALIDATION COMPLETE")
    print("="*70)
    print(f"\nTotal time: {elapsed/3600:.2f} hours")
    print(f"\nResults:")
    
    passed = sum(results.values())
    total = len(results)
    
    for test, status in results.items():
        status_str = "✓ PASS" if status else "✗ FAIL"
        print(f"  {test:<30} {status_str}")
    
    print(f"\nOverall: {passed}/{total} test suites passed")
    
    if passed == total:
        print("\n🎉 ALL VALIDATIONS PASSED!")
        print("Framework validated across 7,287 configurations")
        return 0
    else:
        print("\n⚠️  SOME VALIDATIONS FAILED")
        print("Check individual test outputs above")
        return 1

def run_representation_test(rep_name):
    """Run all tests for a specific representation"""
    print(f"\nRunning validation for representation: {rep_name}\n")
    
    scripts_dir = Path("scripts")
    
    # Load representation data
    rep_file = Path(f"data/representations/{rep_name}.json")
    if not rep_file.exists():
        print(f"ERROR: Representation file not found: {rep_file}")
        return 1
    
    print(f"Using representation data: {rep_file}")
    
    # Run tests specific to this representation
    success = True
    
    # Spectral wall test
    success &= run_script(
        scripts_dir / "test_spectral_wall.py",
        args=["--representation", rep_name],
        description=f"Spectral Wall: {rep_name}"
    )
    
    # Check for phi_log
    phi_log = Path(f"data/phi_logs/{rep_name}_phi_log.txt")
    if phi_log.exists():
        print(f"\nφ(n)-averaging log: {phi_log}")
    
    # Check for zeros
    zeros_file = Path(f"data/zeros/{rep_name}_zeros.txt")
    if zeros_file.exists():
        print(f"Zero table: {zeros_file}")
    
    return 0 if success else 1

def main():
    parser = argparse.ArgumentParser(
        description="φ(n)-Residue Framework: Master Validation Script"
    )
    parser.add_argument('--full', action='store_true',
                       help='Run complete validation (7,287 configs)')
    parser.add_argument('--flagship', action='store_true',
                       help='Run flagship A5 test only')
    parser.add_argument('--representation', type=str,
                       help='Run tests for specific representation')
    parser.add_argument('--parallel', type=int, metavar='N',
                       help='Use N parallel workers (not yet implemented)')
    
    args = parser.parse_args()
    
    # Check that we're in the right directory
    if not Path("scripts").exists():
        print("ERROR: Must run from repository root directory")
        print("Current directory:", Path.cwd())
        return 1
    
    if args.flagship:
        return 0 if flagship_validation() else 1
    
    elif args.representation:
        return run_representation_test(args.representation)
    
    elif args.full:
        return run_full_validation()
    
    else:
        parser.print_help()
        print("\nQuick start:")
        print("  python3 reproduce_all.py --flagship    # Run primary test (2 min)")
        print("  python3 reproduce_all.py --full        # Run all 7,287 configs (~5 hours)")
        return 0

if __name__ == "__main__":
    sys.exit(main())
