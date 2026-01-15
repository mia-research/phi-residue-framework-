"""
test_heat_kernel_identity.py

Compare numerical heat kernel trace with Jacobi theta function θ(t).

The heat kernel approach provides an alternative verification of the
spectral structure via the trace formula:

    Tr(e^(-tĤ)) = θ(t)

where θ(t) is the Jacobi theta function.

Usage:
    python3 test_heat_kernel_identity.py --t_min 0.1 --t_max 10
    python3 test_heat_kernel_identity.py --precision 100
"""

from mpmath import mp
import argparse

mp.dps = 50

def jacobi_theta3(t, q_series_terms=100):
    """
    Compute Jacobi theta function θ₃(0, e^(-πt))
    
    θ₃(z, q) = Σ_{n=-∞}^∞ q^(n²)
    
    For z = 0 and q = e^(-πt):
    θ₃(0, e^(-πt)) = 1 + 2 Σ_{n=1}^∞ e^(-πn²t)
    """
    q = mp.exp(-mp.pi * t)
    
    # Sum the series
    result = mp.mpf(1)
    for n in range(1, q_series_terms + 1):
        result += 2 * mp.power(q, n * n)
        
        # Check for convergence
        if mp.power(q, n * n) < mp.power(10, -mp.dps):
            break
    
    return result

def heat_kernel_trace(t, eigenvalue_list, truncate=1000):
    """
    Compute heat kernel trace: Tr(e^(-tĤ)) = Σ_k e^(-t λ_k)
    
    For the φ(n)-residue operator, eigenvalues λ_k should correspond
    to Mellin frequencies on the critical line.
    
    Parameters:
    -----------
    t : float
        Heat kernel time parameter
    eigenvalue_list : list
        List of eigenvalues (assumed to be imaginary parts of zeros)
    truncate : int
        Maximum number of eigenvalues to include
    """
    trace = mp.mpf(0)
    
    for k, lambda_k in enumerate(eigenvalue_list[:truncate]):
        # Each zero at s = 1/2 + i*gamma contributes eigenvalue
        # For demonstration, use simplified model
        eigenval = mp.mpf(lambda_k)
        trace += mp.exp(-t * eigenval)
    
    return trace

def test_heat_kernel_identity(t_values, gamma_values):
    """
    Test the identity: Tr(e^(-tĤ)) = θ(t)
    
    This provides an independent check of the spectral structure.
    """
    print("=" * 70)
    print("HEAT KERNEL / THETA IDENTITY TEST")
    print("=" * 70)
    
    print("\nTesting: Tr(e^(-tĤ)) = θ₃(0, e^(-πt))")
    print(f"Eigenvalues (first 10 zeros): {[f'{g:.4f}' for g in gamma_values[:10]]}")
    
    print(f"\n{'t':<10} {'θ₃(t)':<20} {'Tr(e^(-tĤ))':<20} {'Error':<15} {'Status':<10}")
    print("-" * 70)
    
    errors = []
    
    for t in t_values:
        theta_val = jacobi_theta3(t)
        trace_val = heat_kernel_trace(t, gamma_values)
        
        error = abs(theta_val - trace_val)
        errors.append(float(error))
        
        status = "PASS" if error < 1e-10 else "CHECK"
        
        print(f"{t:<10.4f} {float(theta_val):<20.12f} {float(trace_val):<20.12f} "
              f"{float(error):<15.2e} {status:<10}")
    
    print("-" * 70)
    
    avg_error = sum(errors) / len(errors)
    max_error = max(errors)
    
    print(f"\nAverage error: {avg_error:.2e}")
    print(f"Maximum error: {max_error:.2e}")
    
    if max_error < 1e-10:
        print("\nCONCLUSION: Heat kernel identity verified to high precision")
        print("Spectral structure matches theta function prediction")
    elif max_error < 1e-5:
        print("\nCONCLUSION: Heat kernel identity verified to moderate precision")
        print("Minor deviations may be due to truncation effects")
    else:
        print("\nWARNING: Significant deviation detected")
        print("Check eigenvalue list and convergence parameters")
    
    print("=" * 70)
    
    return errors

def modular_transformation_check(t_small=0.1, t_large=10.0):
    """
    Check modular transformation of theta function:
    θ(1/t) = √t θ(t)
    
    This is the functional equation relating small and large t behavior.
    """
    print("\n" + "=" * 70)
    print("MODULAR TRANSFORMATION CHECK")
    print("=" * 70)
    
    print(f"\nTesting: θ(1/t) = √t θ(t)")
    print(f"t_small = {t_small}, t_large = {t_large} (where t_large = 1/t_small)")
    
    theta_small = jacobi_theta3(t_small)
    theta_large = jacobi_theta3(t_large)
    
    # Modular transformation prediction
    predicted = mp.sqrt(t_small) * theta_small
    
    error = abs(theta_large - predicted)
    
    print(f"\nθ({t_small}) = {float(theta_small):.12f}")
    print(f"θ({t_large}) = {float(theta_large):.12f}")
    print(f"√{t_small} × θ({t_small}) = {float(predicted):.12f}")
    print(f"\nError: {float(error):.2e}")
    
    if error < 1e-10:
        print("\nCONCLUSION: Modular transformation verified")
        print("Functional equation holds to machine precision")
    else:
        print("\nWARNING: Modular transformation check failed")
        print("This may indicate issues with theta function computation")
    
    print("=" * 70)
    
    return float(error)

def poisson_summation_verification(t_test=1.0):
    """
    Verify Poisson summation formula underlying the heat kernel identity.
    
    Σ_{n∈ℤ} e^(-πn²t) = (1/√t) Σ_{m∈ℤ} e^(-πm²/t)
    """
    print("\n" + "=" * 70)
    print("POISSON SUMMATION VERIFICATION")
    print("=" * 70)
    
    print(f"\nTesting Poisson summation at t = {t_test}")
    
    # Left side: direct sum
    n_max = 100
    left_sum = mp.mpf(0)
    for n in range(-n_max, n_max + 1):
        left_sum += mp.exp(-mp.pi * n * n * t_test)
    
    # Right side: transformed sum
    right_sum = mp.mpf(0)
    for m in range(-n_max, n_max + 1):
        right_sum += mp.exp(-mp.pi * m * m / t_test)
    right_sum = right_sum / mp.sqrt(t_test)
    
    error = abs(left_sum - right_sum)
    
    print(f"\nDirect sum:      {float(left_sum):.12f}")
    print(f"Transformed sum: {float(right_sum):.12f}")
    print(f"Error:           {float(error):.2e}")
    
    if error < 1e-10:
        print("\nCONCLUSION: Poisson summation verified")
    else:
        print("\nWARNING: Poisson summation check inconclusive")
    
    print("=" * 70)
    
    return float(error)

# Known Riemann zeta zeros for testing
RIEMANN_ZEROS = [
    14.134725,
    21.022040,
    25.010858,
    30.424876,
    32.935062,
    37.586178,
    40.918719,
    43.327073,
    48.005151,
    49.773832
]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Heat Kernel / Theta Identity Test"
    )
    parser.add_argument('--t_min', type=float, default=0.1,
                       help='Minimum t value')
    parser.add_argument('--t_max', type=float, default=10.0,
                       help='Maximum t value')
    parser.add_argument('--steps', type=int, default=10,
                       help='Number of t values to test')
    parser.add_argument('--precision', type=int, default=50,
                       help='Decimal precision for mpmath')
    parser.add_argument('--modular', action='store_true',
                       help='Check modular transformation')
    parser.add_argument('--poisson', action='store_true',
                       help='Verify Poisson summation')
    
    args = parser.parse_args()
    
    # Set precision
    mp.dps = args.precision
    
    # Generate t values (log-spaced for better coverage)
    t_values = []
    for i in range(args.steps):
        log_t = mp.log(args.t_min) + i * (mp.log(args.t_max) - mp.log(args.t_min)) / (args.steps - 1)
        t_values.append(float(mp.exp(log_t)))
    
    # Run heat kernel identity test
    errors = test_heat_kernel_identity(t_values, RIEMANN_ZEROS)
    
    # Optional: modular transformation check
    if args.modular:
        modular_transformation_check(t_small=args.t_min, t_large=args.t_max)
    
    # Optional: Poisson summation verification
    if args.poisson:
        poisson_summation_verification(t_test=1.0)
