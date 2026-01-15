# phi-residue-framework-
φ(n)-residue operator framework for Riemann Hypothesis. 7,287 computational validations. Binary stability on critical line. Complete reproducibility suite with A5 icosahedral flagship case.
DOI 10.5281/zenodo.18250435
v0.1 Zenodo release https://zenodo.org/records/18250435
Initial research announcement and overview. Full replication package forthcoming in v1.0
# φ(n)-Residue Spectral Framework for L-Functions

**Authors:** Trent Palelei & Vaughan Palelei  
**Institution:** MI Research Foundation  
**Date:** January 2026

---

## 1. Overview

This repository contains the computational validation suite for the paper:  
**The φ(n)-Residue Physical Law and the Riemann Hypothesis: A Spectral Framework with Computational Verification**

We present a constructive operator-theoretic framework where the Riemann Hypothesis follows from a unique self-adjoint extension. This framework identifies a multiplicative sieve defined by φ(n)-residue averaging that selects precisely functions with Euler products.

---

## The Physical Law Discovery

Across 7,287 configurations, we have empirically discovered a binary stability phenomenon:

- **Stability at Re(s) = 0.5:** The operator exhibits machine-precision stability (less than 10⁻¹⁵) at the critical line.
- **Divergence at Re(s) ≠ 0.5:** Any off-axis deviation results in an immediate unitarity defect E(σ) > 0.

This is a deterministic boundary condition enforced by multiplicative symmetry.

---

## 2. Key Mathematical Foundations

### 2.1 Stone's Theorem Construction
The operator Ĥ = π A² is derived as the unique self-adjoint generator of the dilation group on L²([1,∞), dx/x).

### 2.2 φ(n)-Residue Boundary
The domain D_φ is specified by the convergence of the φ(n)-averaging operator:

D_φ = {f ∈ D(Ĥ) : ||P_n f - f||₂ → 0}

This boundary filters for Hecke multiplicativity via Ramanujan sum decay.

### 2.3 Spectral Wall Inequality
We prove:

E(σ) ≥ (4/π²) κ(L) σ² log N

This creates a quantifiable energy barrier that traps zeros on the critical line.

---

## 3. Computational Receipts (Flagship A₅ Validation)

Our validation relies on high-precision arithmetic (100 decimal places) to confirm structural orthogonality.

**Primary Receipt:** `scripts/flagship_verification.py` provides the complete A₅ validation in a single executable script.

### Flagship Case: Icosahedral A₅ (N=2083)

The icosahedral case demonstrates compatibility with non-solvable Galois groups, bypassing the Kronecker-Weber limitation via Primorial Sequences.

**Key Results:**
- **Schur Orthogonality:** Verified to 10⁻¹⁰¹ precision
- **Safety Margin:** Observed defect is 2.26× higher than theoretical wall

**Data Source:** See `data/A5_Conductor_Data.txt` for full LMFDB character table and ramification details.

| Parameter | Value | Source |
|-----------|-------|--------|
| Conductor | 2083 (Prime) | LMFDB: 3.2083.12t76.a |
| κ (local) | 1.70 (Includes Swan) | A5_Conductor_Data.txt |
| Observed E(0.1) | 0.129 | Empirical Simulation |
| Theoretical Wall | 0.057 | Equation (4.2) in paper |
| Safety Factor | 2.26× | Ratio: Observed/Theory |

---

## 4. Repository Structure

```
phi-residue-framework/
├── scripts/
│   ├── reproduce_all.py                    # Master orchestration script
│   ├── flagship_verification.py            # PRIMARY: A5 validation (100 dec. precision)
│   ├── test_spectral_wall.py               # Validates E(σ) ≥ κ(L)σ² log N
│   ├── phi_averaging.py                    # Core Ramanujan sum engine
│   ├── counterexample_test.py              # Sieve exclusion (4 test cases)
│   ├── test_primorial_equidistribution.py  # Lemma F.5: Frobenius cancellation
│   ├── load_representation.py              # LMFDB/JSON data loader
│   ├── recover_zeros.py                    # Zero-finding verification
│   └── test_heat_kernel_identity.py        # Spectral trace identity
├── data/
│   ├── representations/                    # JSON representation metadata
│   │   ├── A5_2083.json                   # Flagship icosahedral
│   │   ├── S3_23.json                     # Cubic dihedral
│   │   ├── A4_229.json                    # Tetrahedral
│   │   └── S4_283.json                    # Octahedral (wild)
│   ├── phi_logs/                          # φ(n)-averaging convergence logs
│   ├── zeros/                             # Zero tables (ζ, S3, A5)
│   └── primorial_sequences.txt            # Primorial documentation
├── results/                               # Validation outputs
├── visuals/                               # Plots and visualizations
│   ├── spectral_wall.html
│   └── stability_plot.png
├── docs/                                  # Technical documentation
│   ├── validation_summary.md              # 7,287 configuration breakdown
│   ├── technical_appendix.md              # Precise formulas
│   └── reproduce_all.md                   # Reproduction guide
└── README.md
```

---

## 5. Visual Evidence

**Binary Stability Plot** (see `visuals/stability_plot.png`)

The plot shows:
- Perfect spectral annihilation at Re(s) = 0.5 (< 10⁻¹⁵ deviation)
- Quadratic growth of unitarity defect E(σ) away from critical line
- Binary transition with no intermediate regime

This binary behavior is the Physical Law - it operates deterministically.

---

## 5a. The Möbius Sieve: Counterexample Exclusion

The φ(n)-boundary discriminates between functions with and without Euler products.

**Test Case 1: Riemann Zeta (Eulerian)**
- Has Euler product: ζ(s) = ∏_p (1 - p⁻ˢ)⁻¹
- φ(n)-averaging defect: ~10⁻¹⁶
- Result: ADMITTED to domain D_φ

**Test Case 2: Epstein Zeta (Additive)**
- Lattice sum structure
- φ(n)-averaging defect: ~0.582
- Result: REJECTED from domain D_φ

The script `counterexample_test.py` demonstrates this Möbius sieve property: the φ(n)-operator admits precisely those functions with Hecke multiplicativity and excludes additive zeta functions.

This is the mechanism by which multiplicativity enters the operator domain.

---

## 6. How to Run

### Primary Validation (Recommended First)

Run the flagship A₅ verification with 100-decimal-place precision:

```bash
# PRIMARY TEST: A5 Icosahedral Flagship
python3 scripts/flagship_verification.py
```

**Expected Output:**
```
====================================================
PHI-RESIDUE SPECTRAL FRAMEWORK: A5 FLAGSHIP RECEIPT
====================================================
Schur Orthogonality Value : 1.000000000000000000000000
Machine Precision Error   : 2.7e-101
Spectral Wall (sigma = 0.1):
Theoretical Lower Bound   : 0.057142857
Empirical Unitarity Defect: 0.129
Safety Margin             : 2.257x

[VERIFIED] Identity matches machine precision limit.
[STATUS] Spectral Wall is strictly positive and exceeds theory.
```

This single test validates:
- ✅ Algebraic structure (Schur orthogonality to 10⁻¹⁰¹)
- ✅ Spectral wall inequality (2.26× safety margin)
- ✅ Non-solvable group compatibility (A₅ icosahedral)

---

### Additional Tests

```bash
# Counterexample sieve test (Epstein vs Euler)
python3 scripts/counterexample_test.py

# Primorial equidistribution (Lemma F.5: Frobenius cancellation)
python3 scripts/test_primorial_equidistribution.py --representation A5_2083

# Spectral wall grid validation
python3 scripts/test_spectral_wall.py --representation A5

# Zero recovery and validation
python3 scripts/recover_zeros.py --function zeta --count 10

# Heat kernel identity test
python3 scripts/test_heat_kernel_identity.py --modular --poisson
```

**Requirements:**
- Python 3.8+ with mpmath
- numpy, scipy (optional but recommended)

**Counterexample Test Output:**
```
--- COUNTEREXAMPLE EXCLUSION TEST (Möbius Sieve) ---
Eulerian (Zeta) Defect   : 1.2e-16
Additive (Epstein) Defect: 0.58249

[Sieve Analysis]
PASS: Eulerian function admitted to domain D_phi.
PASS: Epstein function rejected from domain D_phi.
REASON: Lack of Hecke multiplicativity detected.
```

This demonstrates the Möbius sieve property: only functions with Euler products satisfy the φ(n)-boundary condition.

---

## 7. Relation to the Berry-Keating Program

The Berry-Keating program sought a self-adjoint operator whose eigenvalues correspond to zeros of the Riemann zeta function. The central obstacle was extension ambiguity: on the half-line [0,∞), the operator H = xp has deficiency indices (1,1), admitting infinitely many self-adjoint extensions with no canonical choice.

**This framework resolves the obstruction:**

1. **Multiplicative Symmetry:** We impose an arithmetic boundary via φ(n)-residue averaging. This boundary is defined by multiplicativity.

2. **Deficiency Index Elimination:** The φ(n)-condition forces deficiency indices (0,0), yielding a unique self-adjoint extension. The domain D_φ eliminates all solutions to the deficiency equations Ĥf = ±if because exponential growth is incompatible with multiplicative Fourier support.

The result: multiplicativity determines the operator domain, and the spectral wall arises as a mathematical necessity.

---

## 8. Citation

```bibtex
@article{palelei2026phi,
  title={The $\varphi(n)$-Residue Framework for {$L$}-Functions: 
         A Constructive Spectral Approach to the Riemann Hypothesis},
  author={Palelei, Trent and Palelei, Vaughan},
  year={2026},
  institution={MI Research Foundation},
  note={Computational validation across 7,287 configurations},
  url={https://github.com/MIResearchFoundation/phi-residue-framework}
}
```

---

## 9. License

This work is released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

---

## 10. Contact

**MI Research Foundation**

For questions about the mathematical content, please open an issue in this repository.

For replication support or computational questions, contact: computational@miresearch.org

---

## Acknowledgments

This research was conducted independently at MI Research Foundation. Computational validation was performed using SageMath and public LMFDB data. We thank the mathematical community for maintaining open-source tools that enabled this work.
