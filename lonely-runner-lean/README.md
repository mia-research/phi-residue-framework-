# Lonely Runner Conjecture — Lean 4 Formalizations

**Status:** Mechanically verified proofs for n = 4, 6. n = 7 in progress.

## What this is

This repository contains Lean 4 formalizations of known analytical proofs
of the Lonely Runner Conjecture for small values of n. These are **not
original proofs** — they are mechanical verifications of results by:

- n = 4: Betke & Wills (1972), Cusick (1982)
- n = 6: Barajas & Serra (2008)

## The conjecture

Consider n runners on a circular track of circumference 1, all starting at
the same point and running at distinct constant speeds. The Lonely Runner
Conjecture asserts that each runner is, at some time, at distance at least
1/n from all other runners.

## Why formalize known results?

- Confirm no gaps or implicit assumptions were missed in the original proofs
- Provide didactic examples of analysis formalization in Lean
- Build infrastructure toward verifying n ≥ 7 or related conjectures

## Structure

```
LRC/
├── Definitions.lean   # The conjecture statement and key definitions
├── Common.lean        # Shared lemmas and tactics
├── N4.lean            # Proof for n = 4
├── N6.lean            # Proof for n = 6
└── N7.lean            # Work in progress
```

## Build

```bash
lake exe cache get
lake build
```

Requires Lean 4 and a compatible mathlib. See `lean-toolchain` for the
pinned version.

## References

- W. Betke, J.M. Wills. "Untere Schranken für zwei diophantische
  Approximations-Funktionen." Monatsh. Math. 76 (1972), 214–217.
- T.W. Cusick. "View-obstruction problems in n-dimensional geometry."
  J. Combin. Theory Ser. A 16 (1974), 1–11.
- J. Barajas, O. Serra. "The lonely runner with seven runners."
  Electron. J. Combin. 15 (2008), #R48.

## License

Apache-2.0. See [LICENSE](LICENSE).
