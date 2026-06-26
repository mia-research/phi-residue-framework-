import Mathlib.Topology.MetricSpace.Basic
import Mathlib.Topology.Order
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic

/-!
# Lonely Runner Conjecture — Definitions

We formalize the Lonely Runner Conjecture for runners on a circular track
of circumference 1.

## Setup

Given `n` runners on a unit circle, each with distinct constant speeds
`v₁, v₂, …, vₙ`, all starting at the origin, the conjecture states that
for each runner `i`, there exists a time `t` such that runner `i` is at
distance at least `1/n` from every other runner.

We use the fractional-part formulation: the position of runner `i` at
time `t` is `{v_i · t}` (fractional part), and "distance on the circle"
is `‖{v_i · t} - {v_j · t}‖` where `‖x‖ = min(x mod 1, 1 - x mod 1)`.
-/

noncomputable section

open Real

/-- Distance on the unit circle `ℝ/ℤ`, i.e. `min({x}, 1 - {x})` where
    `{x}` denotes the fractional part. -/
def circDist (x : ℝ) : ℝ :=
  min (Int.fract x) (1 - Int.fract x)

/-- A runner's position at time `t` given speed `v`, on the unit circle. -/
def runnerPos (v t : ℝ) : ℝ :=
  Int.fract (v * t)

/-- The separation between two runners with speeds `v` and `w` at time `t`. -/
def runnerSep (v w t : ℝ) : ℝ :=
  circDist ((v - w) * t)

/-- The Lonely Runner Conjecture for `n` runners: given `n` distinct speeds,
    for each runner `i` there exists a time `t` such that runner `i` is at
    distance at least `1/n` from all other runners.

    We use the standard reduction: fix one runner at speed 0 and consider
    `n - 1` other runners with distinct nonzero speeds. Then the conjecture
    states that there exists `t` such that all `n - 1` runners are at
    distance at least `1/n` from the origin. -/
def LonelyRunnerConjecture (n : ℕ) : Prop :=
  ∀ (speeds : Fin (n - 1) → ℤ),
    Function.Injective speeds →
    (∀ i, speeds i ≠ 0) →
    ∃ t : ℝ, ∀ i : Fin (n - 1), circDist (speeds i * t) ≥ 1 / n

end
