import LRC.Definitions
import Mathlib.Data.Int.Lemmas
import Mathlib.Analysis.SpecialFunctions.Integrals

/-!
# Common lemmas for Lonely Runner proofs

Shared infrastructure used across the case-specific proofs.
-/

noncomputable section

open Real

lemma circDist_nonneg (x : ℝ) : circDist x ≥ 0 := by
  unfold circDist
  simp [min_def]
  split
  · exact Int.fract_nonneg x
  · linarith [Int.fract_lt_one x]

lemma circDist_le_half (x : ℝ) : circDist x ≤ 1 / 2 := by
  unfold circDist
  simp [min_def]
  split
  · linarith [Int.fract_lt_one x]
  · linarith [Int.fract_nonneg x]

lemma circDist_int (n : ℤ) : circDist (n : ℝ) = 0 := by
  unfold circDist
  simp [Int.fract_intCast]

end
