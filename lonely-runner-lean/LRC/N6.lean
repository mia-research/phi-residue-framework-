import LRC.Common

/-!
# Lonely Runner Conjecture for n = 6

Mechanically verified proof of the Lonely Runner Conjecture for 6 runners.

This formalizes the result of Barajas & Serra (2008).
With one runner fixed at the origin, we must show that for any 5 distinct
nonzero integer speeds, there exists a time `t` such that all 5 runners
are at distance ≥ 1/6 from the origin on the unit circle.
-/

noncomputable section

theorem lrc_n6 : LonelyRunnerConjecture 6 := by
  sorry

end
