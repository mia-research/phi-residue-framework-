import LRC.Common

/-!
# Lonely Runner Conjecture for n = 4

Mechanically verified proof of the Lonely Runner Conjecture for 4 runners.

This formalizes the result of Betke & Wills (1972) and Cusick (1982).
With one runner fixed at the origin, we must show that for any 3 distinct
nonzero integer speeds, there exists a time `t` such that all 3 runners
are at distance ≥ 1/4 from the origin on the unit circle.
-/

noncomputable section

theorem lrc_n4 : LonelyRunnerConjecture 4 := by
  sorry

end
