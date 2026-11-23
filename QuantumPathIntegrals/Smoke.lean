import Mathlib

/-!
# Smoke test

These statements confirm that the pinned toolchain and Mathlib build locally
and in CI. They are not project results.
-/

namespace QuantumPathIntegrals

theorem smoke_add_comm (a b : ℝ) : a + b = b + a := add_comm a b

theorem smoke_exp_pos (x : ℝ) : 0 < Real.exp x := Real.exp_pos x

end QuantumPathIntegrals

#print axioms QuantumPathIntegrals.smoke_add_comm
#print axioms QuantumPathIntegrals.smoke_exp_pos
