import Mathlib

/-!
# Spectral gap of the discrete harmonic Hessian and the stiff Green bound

The eigenvalues of the discrete harmonic Hessian are `(2m/d) (cos θ - cos (k π / N))`. Lemma D
(correspondence 002) and the bound `M ≤ 3T/(8m)` on the stiff Green function (correspondence 011)
both rest on a lower bound for a difference of cosines obtained from Jordan's inequality, and the
second also on the value of `∑ 1/(k² - 1)`. This file proves those two ingredients.

The identification of the eigenvalues and the assembly of the matrix bound are not formalized here.
-/

namespace QuantumPathIntegrals

open Real

/-- Jordan-type lower bound: for `0 ≤ a ≤ b` with `a + b ≤ π`,
`cos a - cos b ≥ 2 (b² - a²) / π²`. -/
theorem cos_sub_cos_ge {a b : ℝ} (ha : 0 ≤ a) (hab : a ≤ b) (hsum : a + b ≤ π) :
    2 * (b ^ 2 - a ^ 2) / π ^ 2 ≤ cos a - cos b := by
  have hπ : 0 < π := pi_pos
  have hπ0 : π ≠ 0 := hπ.ne'
  have h1 : 2 / π * ((a + b) / 2) ≤ sin ((a + b) / 2) :=
    mul_le_sin (by linarith) (by linarith)
  have h2 : 2 / π * ((b - a) / 2) ≤ sin ((b - a) / 2) :=
    mul_le_sin (by linarith) (by linarith)
  have h3 : cos a - cos b = 2 * sin ((a + b) / 2) * sin ((b - a) / 2) := by
    rw [Real.cos_sub_cos]
    have hneg : sin ((a - b) / 2) = - sin ((b - a) / 2) := by
      rw [← sin_neg]; congr 1; ring
    rw [hneg]; ring
  have hl1 : 0 ≤ 2 / π * ((a + b) / 2) := mul_nonneg (by positivity) (by linarith)
  have hl2 : 0 ≤ 2 / π * ((b - a) / 2) := mul_nonneg (by positivity) (by linarith)
  have hs1 : 0 ≤ sin ((a + b) / 2) := le_trans hl1 h1
  have hprod := mul_le_mul h1 h2 hl2 hs1
  have heq : 2 * (b ^ 2 - a ^ 2) / π ^ 2
      = 2 * ((2 / π * ((a + b) / 2)) * (2 / π * ((b - a) / 2))) := by
    field_simp <;> ring
  rw [h3, heq]
  linarith

/-- The gap between the first and the `k`-th eigenvalue, in cosine form:
for `2 ≤ k` and `k + 1 ≤ N`, `cos (π/N) - cos (k π/N) ≥ 2 (k² - 1) / N²`. -/
theorem cos_gap_ge (N k : ℕ) (hk2 : 2 ≤ k) (hkN : k + 1 ≤ N) :
    2 * ((k : ℝ) ^ 2 - 1) / (N : ℝ) ^ 2 ≤ cos (π / N) - cos (k * π / N) := by
  have hNnat : 0 < N := by omega
  have hN : (0 : ℝ) < N := by exact_mod_cast hNnat
  have hk : (1 : ℝ) ≤ k := by exact_mod_cast (by omega : 1 ≤ k)
  have hkN' : (k : ℝ) + 1 ≤ N := by exact_mod_cast hkN
  have hπ : 0 < π := pi_pos
  have ha : 0 ≤ π / N := by positivity
  have hab : π / N ≤ k * π / N := by
    rw [div_le_div_iff_of_pos_right hN]
    nlinarith
  have hsum : π / N + k * π / N ≤ π := by
    rw [← add_div, div_le_iff₀ hN]
    nlinarith
  have key := cos_sub_cos_ge ha hab hsum
  have heq : 2 * ((k * π / N) ^ 2 - (π / N) ^ 2) / π ^ 2 = 2 * ((k : ℝ) ^ 2 - 1) / (N : ℝ) ^ 2 := by
    have hN0 : (N : ℝ) ≠ 0 := hN.ne'
    have hπ0 : π ≠ 0 := hπ.ne'
    field_simp <;> ring
  rw [heq] at key
  exact key

/-- Closed form of the partial sums of `∑ 1 / (k² - 1)` starting at `k = 2`. -/
theorem sum_inv_sq_sub_one_eq (n : ℕ) :
    ∑ j ∈ Finset.range n, (1 : ℝ) / (((j : ℝ) + 2) ^ 2 - 1)
      = 3 / 4 - (2 * (n : ℝ) + 3) / (2 * ((n : ℝ) + 1) * ((n : ℝ) + 2)) := by
  induction n with
  | zero => norm_num
  | succ n ih =>
    rw [Finset.sum_range_succ, ih]
    push_cast
    have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
    have h1 : ((n : ℝ) + 1) ≠ 0 := by positivity
    have h2 : ((n : ℝ) + 2) ≠ 0 := by positivity
    have h3 : ((n : ℝ) + 1 + 1) ≠ 0 := by positivity
    have h4 : ((n : ℝ) + 1 + 2) ≠ 0 := by positivity
    have h5 : ((n : ℝ) + 2) ^ 2 - 1 ≠ 0 := by
      have : (0 : ℝ) < ((n : ℝ) + 2) ^ 2 - 1 := by nlinarith
      exact this.ne'
    field_simp
    ring

/-- The partial sums are bounded by `3/4`, the constant in `M ≤ 3T/(8m)`. -/
theorem sum_inv_sq_sub_one_le (n : ℕ) :
    ∑ j ∈ Finset.range n, (1 : ℝ) / (((j : ℝ) + 2) ^ 2 - 1) ≤ 3 / 4 := by
  rw [sum_inv_sq_sub_one_eq]
  have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  have : (0 : ℝ) ≤ (2 * (n : ℝ) + 3) / (2 * ((n : ℝ) + 1) * ((n : ℝ) + 2)) := by positivity
  linarith

end QuantumPathIntegrals

#print axioms QuantumPathIntegrals.cos_sub_cos_ge
#print axioms QuantumPathIntegrals.cos_gap_ge
#print axioms QuantumPathIntegrals.sum_inv_sq_sub_one_eq
#print axioms QuantumPathIntegrals.sum_inv_sq_sub_one_le
