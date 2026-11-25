import Mathlib

/-!
# Algebra behind the tilted-plane homotopy and the endpoint trace bound

Nova 023 (tilted-plane homotopy). On the plane `x = (1 + i s T) u` the imaginary quadratic form is
`B_s = D + s (A T + T A) - s² T D T`. The homotopy from the real plane (`s = 0`) to the tilted plane
(`s = 1`) is controlled by the identity

  `B_s = (1 - s) D + s B_1 + s (1 - s) T D T`,

which gives `B_s ≥ s B_1` as quadratic forms when `D` and `T D T` are nonnegative and `0 ≤ s ≤ 1`.
This file proves the identity (in any real module and for matrices) and the quadratic form
inequality.

Nova 045 (endpoint trace). The Cauchy–Schwarz constant in `|v_0 + v_N - 2 μ|² ≤ (2/3)(1 - N⁻²) E`
rests on `∑_{j<N} (2 j + 1 - N)² = N (N² - 1) / 3`, proved here over the integers.

Not formalized here: that the form `B_s` is the imaginary part of the phase on the plane, the
regulator, and the endpoint trace identity itself.
-/

namespace QuantumPathIntegrals

open Matrix

/-- The homotopy identity in a real module, with `X` in the role of `A T + T A` and `Y` in the role
of `T D T`. -/
theorem tilt_homotopy_module {M : Type*} [AddCommGroup M] [Module ℝ M] (D X Y : M) (s : ℝ) :
    D + s • X - s ^ 2 • Y = (1 - s) • D + s • (D + X - Y) + (s * (1 - s)) • Y := by
  module

/-- The homotopy identity for matrices: `B_s = (1 - s) D + s B_1 + s (1 - s) T D T`. -/
theorem tilt_homotopy_matrix {n : Type*} [Fintype n] (A T D : Matrix n n ℝ) (s : ℝ) :
    D + s • (A * T + T * A) - s ^ 2 • (T * D * T)
      = (1 - s) • D + s • (D + (A * T + T * A) - T * D * T) + (s * (1 - s)) • (T * D * T) :=
  tilt_homotopy_module D (A * T + T * A) (T * D * T) s

/-- Scalar form of `B_s ≥ s B_1`. -/
theorem tilt_homotopy_scalar {qD qX qY s : ℝ} (hD : 0 ≤ qD) (hY : 0 ≤ qY) (hs0 : 0 ≤ s)
    (hs1 : s ≤ 1) : s * (qD + qX - qY) ≤ qD + s * qX - s ^ 2 * qY := by
  nlinarith [mul_nonneg (sub_nonneg.2 hs1) hD, mul_nonneg (mul_nonneg hs0 (sub_nonneg.2 hs1)) hY]

/-- `B_s ≥ s B_1` as quadratic forms: if `vᵀ D v ≥ 0` and `vᵀ (T D T) v ≥ 0` and `0 ≤ s ≤ 1`, then
`s · vᵀ B_1 v ≤ vᵀ B_s v`. -/
theorem tilt_homotopy_quadratic {n : Type*} [Fintype n] (A T D : Matrix n n ℝ) (v : n → ℝ) {s : ℝ}
    (hD : 0 ≤ v ⬝ᵥ D *ᵥ v) (hY : 0 ≤ v ⬝ᵥ (T * D * T) *ᵥ v) (hs0 : 0 ≤ s) (hs1 : s ≤ 1) :
    s * (v ⬝ᵥ (D + (A * T + T * A) - T * D * T) *ᵥ v)
      ≤ v ⬝ᵥ (D + s • (A * T + T * A) - s ^ 2 • (T * D * T)) *ᵥ v := by
  have h1 : v ⬝ᵥ (D + (A * T + T * A) - T * D * T) *ᵥ v
      = v ⬝ᵥ D *ᵥ v + v ⬝ᵥ (A * T + T * A) *ᵥ v - v ⬝ᵥ (T * D * T) *ᵥ v := by
    simp only [Matrix.sub_mulVec, Matrix.add_mulVec, dotProduct_sub, dotProduct_add]
  have h2 : v ⬝ᵥ (D + s • (A * T + T * A) - s ^ 2 • (T * D * T)) *ᵥ v
      = v ⬝ᵥ D *ᵥ v + s * (v ⬝ᵥ (A * T + T * A) *ᵥ v) - s ^ 2 * (v ⬝ᵥ (T * D * T) *ᵥ v) := by
    simp only [Matrix.sub_mulVec, Matrix.add_mulVec, Matrix.smul_mulVec, dotProduct_sub,
      dotProduct_add, dotProduct_smul, smul_eq_mul]
  rw [h1, h2]
  exact tilt_homotopy_scalar hD hY hs0 hs1

/-- `3 ∑_{j<n} (2 j + 1 - c)² = n (4 n² - 1) - 6 c n² + 3 n c²`. -/
theorem sum_sq_centered (c : ℤ) (n : ℕ) :
    3 * ∑ j ∈ Finset.range n, (2 * (j : ℤ) + 1 - c) ^ 2
      = n * (4 * n ^ 2 - 1) - 6 * c * n ^ 2 + 3 * n * c ^ 2 := by
  induction n with
  | zero => simp
  | succ k ih =>
    rw [Finset.sum_range_succ, mul_add, ih]
    push_cast
    ring

/-- `3 ∑_{j<N} (2 j + 1 - N)² = N (N² - 1)`, the constant of the endpoint trace bound of Nova 045. -/
theorem sum_sq_endpoint_weights (N : ℕ) :
    3 * ∑ j ∈ Finset.range N, (2 * (j : ℤ) + 1 - (N : ℤ)) ^ 2 = (N : ℤ) * ((N : ℤ) ^ 2 - 1) := by
  rw [sum_sq_centered]
  ring

end QuantumPathIntegrals

#print axioms QuantumPathIntegrals.tilt_homotopy_matrix
#print axioms QuantumPathIntegrals.tilt_homotopy_quadratic
#print axioms QuantumPathIntegrals.sum_sq_endpoint_weights
