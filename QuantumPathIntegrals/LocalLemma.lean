import Mathlib

/-!
# Pointwise core of the local lemma (Lemma L, correspondence 008)

On a saddle segment the integrand is `exp (-β) * exp w` with `β = a t² / (2ħ) ≥ 0` and a complex
remainder `w`. Lemma L writes `exp w = 1 + q + (w - q) + (exp w - 1 - w)` and bounds the last term by
`‖w‖² / 2 * exp ‖w‖`. This file proves that bound, which Mathlib has only for `‖w‖ ≤ 1` and without
the factor `1/2`, and the resulting pointwise inequality that is then integrated.

The Gaussian moment integrals that turn the pointwise inequality into the constants `1/8` and `5/24`
are not formalized here.
-/

namespace QuantumPathIntegrals

/-- For `r ≥ 0`, `exp r - 1 - r ≤ r² / 2 * exp r`. -/
theorem exp_sub_one_sub_id_le_sq_mul_exp {r : ℝ} (hr : 0 ≤ r) :
    Real.exp r - 1 - r ≤ r ^ 2 / 2 * Real.exp r := by
  have hderiv : ∀ s : ℝ, HasDerivAt (fun s : ℝ => s * s / 2 - 1 + Real.exp (-s) * (1 + s))
      (s - s * Real.exp (-s)) s := by
    intro s
    have hid : HasDerivAt (fun s : ℝ => s) 1 s := hasDerivAt_id' s
    have hsq : HasDerivAt (fun s : ℝ => s * s) (1 * s + s * 1) s := hid.mul hid
    have h1 : HasDerivAt (fun s : ℝ => s * s / 2 - 1) ((1 * s + s * 1) / 2) s :=
      (hsq.div_const 2).sub_const 1
    have hneg : HasDerivAt (fun s : ℝ => -s) (-1) s := hid.neg
    have hexp : HasDerivAt (fun s : ℝ => Real.exp (-s)) (Real.exp (-s) * (-1)) s := hneg.exp
    have hlin : HasDerivAt (fun s : ℝ => 1 + s) 1 s := hid.const_add 1
    have hprod := hexp.mul hlin
    have hG := h1.add hprod
    convert hG using 1
    ring
  have hmono : MonotoneOn (fun s : ℝ => s * s / 2 - 1 + Real.exp (-s) * (1 + s)) (Set.Ici 0) := by
    apply monotoneOn_of_deriv_nonneg (convex_Ici 0)
    · exact (continuous_iff_continuousAt.mpr fun s => (hderiv s).continuousAt).continuousOn
    · exact fun s _ => (hderiv s).differentiableAt.differentiableWithinAt
    · intro s hs
      rw [interior_Ici] at hs
      have hs0 : 0 < s := hs
      rw [(hderiv s).deriv]
      have hle : Real.exp (-s) ≤ 1 := by
        rw [Real.exp_le_one_iff]; linarith
      nlinarith
  have h0 := hmono (Set.mem_Ici.mpr le_rfl) (Set.mem_Ici.mpr hr) hr
  have e0 : (0 : ℝ) * 0 / 2 - 1 + Real.exp (-0) * (1 + 0) = 0 := by norm_num
  beta_reduce at h0
  rw [e0] at h0
  have hE : Real.exp r * Real.exp (-r) = 1 := by
    rw [← Real.exp_add]; simp
  have key : 0 ≤ Real.exp r * (r * r / 2 - 1 + Real.exp (-r) * (1 + r)) :=
    mul_nonneg (Real.exp_pos r).le h0
  have expand : Real.exp r * (r * r / 2 - 1 + Real.exp (-r) * (1 + r))
      = r ^ 2 / 2 * Real.exp r - Real.exp r + (Real.exp r * Real.exp (-r)) * (1 + r) := by ring
  rw [expand, hE] at key
  linarith

/-- Second-order remainder of the complex exponential, valid for every `w`. -/
theorem norm_exp_sub_one_sub_id_le_sq (w : ℂ) :
    ‖Complex.exp w - 1 - w‖ ≤ ‖w‖ ^ 2 / 2 * Real.exp ‖w‖ := by
  have h := Complex.norm_exp_sub_sum_le_exp_norm_sub_sum w 2
  have hs1 : ∑ m ∈ Finset.range 2, w ^ m / (m.factorial : ℂ) = 1 + w := by
    simp [Finset.sum_range_succ]
  have hs2 : ∑ m ∈ Finset.range 2, ‖w‖ ^ m / (m.factorial : ℝ) = 1 + ‖w‖ := by
    simp [Finset.sum_range_succ]
  rw [hs1, hs2] at h
  have h2 := exp_sub_one_sub_id_le_sq_mul_exp (norm_nonneg w)
  have e : Complex.exp w - 1 - w = Complex.exp w - (1 + w) := by ring
  rw [e]
  linarith

/-- The pointwise inequality behind Lemma L. With `β ≥ 0` the Gaussian exponent, `w` the complex
remainder, `q` its odd cubic part, `‖w - q‖ ≤ c4` and `‖w‖ ≤ κ β`:
`exp (-β) ‖exp w - 1 - q‖ ≤ exp (-β) c4 + ‖w‖² / 2 * exp (-(1 - κ) β)`. -/
theorem local_pointwise (w q : ℂ) (β κ c4 : ℝ) (hw : ‖w‖ ≤ κ * β) (hq : ‖w - q‖ ≤ c4) :
    Real.exp (-β) * ‖Complex.exp w - 1 - q‖
      ≤ Real.exp (-β) * c4 + ‖w‖ ^ 2 / 2 * Real.exp (-(1 - κ) * β) := by
  have h1 := norm_exp_sub_one_sub_id_le_sq w
  have h2 : ‖Complex.exp w - 1 - q‖ ≤ ‖Complex.exp w - 1 - w‖ + ‖w - q‖ := by
    have hsplit : Complex.exp w - 1 - q = (Complex.exp w - 1 - w) + (w - q) := by ring
    rw [hsplit]; exact norm_add_le _ _
  have h3 : Real.exp ‖w‖ ≤ Real.exp (κ * β) := Real.exp_le_exp.mpr hw
  have he : Real.exp (-β) * Real.exp (κ * β) = Real.exp (-(1 - κ) * β) := by
    rw [← Real.exp_add]; congr 1; ring
  have hpos : 0 ≤ Real.exp (-β) := (Real.exp_pos _).le
  have hw2 : 0 ≤ ‖w‖ ^ 2 / 2 := by positivity
  have h4 : ‖w‖ ^ 2 / 2 * Real.exp ‖w‖ ≤ ‖w‖ ^ 2 / 2 * Real.exp (κ * β) :=
    mul_le_mul_of_nonneg_left h3 hw2
  have h5 : ‖Complex.exp w - 1 - q‖ ≤ c4 + ‖w‖ ^ 2 / 2 * Real.exp (κ * β) := by linarith
  calc Real.exp (-β) * ‖Complex.exp w - 1 - q‖
      ≤ Real.exp (-β) * (c4 + ‖w‖ ^ 2 / 2 * Real.exp (κ * β)) := mul_le_mul_of_nonneg_left h5 hpos
    _ = Real.exp (-β) * c4 + ‖w‖ ^ 2 / 2 * (Real.exp (-β) * Real.exp (κ * β)) := by ring
    _ = Real.exp (-β) * c4 + ‖w‖ ^ 2 / 2 * Real.exp (-(1 - κ) * β) := by rw [he]

end QuantumPathIntegrals

#print axioms QuantumPathIntegrals.exp_sub_one_sub_id_le_sq_mul_exp
#print axioms QuantumPathIntegrals.norm_exp_sub_one_sub_id_le_sq
#print axioms QuantumPathIntegrals.local_pointwise
