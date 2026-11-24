import Mathlib

/-!
# Second order remainder of the complex exponential with the factor `max 1 (exp (re w))`

Lemma E-local (correspondence 049) uses, for complex `w`,

  `‖exp w - 1 - w‖ ≤ ‖w‖² / 2 * max 1 (exp (re w))`.

`LocalLemma.lean` has the weaker bound with `exp ‖w‖` in place of `max 1 (exp (re w))`, which is not
sufficient there: on the steepest descent plane `re w` is bounded by a small quadratic form while
`‖w‖` is not. The proof below applies the mean value inequality twice along the segment
`t ↦ exp (t w)`, `0 ≤ t ≤ 1`.

The Gaussian moment inequalities and the change of Gaussian measure used in Lemma E-local are not
formalized here.
-/

namespace QuantumPathIntegrals

open Set

/-- For `0 ≤ t ≤ 1`, `exp (t x) ≤ max 1 (exp x)`. -/
theorem exp_mul_le_max_one_exp {x t : ℝ} (ht0 : 0 ≤ t) (ht1 : t ≤ 1) :
    Real.exp (t * x) ≤ max 1 (Real.exp x) := by
  rcases le_total 0 x with hx | hx
  · have h : t * x ≤ x := by nlinarith [mul_nonneg (sub_nonneg.mpr ht1) hx]
    exact (Real.exp_le_exp.mpr h).trans (le_max_right _ _)
  · have h : t * x ≤ 0 := mul_nonpos_of_nonneg_of_nonpos ht0 hx
    exact (Real.exp_le_one_iff.mpr h).trans (le_max_left _ _)

/-- Derivative of `s ↦ s w` for real `s` and complex `w`. -/
theorem hasDerivAt_ofReal_mul_const (w : ℂ) (t : ℝ) :
    HasDerivAt (fun s : ℝ => (s : ℂ) * w) w t := by
  have h1 : HasDerivAt (fun s : ℝ => (s : ℂ)) (1 : ℂ) t := by
    simpa using (hasDerivAt_id t).ofReal_comp
  simpa using h1.mul_const w

/-- Derivative of `s ↦ exp (s w)` for real `s` and complex `w`. -/
theorem hasDerivAt_cexp_ofReal_mul (w : ℂ) (t : ℝ) :
    HasDerivAt (fun s : ℝ => Complex.exp ((s : ℂ) * w)) (Complex.exp ((t : ℂ) * w) * w) t :=
  (hasDerivAt_ofReal_mul_const w t).cexp

/-- For `0 ≤ t ≤ 1`, `‖exp (t w)‖ ≤ max 1 (exp (re w))`. -/
theorem norm_cexp_ofReal_mul_le (w : ℂ) {t : ℝ} (ht0 : 0 ≤ t) (ht1 : t ≤ 1) :
    ‖Complex.exp ((t : ℂ) * w)‖ ≤ max 1 (Real.exp w.re) := by
  rw [Complex.norm_exp, Complex.re_ofReal_mul]
  exact exp_mul_le_max_one_exp ht0 ht1

/-- For `0 ≤ t ≤ 1`, `‖exp (t w) - 1‖ ≤ ‖w‖ max 1 (exp (re w)) t`. -/
theorem norm_cexp_ofReal_mul_sub_one_le (w : ℂ) {t : ℝ} (ht : t ∈ Icc (0 : ℝ) 1) :
    ‖Complex.exp ((t : ℂ) * w) - 1‖ ≤ ‖w‖ * max 1 (Real.exp w.re) * t := by
  have hd : ∀ s : ℝ, HasDerivAt (fun s : ℝ => Complex.exp ((s : ℂ) * w) - 1)
      (Complex.exp ((s : ℂ) * w) * w) s := fun s => (hasDerivAt_cexp_ofReal_mul w s).sub_const 1
  have hB : ∀ s : ℝ, HasDerivAt (fun s : ℝ => ‖w‖ * max 1 (Real.exp w.re) * s)
      (‖w‖ * max 1 (Real.exp w.re) * 1) s := fun s =>
    (hasDerivAt_id' s).const_mul (‖w‖ * max 1 (Real.exp w.re))
  have hcont : ContinuousOn (fun s : ℝ => Complex.exp ((s : ℂ) * w) - 1) (Icc 0 1) :=
    (continuous_iff_continuousAt.2 fun s => (hd s).continuousAt).continuousOn
  have h0 : ‖Complex.exp (((0 : ℝ) : ℂ) * w) - 1‖ ≤ ‖w‖ * max 1 (Real.exp w.re) * 0 := by simp
  have hbound : ∀ s ∈ Ico (0 : ℝ) 1,
      ‖Complex.exp ((s : ℂ) * w) * w‖ ≤ ‖w‖ * max 1 (Real.exp w.re) * 1 := by
    intro s hs
    calc ‖Complex.exp ((s : ℂ) * w) * w‖ = ‖w‖ * ‖Complex.exp ((s : ℂ) * w)‖ := by
          rw [norm_mul, mul_comm]
      _ ≤ ‖w‖ * max 1 (Real.exp w.re) :=
          mul_le_mul_of_nonneg_left (norm_cexp_ofReal_mul_le w hs.1 hs.2.le) (norm_nonneg w)
      _ = ‖w‖ * max 1 (Real.exp w.re) * 1 := by ring
  exact image_norm_le_of_norm_deriv_right_le_deriv_boundary hcont
    (fun s _ => (hd s).hasDerivWithinAt) h0 hB hbound ht

/-- For complex `w`, `‖exp w - 1 - w‖ ≤ ‖w‖² / 2 * max 1 (exp (re w))`. -/
theorem norm_cexp_sub_one_sub_id_le (w : ℂ) :
    ‖Complex.exp w - 1 - w‖ ≤ ‖w‖ ^ 2 / 2 * max 1 (Real.exp w.re) := by
  have hd : ∀ s : ℝ, HasDerivAt (fun s : ℝ => Complex.exp ((s : ℂ) * w) - 1 - (s : ℂ) * w)
      (Complex.exp ((s : ℂ) * w) * w - w) s := by
    intro s
    have h := ((hasDerivAt_cexp_ofReal_mul w s).sub_const 1).sub (hasDerivAt_ofReal_mul_const w s)
    exact h
  have hB : ∀ s : ℝ, HasDerivAt (fun s : ℝ => ‖w‖ ^ 2 * max 1 (Real.exp w.re) / 2 * (s * s))
      (‖w‖ ^ 2 * max 1 (Real.exp w.re) / 2 * (1 * s + s * 1)) s := by
    intro s
    have hid : HasDerivAt (fun s : ℝ => s) 1 s := hasDerivAt_id' s
    have hsq : HasDerivAt (fun s : ℝ => s * s) (1 * s + s * 1) s := hid.mul hid
    exact hsq.const_mul (‖w‖ ^ 2 * max 1 (Real.exp w.re) / 2)
  have hcont : ContinuousOn (fun s : ℝ => Complex.exp ((s : ℂ) * w) - 1 - (s : ℂ) * w) (Icc 0 1) :=
    (continuous_iff_continuousAt.2 fun s => (hd s).continuousAt).continuousOn
  have h0 : ‖Complex.exp (((0 : ℝ) : ℂ) * w) - 1 - ((0 : ℝ) : ℂ) * w‖
      ≤ ‖w‖ ^ 2 * max 1 (Real.exp w.re) / 2 * (0 * 0) := by simp
  have hbound : ∀ s ∈ Ico (0 : ℝ) 1, ‖Complex.exp ((s : ℂ) * w) * w - w‖
      ≤ ‖w‖ ^ 2 * max 1 (Real.exp w.re) / 2 * (1 * s + s * 1) := by
    intro s hs
    have h1 := norm_cexp_ofReal_mul_sub_one_le w (Ico_subset_Icc_self hs)
    have hfac : Complex.exp ((s : ℂ) * w) * w - w = w * (Complex.exp ((s : ℂ) * w) - 1) := by ring
    rw [hfac, norm_mul]
    calc ‖w‖ * ‖Complex.exp ((s : ℂ) * w) - 1‖ ≤ ‖w‖ * (‖w‖ * max 1 (Real.exp w.re) * s) :=
          mul_le_mul_of_nonneg_left h1 (norm_nonneg w)
      _ = ‖w‖ ^ 2 * max 1 (Real.exp w.re) / 2 * (1 * s + s * 1) := by ring
  have key : ‖Complex.exp (((1 : ℝ) : ℂ) * w) - 1 - ((1 : ℝ) : ℂ) * w‖
      ≤ ‖w‖ ^ 2 * max 1 (Real.exp w.re) / 2 * (1 * 1) :=
    image_norm_le_of_norm_deriv_right_le_deriv_boundary hcont
      (fun s _ => (hd s).hasDerivWithinAt) h0 hB hbound (right_mem_Icc.2 zero_le_one)
  have e1 : Complex.exp (((1 : ℝ) : ℂ) * w) - 1 - ((1 : ℝ) : ℂ) * w = Complex.exp w - 1 - w := by
    simp
  rw [e1] at key
  calc ‖Complex.exp w - 1 - w‖ ≤ ‖w‖ ^ 2 * max 1 (Real.exp w.re) / 2 * (1 * 1) := key
    _ = ‖w‖ ^ 2 / 2 * max 1 (Real.exp w.re) := by ring

/-- The form used in Lemma E-local: if `re w ≤ b` then
`‖exp w - 1 - w‖ ≤ ‖w‖² / 2 * max 1 (exp b)`. -/
theorem norm_cexp_sub_one_sub_id_le_of_re_le {w : ℂ} {b : ℝ} (hb : w.re ≤ b) :
    ‖Complex.exp w - 1 - w‖ ≤ ‖w‖ ^ 2 / 2 * max 1 (Real.exp b) := by
  refine (norm_cexp_sub_one_sub_id_le w).trans ?_
  have hmax : max 1 (Real.exp w.re) ≤ max 1 (Real.exp b) :=
    max_le_max le_rfl (Real.exp_le_exp.mpr hb)
  exact mul_le_mul_of_nonneg_left hmax (by positivity)

end QuantumPathIntegrals

#print axioms QuantumPathIntegrals.norm_cexp_sub_one_sub_id_le
#print axioms QuantumPathIntegrals.norm_cexp_sub_one_sub_id_le_of_re_le
