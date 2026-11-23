import Mathlib

/-!
# The three-term recurrence behind the discrete harmonic determinant

For the time-sliced harmonic oscillator the interior Hessian is
`(m/Δt) · tridiag(-1, 2 cos θ, -1)`. Its scaled leading principal minors `D n`
satisfy `D (n+2) = 2 cos θ · D (n+1) - D n` with `D 0 = 1`, `D 1 = 2 cos θ`.
This file proves the closed form `D n · sin θ = sin ((n+1) θ)`, which gives
`det = (m/Δt)^(N-1) · sin (N θ) / sin θ` and locates the discrete caustics at
`N θ ∈ π ℤ`.

This is the algebraic core of Lemma A (correspondence 002, Section 3). It does
not yet connect `D` to `Matrix.det`; that link is a separate statement.
-/

namespace QuantumPathIntegrals

open Real

/-- Closed form of the Chebyshev-type recurrence for the scaled determinants. -/
theorem scaledDet_mul_sin (θ : ℝ) (D : ℕ → ℝ)
    (h0 : D 0 = 1) (h1 : D 1 = 2 * cos θ)
    (hrec : ∀ n : ℕ, D (n + 2) = 2 * cos θ * D (n + 1) - D n) (n : ℕ) :
    D n * sin θ = sin (((n : ℝ) + 1) * θ) := by
  have key : ∀ k : ℕ, D k * sin θ = sin (((k : ℝ) + 1) * θ) ∧
      D (k + 1) * sin θ = sin (((k : ℝ) + 2) * θ) := by
    intro k
    induction k with
    | zero =>
      constructor
      · rw [h0]; norm_num
      · rw [show (0 : ℕ) + 1 = 1 from rfl, h1]
        norm_num
        rw [sin_two_mul]; ring
    | succ k ih =>
      obtain ⟨ih1, ih2⟩ := ih
      constructor
      · rw [ih2]; congr 1; push_cast; ring
      · have hstep : D (k + 1 + 1) = 2 * cos θ * D (k + 1) - D k := hrec k
        have e1 : (((k + 1 : ℕ) : ℝ) + 2) * θ = ((k : ℝ) + 2) * θ + θ := by
          push_cast; ring
        have e2 : ((k : ℝ) + 1) * θ = ((k : ℝ) + 2) * θ - θ := by ring
        calc D (k + 1 + 1) * sin θ
            = 2 * cos θ * (D (k + 1) * sin θ) - D k * sin θ := by rw [hstep]; ring
          _ = 2 * cos θ * sin (((k : ℝ) + 2) * θ) - sin (((k : ℝ) + 1) * θ) := by
              rw [ih2, ih1]
          _ = sin ((((k + 1 : ℕ) : ℝ) + 2) * θ) := by
              rw [e1, e2, sin_add, sin_sub]; ring
  exact (key n).1

/-- Away from `sin θ = 0` the scaled determinant is `sin ((n+1) θ) / sin θ`. -/
theorem scaledDet_eq (θ : ℝ) (D : ℕ → ℝ)
    (h0 : D 0 = 1) (h1 : D 1 = 2 * cos θ)
    (hrec : ∀ n : ℕ, D (n + 2) = 2 * cos θ * D (n + 1) - D n)
    (hθ : sin θ ≠ 0) (n : ℕ) :
    D n = sin (((n : ℝ) + 1) * θ) / sin θ := by
  rw [eq_div_iff hθ]
  exact scaledDet_mul_sin θ D h0 h1 hrec n

/-- The scaled determinant vanishes exactly at the discrete caustics. -/
theorem scaledDet_eq_zero_iff (θ : ℝ) (D : ℕ → ℝ)
    (h0 : D 0 = 1) (h1 : D 1 = 2 * cos θ)
    (hrec : ∀ n : ℕ, D (n + 2) = 2 * cos θ * D (n + 1) - D n)
    (hθ : sin θ ≠ 0) (n : ℕ) :
    D n = 0 ↔ sin (((n : ℝ) + 1) * θ) = 0 := by
  rw [scaledDet_eq θ D h0 h1 hrec hθ n, div_eq_zero_iff]
  constructor
  · rintro (h | h)
    · exact h
    · exact absurd h hθ
  · intro h; exact Or.inl h

end QuantumPathIntegrals

#print axioms QuantumPathIntegrals.scaledDet_mul_sin
#print axioms QuantumPathIntegrals.scaledDet_eq
#print axioms QuantumPathIntegrals.scaledDet_eq_zero_iff
