import Mathlib

/-!
# Assembly of the selection certificate

Theorem S of correspondence 008 has three analytic inputs (the contour lemma, the local lemma on each
saddle segment, and the tail bound) and one purely logical step that combines them. This file
formalizes the logical step and the necessity statement in the form accepted in correspondence 010.

The analytic inputs are hypotheses here. They are not proved in Lean.
-/

namespace QuantumPathIntegrals

/-- If the integral splits into saddle-segment pieces `S i` plus a tail, each piece is within `E i`
of its leading term `L i`, and the tail is at most `T`, then the integral is within `∑ E i + T` of
the sum of the leading terms. -/
theorem certificate_sum {ι : Type*} (s : Finset ι) (I Tl : ℂ) (S L : ι → ℂ) (E : ι → ℝ) (T : ℝ)
    (hI : I = ∑ i ∈ s, S i + Tl) (hloc : ∀ i ∈ s, ‖S i - L i‖ ≤ E i) (hT : ‖Tl‖ ≤ T) :
    ‖I - ∑ i ∈ s, L i‖ ≤ ∑ i ∈ s, E i + T := by
  have hsplit : I - ∑ i ∈ s, L i = ∑ i ∈ s, (S i - L i) + Tl := by
    rw [hI, Finset.sum_sub_distrib]; ring
  rw [hsplit]
  calc ‖∑ i ∈ s, (S i - L i) + Tl‖
      ≤ ‖∑ i ∈ s, (S i - L i)‖ + ‖Tl‖ := norm_add_le _ _
    _ ≤ ∑ i ∈ s, ‖S i - L i‖ + ‖Tl‖ := add_le_add (norm_sum_le _ _) le_rfl
    _ ≤ ∑ i ∈ s, E i + T := add_le_add (Finset.sum_le_sum hloc) hT

/-- Two retained saddles, the case used in correspondences 008 and 010. -/
theorem certificate_two (I S₁ S₂ Tl L₁ L₂ : ℂ) (E₁ E₂ T : ℝ)
    (hI : I = S₁ + S₂ + Tl) (h₁ : ‖S₁ - L₁‖ ≤ E₁) (h₂ : ‖S₂ - L₂‖ ≤ E₂) (hT : ‖Tl‖ ≤ T) :
    ‖I - (L₁ + L₂)‖ ≤ E₁ + E₂ + T := by
  have hsplit : I - (L₁ + L₂) = (S₁ - L₁) + (S₂ - L₂) + Tl := by rw [hI]; ring
  rw [hsplit]
  have h3 := norm_add_le ((S₁ - L₁) + (S₂ - L₂)) Tl
  have h4 := norm_add_le (S₁ - L₁) (S₂ - L₂)
  linarith

/-- Necessity in the qualified form: if the two-term formula has error at most `B`, then the
approximation by the first leading term alone has error at least `‖L₂‖ - B`. -/
theorem first_term_alone_error_ge (I L₁ L₂ : ℂ) (B : ℝ) (h : ‖I - (L₁ + L₂)‖ ≤ B) :
    ‖L₂‖ - B ≤ ‖I - L₁‖ := by
  have h1 : L₂ = (I - L₁) - (I - (L₁ + L₂)) := by ring
  have h2 : ‖L₂‖ ≤ ‖I - L₁‖ + ‖I - (L₁ + L₂)‖ := by
    calc ‖L₂‖ = ‖(I - L₁) - (I - (L₁ + L₂))‖ := congrArg norm h1
      _ ≤ ‖I - L₁‖ + ‖I - (L₁ + L₂)‖ := norm_sub_le _ _
  linarith

end QuantumPathIntegrals

#print axioms QuantumPathIntegrals.certificate_sum
#print axioms QuantumPathIntegrals.certificate_two
#print axioms QuantumPathIntegrals.first_term_alone_error_ge
