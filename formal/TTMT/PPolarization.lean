/-
Proposition AB.A.9 (docs/proposal/Supplementary_Evidence_Appendix_AB.md; first stated as
lemma F3 of the 2026-11-20 formulation audit).

In the far-field operator the P amplitude of a moment tensor `M` at station i is
`g_iP(M) = u_i · (n_srcᵀ M n_src) · c_iP`: a FIXED vector `u_i` (the ray direction at the
receiver, or with a free surface a fixed vector in the radial-vertical plane) times a scalar
that is the only place where `M` enters. Hence `g¹_iP × g²_iP = 0` for every pair of
tensors and `κ_iP = ‖g¹_iP × g²_iP‖² = 0` in Corollary AB.A.6.

Formal hypotheses: none beyond the form of the operator. `u`, `n_src`, `c` and the two
matrices are arbitrary (M need not be symmetric, `n_src` need not be unit).
-/
import Mathlib.Basic.Real.Basic
import Mathlib.LinearAlgebra.CrossProduct

open Matrix

namespace TTMT

/-- Far-field P amplitude: `(c · n_srcᵀ M n_src) • u`. -/
def pAmplitude (u nsrc : Fin 3 → ℝ) (c : ℝ) (M : Matrix (Fin 3) (Fin 3) ℝ) : Fin 3 → ℝ :=
  (c * (nsrc ⬝ᵥ (M *ᵥ nsrc))) • u

/-- The P amplitude is a scalar multiple of the fixed vector `u`, whatever the tensor. -/
theorem pAmplitude_parallel (u nsrc : Fin 3 → ℝ) (c : ℝ) (M : Matrix (Fin 3) (Fin 3) ℝ) :
    ∃ s : ℝ, pAmplitude u nsrc c M = s • u :=
  ⟨_, rfl⟩

/-- **Proposition AB.A.9.** Two P amplitudes at one station have zero cross product. -/
theorem pAmplitude_cross_eq_zero (u nsrc : Fin 3 → ℝ) (c : ℝ)
    (M₁ M₂ : Matrix (Fin 3) (Fin 3) ℝ) :
    pAmplitude u nsrc c M₁ ⨯₃ pAmplitude u nsrc c M₂ = 0 := by
  simp [pAmplitude, cross_self]

/-- Corollary: the mechanism-diversity term of the P arrival, `‖g¹ × g²‖²`, is zero. -/
theorem kappa_P_eq_zero (u nsrc : Fin 3 → ℝ) (c : ℝ) (M₁ M₂ : Matrix (Fin 3) (Fin 3) ℝ) :
    let w := pAmplitude u nsrc c M₁ ⨯₃ pAmplitude u nsrc c M₂
    w ⬝ᵥ w = 0 := by
  simp only [pAmplitude_cross_eq_zero, dotProduct_zero]

end TTMT
