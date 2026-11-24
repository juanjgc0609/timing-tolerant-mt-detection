/-
Deterministic core of Proposition E.2 (theory of the M6 "one against two" test; project
file teoria-fase5/teoria-m6-prueba.md, not yet in the repository appendix).

Under one event (the signal μ₀ lies in both classes) and at FIXED times, the M6 gain
`gain = RSS₁ − RSS₂` equals the noise quadratic form `εᵀ (P₂ − P₁) ε` and does not depend
on the signal. This is what makes the null law of the gain pivotal (false-alarm control
of M6 at fixed times).

Formal hypotheses: `P₁`, `P₂` symmetric and idempotent (orthogonal projectors) and
`P₁ μ = μ`, `P₂ μ = μ`. Neither class is assumed to contain the other (in the code the
free one-event class has dimension 48 and the pair class 12). The noise `ε` is arbitrary,
so the identity holds for any covariance C. NOT formalized here: the probabilistic part of
E.2 (mean tr(C P₂) − tr(C P₁), Laurent–Massart tail) and the inflation from the search over
times, which E.2 itself does not cover.
-/
import Mathlib.Basic.Real.Basic
import Mathlib.Data.Matrix.Mul
import Mathlib.LinearAlgebra.Matrix.Symmetric

open Matrix

namespace TTMT

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- Residual sum of squares after projecting with `P`. -/
def rss (P : Matrix ι ι ℝ) (y : ι → ℝ) : ℝ :=
  ((1 - P) *ᵥ y) ⬝ᵥ ((1 - P) *ᵥ y)

/-- For an orthogonal projector, `RSS = yᵀ (I − P) y`. -/
theorem rss_eq_quadForm (P : Matrix ι ι ℝ) (hs : Pᵀ = P) (hi : P * P = P) (y : ι → ℝ) :
    rss P y = y ⬝ᵥ ((1 - P) *ᵥ y) := by
  have hsym : (1 - P)ᵀ = 1 - P := by rw [transpose_sub, transpose_one, hs]
  have hid : (1 - P) * (1 - P) = 1 - P := by
    rw [sub_mul, mul_sub, mul_sub, one_mul, mul_one, one_mul, hi, sub_self, sub_zero]
  calc rss P y = ((1 - P) *ᵥ y) ᵥ* (1 - P) ⬝ᵥ y := dotProduct_mulVec _ _ _
    _ = (1 - P)ᵀ *ᵥ ((1 - P) *ᵥ y) ⬝ᵥ y := by rw [mulVec_transpose]
    _ = (1 - P) *ᵥ y ⬝ᵥ y := by rw [hsym, mulVec_mulVec, hid]
    _ = y ⬝ᵥ ((1 - P) *ᵥ y) := dotProduct_comm _ _

/-- A signal in the class leaves no residual: `(I − P)(μ + ε) = (I − P) ε`. -/
theorem resid_signal_free (P : Matrix ι ι ℝ) (μ ε : ι → ℝ) (hμ : P *ᵥ μ = μ) :
    (1 - P) *ᵥ (μ + ε) = (1 - P) *ᵥ ε := by
  rw [mulVec_add, sub_mulVec, one_mulVec, hμ, sub_self, zero_add]

/-- **Proposition E.2, deterministic part.** Under one event (μ in both classes) the M6
gain is `εᵀ (P₂ − P₁) ε`, with no dependence on μ. -/
theorem gain_null_eq (P₁ P₂ : Matrix ι ι ℝ) (hs₁ : P₁ᵀ = P₁) (hi₁ : P₁ * P₁ = P₁)
    (hs₂ : P₂ᵀ = P₂) (hi₂ : P₂ * P₂ = P₂) (μ ε : ι → ℝ)
    (hμ₁ : P₁ *ᵥ μ = μ) (hμ₂ : P₂ *ᵥ μ = μ) :
    rss P₁ (μ + ε) - rss P₂ (μ + ε) = ε ⬝ᵥ ((P₂ - P₁) *ᵥ ε) := by
  have h₁ : rss P₁ (μ + ε) = rss P₁ ε := by unfold rss; rw [resid_signal_free P₁ μ ε hμ₁]
  have h₂ : rss P₂ (μ + ε) = rss P₂ ε := by unfold rss; rw [resid_signal_free P₂ μ ε hμ₂]
  rw [h₁, h₂, rss_eq_quadForm P₁ hs₁ hi₁, rss_eq_quadForm P₂ hs₂ hi₂, ← dotProduct_sub,
    ← sub_mulVec]
  congr 2
  rw [sub_sub_sub_cancel_left]

/-- Corollary: any two signals of the common class give the same gain with the same noise
(the null law of the gain is the same for every event). -/
theorem gain_null_pivotal (P₁ P₂ : Matrix ι ι ℝ) (hs₁ : P₁ᵀ = P₁) (hi₁ : P₁ * P₁ = P₁)
    (hs₂ : P₂ᵀ = P₂) (hi₂ : P₂ * P₂ = P₂) (μ μ' ε : ι → ℝ)
    (hμ₁ : P₁ *ᵥ μ = μ) (hμ₂ : P₂ *ᵥ μ = μ) (hμ₁' : P₁ *ᵥ μ' = μ') (hμ₂' : P₂ *ᵥ μ' = μ') :
    rss P₁ (μ + ε) - rss P₂ (μ + ε) = rss P₁ (μ' + ε) - rss P₂ (μ' + ε) := by
  rw [gain_null_eq P₁ P₂ hs₁ hi₁ hs₂ hi₂ μ ε hμ₁ hμ₂,
    gain_null_eq P₁ P₂ hs₁ hi₁ hs₂ hi₂ μ' ε hμ₁' hμ₂']

end TTMT
