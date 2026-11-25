/-
Mean of the M6 null gain (Proposition E.2, first claim; teoria-fase5/teoria-m6-prueba.md).

- `integral_quadForm`: for any noise with second moments `E[ε_i ε_j] = C_ij` (any law,
  Gaussian or not), `E[εᵀ A ε] = tr(A Cᵀ)`.
- `integral_gain_null`: under one event, at fixed times, `E[gain] = tr((P₂ − P₁) Cᵀ)`
  (`= tr(C P₂) − tr(C P₁)` for symmetric C), combining with `gain_null_eq`.
- `trace_mul_sum_vecMulVec`: with `P = ∑_k q_k q_kᵀ` (orthonormal basis of the class),
  `tr(C P) = ∑_k q_kᵀ C q_k`; `trace_eq_lambda_mul_card`: if every basis vector gives the
  same value `λ_r` (columns `e_c ⊗ r(· − t)` of separated pulses in stationary noise), then
  `tr(C P) = λ_r · dim P`, so `m₀ = λ_r (dim P₂ − dim P₁)`.

Hypothesis left as such: that each column gives exactly `λ_r = rᵀCr / rᵀr` (stationarity of
the noise, pulses separated and away from the window edges); it is an input, not derived.
The tail bound of E.2 (Laurent–Massart) is not formalized here.
-/
import TTMT.M6NullGain
import Mathlib.MeasureTheory.Integral.Bochner.Basic
import Mathlib.LinearAlgebra.Matrix.Trace

open Matrix MeasureTheory Finset

namespace TTMT

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

omit [DecidableEq ι] in
theorem integral_quadForm {Ω : Type*} [MeasurableSpace Ω] (μ : Measure Ω)
    (ε : Ω → ι → ℝ) (A C : Matrix ι ι ℝ)
    (hint : ∀ i j, Integrable (fun ω => ε ω i * ε ω j) μ)
    (hC : ∀ i j, ∫ ω, ε ω i * ε ω j ∂μ = C i j) :
    ∫ ω, ε ω ⬝ᵥ (A *ᵥ ε ω) ∂μ = trace (A * Cᵀ) := by
  have h : ∀ ω, ε ω ⬝ᵥ (A *ᵥ ε ω) = ∑ i, ∑ j, A i j * (ε ω i * ε ω j) := by
    intro ω
    simp only [dotProduct, mulVec, Finset.mul_sum]
    refine sum_congr rfl fun i _ => sum_congr rfl fun j _ => by ring
  simp_rw [h]
  rw [integral_finsetSum univ fun i _ =>
    integrable_finsetSum univ fun j _ => (hint i j).const_mul (A i j)]
  have : ∀ i, ∫ ω, ∑ j, A i j * (ε ω i * ε ω j) ∂μ = ∑ j, A i j * C i j := by
    intro i
    rw [integral_finsetSum univ fun j _ => (hint i j).const_mul (A i j)]
    exact sum_congr rfl fun j _ => by rw [integral_const_mul, hC]
  simp_rw [this]
  simp [trace, Matrix.mul_apply, transpose_apply]

/-- **E.2, mean.** Under one event at fixed times, `E[gain] = tr((P₂ − P₁) Cᵀ)`. -/
theorem integral_gain_null {Ω : Type*} [MeasurableSpace Ω] (μ' : Measure Ω)
    (P₁ P₂ : Matrix ι ι ℝ) (hs₁ : P₁ᵀ = P₁) (hi₁ : P₁ * P₁ = P₁) (hs₂ : P₂ᵀ = P₂)
    (hi₂ : P₂ * P₂ = P₂) (μ : ι → ℝ) (hμ₁ : P₁ *ᵥ μ = μ) (hμ₂ : P₂ *ᵥ μ = μ)
    (ε : Ω → ι → ℝ) (C : Matrix ι ι ℝ)
    (hint : ∀ i j, Integrable (fun ω => ε ω i * ε ω j) μ')
    (hC : ∀ i j, ∫ ω, ε ω i * ε ω j ∂μ' = C i j) :
    ∫ ω, (rss P₁ (μ + ε ω) - rss P₂ (μ + ε ω)) ∂μ' = trace ((P₂ - P₁) * Cᵀ) := by
  simp_rw [gain_null_eq P₁ P₂ hs₁ hi₁ hs₂ hi₂ μ _ hμ₁ hμ₂]
  exact integral_quadForm μ' ε _ C hint hC

/-- `tr(C q qᵀ) = qᵀ C q`, summed over a basis. -/
theorem trace_mul_sum_vecMulVec {m : Type*} [Fintype m] (C : Matrix ι ι ℝ) (q : m → ι → ℝ) :
    trace (C * ∑ k, vecMulVec (q k) (q k)) = ∑ k, q k ⬝ᵥ (C *ᵥ q k) := by
  rw [Matrix.mul_sum, trace_sum]
  refine sum_congr rfl fun k _ => ?_
  simp only [trace, Matrix.mul_apply, vecMulVec_apply, Matrix.diag, dotProduct, mulVec, Finset.mul_sum]
  refine sum_congr rfl fun i _ => sum_congr rfl fun j _ => by ring

/-- If every basis vector of the class gives the same `λ_r`, `tr(C P) = λ_r · dim P`. -/
theorem trace_eq_lambda_mul_card {m : Type*} [Fintype m] (C : Matrix ι ι ℝ)
    (q : m → ι → ℝ) (lam : ℝ) (h : ∀ k, q k ⬝ᵥ (C *ᵥ q k) = lam) :
    trace (C * ∑ k, vecMulVec (q k) (q k)) = lam * Fintype.card m := by
  rw [trace_mul_sum_vecMulVec, sum_congr rfl fun k _ => h k, sum_const, card_univ,
    nsmul_eq_mul, mul_comm]

end TTMT
