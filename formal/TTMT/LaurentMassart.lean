/-
Laurent–Massart tail bounds for weighted chi-square sums (Laurent and Massart, 2000, Lemma 1),
used in the tail of Proposition E.2 (teoria-fase5/teoria-m6-prueba.md).

Let `X_j` be independent standard Gaussians, `λ_j ≥ 0`, `v = ∑ λ_j²`, `a ≥ max λ_j`, and
`Z = ∑_j λ_j (X_j² − 1)`. For every `τ > 0` (with `v > 0`):
- `laurent_massart_upper`: `P(Z ≥ 2√(vτ) + 2aτ) ≤ e^{−τ}`;
- `laurent_massart_lower`: `P(Z ≤ −2√(vτ)) ≤ e^{−τ}`.

Proof (Chernoff): `E e^{s λ(X²−1)} = e^{−sλ}/√(1−2sλ)`, bounded by `exp(s²λ²/(1−2sλ))`
(`neg_log_one_sub_le`) for the upper tail and by `exp(s²λ²)` (`log_one_add_ge`) for the lower
one; independence multiplies the MGFs; the choices `s = √τ/(√v + 2a√τ)` and `s = √τ/√v`
make the Chernoff exponent exactly `−τ`.
-/
import TTMT.LogBounds
import Mathlib.Probability.Moments.Basic
import Mathlib.Probability.Distributions.Gaussian.Real
import Mathlib.Analysis.SpecialFunctions.Gaussian.GaussianIntegral

open MeasureTheory ProbabilityTheory Real Finset

namespace TTMT

/-- `E[exp(s X²)] = 1/√(1−2s)` for `X ~ N(0,1)` and `s < 1/2`. -/
theorem integral_exp_mul_sq_gaussian {s : ℝ} (hs : s < 1 / 2) :
    ∫ x, exp (s * x ^ 2) ∂(gaussianReal 0 1) = 1 / √(1 - 2 * s) := by
  rw [integral_gaussianReal_eq_integral_smul one_ne_zero]
  have h : ∀ x : ℝ, gaussianPDFReal 0 1 x • exp (s * x ^ 2) =
      (√(2 * π))⁻¹ * exp (-(1 / 2 - s) * x ^ 2) := by
    intro x
    simp only [gaussianPDFReal, NNReal.coe_one, mul_one, sub_zero, smul_eq_mul]
    rw [mul_assoc, ← exp_add]
    congr 2; ring
  simp_rw [h]
  rw [integral_const_mul, integral_gaussian]
  have hb : 0 < 1 / 2 - s := by linarith
  have h2 : 0 < 1 - 2 * s := by linarith
  rw [← Real.sqrt_inv, ← Real.sqrt_mul (by positivity), one_div (√(1 - 2 * s)), ← Real.sqrt_inv]
  congr 1
  field_simp

theorem integrable_exp_mul_sq_gaussian {s : ℝ} (hs : s < 1 / 2) :
    Integrable (fun x => exp (s * x ^ 2)) (gaussianReal 0 1) := by
  refine Integrable.of_integral_ne_zero ?_
  rw [integral_exp_mul_sq_gaussian hs]
  have : 0 < 1 - 2 * s := by linarith
  positivity

/-- The MGF of `λ(X² − 1)` at `s` with `sλ < 1/2`. -/
theorem integral_exp_weighted_chisq {s lam : ℝ} (h : s * lam < 1 / 2) :
    ∫ x, exp (s * (lam * (x ^ 2 - 1))) ∂(gaussianReal 0 1) =
      exp (-(s * lam)) * (1 / √(1 - 2 * (s * lam))) := by
  have e : ∀ x : ℝ, exp (s * (lam * (x ^ 2 - 1))) = exp (-(s * lam)) * exp ((s * lam) * x ^ 2) := by
    intro x; rw [← exp_add]; congr 1; ring
  simp_rw [e]
  rw [integral_const_mul, integral_exp_mul_sq_gaussian h]

/-- Upper-tail MGF bound: `e^{−y}/√(1−2y) ≤ exp(y²/(1−2y))` for `0 ≤ y < 1/2`. -/
theorem mgf_bound_upper {y : ℝ} (hy0 : 0 ≤ y) (hy1 : y < 1 / 2) :
    exp (-y) * (1 / √(1 - 2 * y)) ≤ exp (y ^ 2 / (1 - 2 * y)) := by
  have hpos : 0 < 1 - 2 * y := by linarith
  have hs : 1 / √(1 - 2 * y) = exp (-(1 / 2) * log (1 - 2 * y)) := by
    rw [Real.sqrt_eq_rpow, Real.rpow_def_of_pos hpos, one_div, ← Real.exp_neg]
    congr 1; ring
  rw [hs, ← exp_add, exp_le_exp]
  have := neg_log_one_sub_le (x := 2 * y) (by linarith) (by linarith)
  have e : (2 * y) ^ 2 / (2 * (1 - 2 * y)) = 2 * (y ^ 2 / (1 - 2 * y)) := by
    field_simp
  rw [e] at this
  linarith

/-- Lower-tail MGF bound: `e^{y}/√(1+2y) ≤ exp(y²)` for `0 ≤ y`. -/
theorem mgf_bound_lower {y : ℝ} (hy0 : 0 ≤ y) :
    exp (-(-y)) * (1 / √(1 - 2 * (-y))) ≤ exp (y ^ 2) := by
  have hpos : 0 < 1 - 2 * (-y) := by linarith
  have hs : 1 / √(1 - 2 * (-y)) = exp (-(1 / 2) * log (1 - 2 * (-y))) := by
    rw [Real.sqrt_eq_rpow, Real.rpow_def_of_pos hpos, one_div, ← Real.exp_neg]
    congr 1; ring
  rw [hs, ← exp_add, exp_le_exp]
  have := log_one_add_ge (u := 2 * y) (by linarith)
  have e : 1 - 2 * (-y) = 1 + 2 * y := by ring
  rw [e]
  nlinarith

variable {Ω ι : Type*} [MeasurableSpace Ω] {P : Measure Ω} [IsProbabilityMeasure P] [Fintype ι]

/-- MGF of the weighted chi-square sum, for any `s` with `s λ_j < 1/2` for all `j`. -/
theorem mgf_weighted_chisq_sum (X : ι → Ω → ℝ) (hmeas : ∀ j, Measurable (X j))
    (hind : iIndepFun X P) (hlaw : ∀ j, P.map (X j) = gaussianReal 0 1) (lam : ι → ℝ) {s : ℝ}
    (hs : ∀ j, s * lam j < 1 / 2) :
    Integrable (fun ω => exp (s * ∑ j, lam j * (X j ω ^ 2 - 1))) P ∧
    mgf (fun ω => ∑ j, lam j * (X j ω ^ 2 - 1)) P s =
      ∏ j, (exp (-(s * lam j)) * (1 / √(1 - 2 * (s * lam j)))) := by
  set Y : ι → Ω → ℝ := fun j ω => lam j * (X j ω ^ 2 - 1)
  have hYm : ∀ j, Measurable (Y j) := fun j => by
    simp only [Y]; exact measurable_const.mul (((hmeas j).pow_const 2).sub measurable_const)
  have hYind : iIndepFun Y P :=
    hind.comp (fun j x => lam j * (x ^ 2 - 1)) (fun j => by fun_prop)
  have hint1 : ∀ j, Integrable (fun ω => exp (s * Y j ω)) P := by
    intro j
    have hg : Integrable (fun x => exp (s * (lam j * (x ^ 2 - 1)))) (P.map (X j)) := by
      rw [hlaw j]
      refine Integrable.of_integral_ne_zero ?_
      rw [integral_exp_weighted_chisq (hs j)]
      have : 0 < 1 - 2 * (s * lam j) := by linarith [hs j]
      positivity
    exact (integrable_map_measure hg.aestronglyMeasurable (hmeas j).aemeasurable).mp hg
  have hmgf1 : ∀ j, mgf (Y j) P s = exp (-(s * lam j)) * (1 / √(1 - 2 * (s * lam j))) := by
    intro j
    rw [mgf, ← integral_exp_weighted_chisq (hs j), ← hlaw j,
      integral_map (hmeas j).aemeasurable (by fun_prop)]
  have hsum : (fun ω => ∑ j, lam j * (X j ω ^ 2 - 1)) = ∑ j, Y j := by
    funext ω; simp [Y, Finset.sum_apply]
  refine ⟨?_, ?_⟩
  · have := hYind.integrable_exp_mul_sum hYm (s := univ) (fun j _ => hint1 j)
    simpa [hsum, Finset.sum_apply] using this
  · rw [hsum, hYind.mgf_sum hYm]
    exact prod_congr rfl fun j _ => hmgf1 j

/-- **Laurent–Massart, upper tail.** -/
theorem laurent_massart_upper (X : ι → Ω → ℝ) (hmeas : ∀ j, Measurable (X j))
    (hind : iIndepFun X P) (hlaw : ∀ j, P.map (X j) = gaussianReal 0 1) (lam : ι → ℝ)
    (hlam : ∀ j, 0 ≤ lam j) {a : ℝ} (ha : ∀ j, lam j ≤ a) (hv : 0 < ∑ j, lam j ^ 2)
    {τ : ℝ} (hτ : 0 < τ) :
    P.real {ω | 2 * √((∑ j, lam j ^ 2) * τ) + 2 * a * τ ≤ ∑ j, lam j * (X j ω ^ 2 - 1)} ≤
      exp (-τ) := by
  set v := ∑ j, lam j ^ 2
  have ha0 : 0 ≤ a := by
    obtain ⟨j, -, hj⟩ : ∃ j ∈ univ, 0 < lam j ^ 2 := by
      by_contra h; push Not at h
      exact absurd hv (not_lt.mpr (sum_nonpos fun j hj => h j hj))
    exact (hlam j).trans (ha j)
  set q := √v
  set w := √τ
  have hq : 0 < q := Real.sqrt_pos.mpr hv
  have hw : 0 < w := Real.sqrt_pos.mpr hτ
  have hden : 0 < q + 2 * a * w := by positivity
  set s := w / (q + 2 * a * w)
  have hs0 : 0 ≤ s := by positivity
  have h2as : 1 - 2 * a * s = q / (q + 2 * a * w) := by
    simp only [s]; field_simp; ring
  have h2as_pos : 0 < 1 - 2 * a * s := by rw [h2as]; positivity
  have hsl : ∀ j, s * lam j < 1 / 2 := by
    intro j; nlinarith [mul_le_mul_of_nonneg_left (ha j) hs0]
  obtain ⟨hint, hmgf⟩ := mgf_weighted_chisq_sum X hmeas hind hlaw lam hsl
  have hch := measure_ge_le_exp_mul_mgf (X := fun ω => ∑ j, lam j * (X j ω ^ 2 - 1))
    (2 * √(v * τ) + 2 * a * τ) hs0 hint
  refine hch.trans ?_
  rw [hmgf]
  -- each factor ≤ exp(s²λ²/(1−2as))
  have hfac : ∀ j, exp (-(s * lam j)) * (1 / √(1 - 2 * (s * lam j))) ≤
      exp ((s * lam j) ^ 2 / (1 - 2 * a * s)) := by
    intro j
    refine (mgf_bound_upper (mul_nonneg hs0 (hlam j)) (hsl j)).trans (exp_le_exp.mpr ?_)
    have h1 : 1 - 2 * a * s ≤ 1 - 2 * (s * lam j) := by
      nlinarith [mul_le_mul_of_nonneg_left (ha j) hs0]
    exact div_le_div_of_nonneg_left (sq_nonneg _) h2as_pos h1
  have hprod : ∏ j, (exp (-(s * lam j)) * (1 / √(1 - 2 * (s * lam j)))) ≤
      exp (s ^ 2 * v / (1 - 2 * a * s)) := by
    calc ∏ j, (exp (-(s * lam j)) * (1 / √(1 - 2 * (s * lam j))))
        ≤ ∏ j, exp ((s * lam j) ^ 2 / (1 - 2 * a * s)) :=
          prod_le_prod₀ (fun j _ => by
            have : 0 < 1 - 2 * (s * lam j) := by linarith [hsl j]
            positivity) fun j _ => hfac j
      _ = exp (s ^ 2 * v / (1 - 2 * a * s)) := by
          rw [← exp_sum]; congr 1
          rw [← sum_div]; congr 1
          simp only [v, mul_sum, mul_pow]
  calc exp (-s * (2 * √(v * τ) + 2 * a * τ)) *
        ∏ j, (exp (-(s * lam j)) * (1 / √(1 - 2 * (s * lam j))))
      ≤ exp (-s * (2 * √(v * τ) + 2 * a * τ)) * exp (s ^ 2 * v / (1 - 2 * a * s)) :=
        mul_le_mul_of_nonneg_left hprod (exp_pos _).le
    _ = exp (-τ) := by
        rw [← exp_add]; congr 1
        have hvt : √(v * τ) = q * w := Real.sqrt_mul hv.le τ
        have hqq : q ^ 2 = v := Real.sq_sqrt hv.le
        have hww : w ^ 2 = τ := Real.sq_sqrt hτ.le
        rw [hvt, h2as, ← hqq, ← hww]
        simp only [s]
        field_simp
        ring

/-- **Laurent–Massart, lower tail.** -/
theorem laurent_massart_lower (X : ι → Ω → ℝ) (hmeas : ∀ j, Measurable (X j))
    (hind : iIndepFun X P) (hlaw : ∀ j, P.map (X j) = gaussianReal 0 1) (lam : ι → ℝ)
    (hlam : ∀ j, 0 ≤ lam j) (hv : 0 < ∑ j, lam j ^ 2) {τ : ℝ} (hτ : 0 < τ) :
    P.real {ω | ∑ j, lam j * (X j ω ^ 2 - 1) ≤ -(2 * √((∑ j, lam j ^ 2) * τ))} ≤
      exp (-τ) := by
  set v := ∑ j, lam j ^ 2
  set q := √v
  set w := √τ
  have hq : 0 < q := Real.sqrt_pos.mpr hv
  have hw : 0 < w := Real.sqrt_pos.mpr hτ
  set s := w / q
  have hs0 : 0 ≤ s := by positivity
  have hsl : ∀ j, (-s) * lam j < 1 / 2 := by
    intro j; nlinarith [mul_nonneg hs0 (hlam j)]
  obtain ⟨hint, hmgf⟩ := mgf_weighted_chisq_sum X hmeas hind hlaw lam hsl
  have hch := measure_le_le_exp_mul_mgf (X := fun ω => ∑ j, lam j * (X j ω ^ 2 - 1))
    (-(2 * √(v * τ))) (neg_nonpos.mpr hs0) hint
  refine hch.trans ?_
  rw [hmgf]
  have hfac : ∀ j, exp (-(-s * lam j)) * (1 / √(1 - 2 * (-s * lam j))) ≤
      exp ((s * lam j) ^ 2) := by
    intro j
    have := mgf_bound_lower (y := s * lam j) (mul_nonneg hs0 (hlam j))
    simpa [neg_mul] using this
  have hprod : ∏ j, (exp (-(-s * lam j)) * (1 / √(1 - 2 * (-s * lam j)))) ≤
      exp (s ^ 2 * v) := by
    calc ∏ j, (exp (-(-s * lam j)) * (1 / √(1 - 2 * (-s * lam j))))
        ≤ ∏ j, exp ((s * lam j) ^ 2) :=
          prod_le_prod₀ (fun j _ => by
            have : 0 < 1 - 2 * (-s * lam j) := by linarith [hsl j]
            positivity) fun j _ => hfac j
      _ = exp (s ^ 2 * v) := by
          rw [← exp_sum]; congr 1
          simp only [v, mul_sum, mul_pow]
  calc exp (-(-s) * -(2 * √(v * τ))) *
        ∏ j, (exp (-(-s * lam j)) * (1 / √(1 - 2 * (-s * lam j))))
      ≤ exp (-(-s) * -(2 * √(v * τ))) * exp (s ^ 2 * v) :=
        mul_le_mul_of_nonneg_left hprod (exp_pos _).le
    _ = exp (-τ) := by
        rw [← exp_add]; congr 1
        have hvt : √(v * τ) = q * w := Real.sqrt_mul hv.le τ
        have hqq : q ^ 2 = v := Real.sq_sqrt hv.le
        have hww : w ^ 2 = τ := Real.sq_sqrt hτ.le
        rw [hvt, ← hqq, ← hww]
        simp only [s]
        field_simp
        ring

end TTMT
