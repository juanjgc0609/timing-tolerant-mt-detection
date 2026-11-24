/-
Controls for the second set of formal results (rule 40): radial maps, the Euclidean
reading, Lemma AA.9.1, and the mean of the M6 null gain. Negative controls drop ONE named
hypothesis and prove the conclusion fails on an explicit counterexample. Positive controls
are computed directly, or instantiate a theorem on a concrete configuration to show its
hypotheses can all hold at once (the theorem is not vacuous).

Hypotheses with a negative control:
- `radial_lipschitz`: `r(t) ≤ t` (`r(t) = 2t` doubles distances); `r` 1-Lipschitz (a step
  radius with `0 ≤ r(t) ≤ t` jumps).
- `norm_adjoint_le`: `‖Q y‖ ≤ ‖y‖` (`Q = 2·id` doubles norms).
-/
import TTMT.CappedScore
import TTMT.M6NullMean
import Mathlib.MeasureTheory.Integral.Bochner.Basic

open Matrix MeasureTheory Finset Real

namespace TTMT.Controls2

/-- Negative control, `r(t) ≤ t` dropped: with `r(t) = 2t` the radial map doubles distances. -/
theorem radial_needs_le_self :
    ¬ (‖radial (fun t => 2 * t) (1 : ℝ) - radial (fun t => 2 * t) 0‖ ≤ ‖(1 : ℝ) - 0‖) := by
  simp [radial]

/-- Negative control, Lipschitz radius dropped: the step `r(t) = 1` for `t ≥ 1`, `0` below,
satisfies `0 ≤ r(t) ≤ t` but maps `1` and `0.9` (distance 0.1) to `1` and `0` (distance 1). -/
theorem radial_needs_lipschitz :
    ¬ (‖radial (fun t : ℝ => if 1 ≤ t then 1 else 0) (1 : ℝ) -
        radial (fun t : ℝ => if 1 ≤ t then 1 else 0) (0.9 : ℝ)‖ ≤ ‖(1 : ℝ) - 0.9‖) := by
  norm_num [radial]

/-- Positive control: the cap at `c = 1` sends `3 ↦ 1` and `−1/2 ↦ −1/2` (computed), and the
distance `3/2` is at most `7/2`. -/
example : radial (fun t => min t 1) (3 : ℝ) = 1 ∧
    radial (fun t => min t 1) (-1 / 2 : ℝ) = -1 / 2 := by
  constructor <;> · norm_num [radial]

example : ‖radial (fun t => min t 1) (3 : ℝ) - radial (fun t => min t 1) (-1 / 2 : ℝ)‖ ≤
    ‖(3 : ℝ) - (-1 / 2)‖ := radial_lipschitz (goodRadius_cap zero_le_one) _ _

/-- Positive control for the logarithmic radius: `ℓ_c(0) = 0` (computed) and the theorem's
three conditions hold for `c = 6` (LOG6). -/
example : logRadius 6 0 = 0 := by simp [logRadius]

example : GoodRadius (logRadius 6) := goodRadius_log (by norm_num)

/-- Negative control for Lemma AA.9.1, `‖Q y‖ ≤ ‖y‖` dropped: `Q = 2·id` on `ℝ` has
`Q† = 2·id`, and `‖Q† 1‖ = 2 > 1`. -/
theorem adjoint_needs_contraction :
    ¬ (‖ContinuousLinearMap.adjoint ((2 : ℝ) • ContinuousLinearMap.id ℝ ℝ) (1 : ℝ)‖ ≤ ‖(1 : ℝ)‖) := by
  rw [LinearIsometryEquiv.map_smulₛₗ, ContinuousLinearMap.adjoint_id]
  simp

/-- Positive control for orthonormal-or-zero columns: columns `e₀` and `0` in `ℝ²`. -/
example (y : EuclideanSpace ℝ (Fin 2)) :
    ‖∑ j, y j • (![EuclideanSpace.single 0 1, 0] j : EuclideanSpace ℝ (Fin 2))‖ ≤ ‖y‖ :=
  norm_sum_smul_le _ (by
    intro i j hij
    fin_cases i <;> fin_cases j <;> simp_all) (by
    intro j; fin_cases j <;> simp) y

/-- Positive control, non-vacuity of `tailFunction_lipschitz`: two readings of two blocks,
each row reading one distinct noise sample with weight 1, the cap at 6, and `Q = 0`. -/
example : LipschitzWith 1 (tailFunction (σ := Fin 2) (κ := Fin 2) (ι := Fin 2 × (Fin 2 × Fin 3))
    (fun s r => {(s, r)}) (fun _ _ _ => 1) (fun t => min t 6)
    (0 : EuclideanSpace ℝ (Fin 6) →L[ℝ] R (Fin 2))) :=
  tailFunction_lipschitz _ (fun s => by
      intro a _ b _ hab
      simp only [Function.onFun, Finset.disjoint_singleton]
      intro h; exact hab (Prod.ext_iff.mp h).2 |>.elim)
    _ (fun _ _ => by simp) (goodRadius_cap (by norm_num)) _ (fun y => by simp)

/-- Positive control for `integral_quadForm`: a deterministic noise `ε = e` (Dirac law on a
point) gives `E[εᵀAε] = eᵀAe`, and with `C = e eᵀ = [[1,2],[2,4]]` the theorem's `tr(A Cᵀ)` is the same
number (here `e = (1, 2)`, `A = [[1, 0], [0, 3]]`: both sides are `13`). -/
example : trace (!![(1 : ℝ), 0; 0, 3] * !![(1 : ℝ), 2; 2, 4]ᵀ) = 13 ∧
    ![(1 : ℝ), 2] ⬝ᵥ (!![(1 : ℝ), 0; 0, 3] *ᵥ ![1, 2]) = 13 := by
  constructor
  · simp [trace, Fin.sum_univ_two, vecMul, dotProduct]; norm_num
  · simp [dotProduct, mulVec, Fin.sum_univ_two]; norm_num

example : ∫ _ω, ![(1 : ℝ), 2] ⬝ᵥ (!![(1 : ℝ), 0; 0, 3] *ᵥ ![1, 2]) ∂(Measure.dirac ()) =
    trace (!![(1 : ℝ), 0; 0, 3] * !![(1 : ℝ), 2; 2, 4]ᵀ) :=
  integral_quadForm (Measure.dirac ()) (fun _ => ![1, 2]) _ _
    (fun _ _ => integrable_const _) (fun i j => by fin_cases i <;> fin_cases j <;> norm_num)

/-- Positive control for `trace_mul_sum_vecMulVec`: one basis vector `e₀` gives `C₀₀`. -/
example (C : Matrix (Fin 2) (Fin 2) ℝ) :
    trace (C * ∑ _k : Fin 1, vecMulVec ![1, 0] ![1, 0]) = C 0 0 := by
  simp [trace, Matrix.mul_apply, Fin.sum_univ_two]

end TTMT.Controls2
