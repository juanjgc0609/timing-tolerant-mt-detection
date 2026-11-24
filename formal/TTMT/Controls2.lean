/-
Controls for the second set of formal results (rule 40): radial maps, the Euclidean
reading, Lemma AA.9.1, and the mean of the M6 null gain. Negative controls drop ONE named
hypothesis and prove the conclusion fails on an explicit counterexample. Positive controls
are computed directly, or instantiate a theorem on a concrete configuration to show its
hypotheses can all hold at once (the theorem is not vacuous).

Hypotheses with a negative control:
- `radial_lipschitz`: `r(t) ≤ t` (`r(t) = t + 1`, nonnegative and 1-Lipschitz); `r`
  1-Lipschitz (a step radius with `0 ≤ r(t) ≤ t` jumps).
- `influence_le` (AA.9.2): "differ only in block b" (two blocks differ).
- `cap_score_le` (AA.9.3) and `norm_adjoint_le`: `‖Q y‖ ≤ ‖y‖` (`Q = 2·id`).
- `cap_contamination` (AA.9.5(c)): "equal outside D" (D empty, blocks differ).
- `integral_gain_null` (E.2 mean): signal in the pair class.
- `tailFunction_lipschitz`: exempt; each of its hypotheses already has a negative control in
  the lemma it feeds (disjoint supports and row norm: `reading_contraction`; radius:
  `radial_lipschitz`; `‖Q y‖ ≤ ‖y‖`: `norm_adjoint_le`).
-/
import TTMT.CappedScore
import TTMT.M6NullMean
import TTMT.Superposition
import Mathlib.MeasureTheory.Integral.Bochner.Basic

open Matrix MeasureTheory Finset Real

namespace TTMT.Controls2

/-- Negative control, only `r(t) ≤ t` dropped: `r(t) = t + 1` is nonnegative and 1-Lipschitz,
but `φ(1) = 2`, `φ(−1) = −2`: distance 4 against 2. -/
theorem radial_needs_le_self :
    ¬ (‖radial (fun t => t + 1) (1 : ℝ) - radial (fun t => t + 1) (-1)‖ ≤ ‖(1 : ℝ) - (-1)‖) := by
  norm_num [radial]

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

/-- Positive control, non-vacuity of `tailFunction_lipschitz` with a nonzero `Q` (`Q = id`):
two readings of two blocks, each row reading one distinct noise sample with weight 1, and the
cap at 6. -/
example : LipschitzWith 1 (tailFunction (σ := Fin 2) (κ := Fin 2) (ι := Fin 2 × (Fin 2 × Fin 3))
    (fun s r => {(s, r)}) (fun _ _ _ => 1) (fun t => min t 6)
    (ContinuousLinearMap.id ℝ (R (Fin 2)))) :=
  tailFunction_lipschitz _ (fun s => by
      intro a _ b _ hab
      simp only [Function.onFun, Finset.disjoint_singleton]
      intro h; exact hab (Prod.ext_iff.mp h).2 |>.elim)
    _ (fun _ _ => by simp) (goodRadius_cap (by norm_num)) _ (fun y => le_rfl)

/-! Concrete read vectors: every block equal to `a e₁`. -/

/-- The read vector with `n` blocks, each `a e₁`. -/
noncomputable def blk (κ : Type*) [Fintype κ] (a : ℝ) : R κ :=
  WithLp.toLp 2 fun _ => PiLp.single 2 (0 : Fin 3) a

theorem norm_blockRadial_blk {κ : Type*} [Fintype κ] {r : ℝ → ℝ} (hr : GoodRadius r) (a : ℝ) :
    ‖blockRadial r (blk κ a)‖ = √(Fintype.card κ * r |a| ^ 2) := by
  have h := norm_blockRadial_sq hr (blk κ a)
  simp only [blk, PiLp.toLp_apply, PiLp.norm_single, Real.norm_eq_abs,
    Finset.sum_const, Finset.card_univ, nsmul_eq_mul] at h
  rw [← h, Real.sqrt_sq (norm_nonneg _)]
  rfl

theorem adjoint_id_apply {κ : Type*} [Fintype κ] (v : R κ) :
    ContinuousLinearMap.adjoint (ContinuousLinearMap.id ℝ (R κ)) v = v := by
  rw [ContinuousLinearMap.adjoint_id]; rfl

/-- Positive control for AA.9.2 (tight case, computed): one block `3e₁` against `0`, identity
radius: the score moves by exactly `3 = r(3) + r(0)`. -/
example : |‖ContinuousLinearMap.adjoint (ContinuousLinearMap.id ℝ (R (Fin 1)))
      (blockRadial id (blk (Fin 1) 3))‖ -
    ‖ContinuousLinearMap.adjoint (ContinuousLinearMap.id ℝ (R (Fin 1)))
      (blockRadial id (blk (Fin 1) 0))‖| = 3 := by
  rw [adjoint_id_apply, adjoint_id_apply, norm_blockRadial_blk goodRadius_id,
    norm_blockRadial_blk goodRadius_id]
  norm_num
  rw [abs_of_nonneg (Real.sqrt_nonneg _), show (9 : ℝ) = 3 ^ 2 by norm_num,
    Real.sqrt_sq (by norm_num)]

/-- Negative control for AA.9.2, "differ only in block b" dropped: two blocks `3e₁` against
`0`; the score moves by `√18 > 3 = r(‖v_0‖) + r(‖v'_0‖)`. -/
theorem influence_needs_one_block :
    ¬ (|‖ContinuousLinearMap.adjoint (ContinuousLinearMap.id ℝ (R (Fin 2)))
        (blockRadial id (blk (Fin 2) 3))‖ -
      ‖ContinuousLinearMap.adjoint (ContinuousLinearMap.id ℝ (R (Fin 2)))
        (blockRadial id (blk (Fin 2) 0))‖| ≤ id ‖(blk (Fin 2) 3) 0‖ + id ‖(blk (Fin 2) 0) 0‖) := by
  rw [adjoint_id_apply, adjoint_id_apply, norm_blockRadial_blk goodRadius_id,
    norm_blockRadial_blk goodRadius_id]
  simp only [blk, PiLp.toLp_apply, PiLp.norm_single, id]
  norm_num
  rw [abs_of_nonneg (Real.sqrt_nonneg _), Real.lt_sqrt (by norm_num)]
  norm_num

/-- Positive control for AA.9.3 (tight, computed): two saturated blocks `5e₁`, cap 1, `Q = id`:
`σ² = 2 = n_b c²`. -/
example : ‖ContinuousLinearMap.adjoint (ContinuousLinearMap.id ℝ (R (Fin 2)))
    (blockRadial (fun t => min t 1) (blk (Fin 2) 5))‖ ^ 2 = Fintype.card (Fin 2) * 1 ^ 2 := by
  rw [adjoint_id_apply, norm_blockRadial_blk (goodRadius_cap zero_le_one),
    Real.sq_sqrt (by positivity)]
  norm_num

/-- Negative control for AA.9.3, `‖Q y‖ ≤ ‖y‖` dropped: `Q = 2·id` gives `σ² = 8 > 2`. -/
theorem cap_score_needs_contraction :
    ¬ (‖ContinuousLinearMap.adjoint ((2 : ℝ) • ContinuousLinearMap.id ℝ (R (Fin 2)))
      (blockRadial (fun t => min t 1) (blk (Fin 2) 5))‖ ^ 2 ≤ Fintype.card (Fin 2) * 1 ^ 2) := by
  rw [LinearIsometryEquiv.map_smulₛₗ, ContinuousLinearMap.adjoint_id]
  rw [_root_.smul_apply, starRingEnd_apply, star_trivial,
    ContinuousLinearMap.id_apply, norm_smul, Real.norm_eq_abs, abs_two, norm_blockRadial_blk (goodRadius_cap zero_le_one), mul_pow, Real.sq_sqrt (by positivity)]
  norm_num

/-- Positive control for AA.9.5(c) (computed): one contaminated block `5e₁` against `0`, cap 1:
the score moves by `1 ≤ 2c√1 = 2`. -/
example : |‖ContinuousLinearMap.adjoint (ContinuousLinearMap.id ℝ (R (Fin 1)))
      (blockRadial (fun t => min t 1) (blk (Fin 1) 5))‖ -
    ‖ContinuousLinearMap.adjoint (ContinuousLinearMap.id ℝ (R (Fin 1)))
      (blockRadial (fun t => min t 1) (blk (Fin 1) 0))‖| = 1 ∧
    2 * 1 * √((Finset.univ : Finset (Fin 1)).card) = 2 := by
  rw [adjoint_id_apply, adjoint_id_apply, norm_blockRadial_blk (goodRadius_cap zero_le_one),
    norm_blockRadial_blk (goodRadius_cap zero_le_one)]
  norm_num

/-- Negative control for AA.9.5(c), "equal outside D" dropped: `D = ∅` but the blocks differ;
the bound would be `0` and the score moves by `1`. -/
theorem contamination_needs_equal_outside :
    ¬ (|‖ContinuousLinearMap.adjoint (ContinuousLinearMap.id ℝ (R (Fin 1)))
        (blockRadial (fun t => min t 1) (blk (Fin 1) 5))‖ -
      ‖ContinuousLinearMap.adjoint (ContinuousLinearMap.id ℝ (R (Fin 1)))
        (blockRadial (fun t => min t 1) (blk (Fin 1) 0))‖| ≤
      2 * 1 * √((∅ : Finset (Fin 1)).card)) := by
  rw [adjoint_id_apply, adjoint_id_apply, norm_blockRadial_blk (goodRadius_cap zero_le_one),
    norm_blockRadial_blk (goodRadius_cap zero_le_one)]
  norm_num

/-- Positive control for AB.A.2 margin: `G θ m = θ m` on `ℝ`; the coincident pair `2·1 + 2·3`
lies at distance 0 from the class. -/
example : Metric.infDist ((LinearMap.lsmul ℝ ℝ 2) 1 + (LinearMap.lsmul ℝ ℝ 2) 3)
    (oneEventClass (fun θ : ℝ => LinearMap.lsmul ℝ ℝ θ)) = 0 :=
  coincident_pair_margin_zero _ 2 1 3

/-- Positive control for the E.2 mean (computed): Dirac noise `ε = (1, 2)`, `P₁ = diag(1,0)`,
`P₂ = I`, `μ = (3, 0)` in both classes, `C = εεᵀ`. Both the integral of the gain and
`tr((P₂−P₁)Cᵀ)` equal `4`. -/
example : ∫ _ω, (rss !![(1 : ℝ), 0; 0, 0] (![3, 0] + ![1, 2]) -
      rss (1 : Matrix (Fin 2) (Fin 2) ℝ) (![3, 0] + ![1, 2])) ∂(Measure.dirac ()) = 4 ∧
    trace ((1 - !![(1 : ℝ), 0; 0, 0]) * !![(1 : ℝ), 2; 2, 4]ᵀ) = 4 := by
  constructor
  · simp [rss, dotProduct, mulVec, Fin.sum_univ_two]; norm_num
  · simp [trace, Fin.sum_univ_two, Matrix.mul_apply]

/-- Negative control for the E.2 mean, "signal in the pair class" dropped: `P₁ = I`, `P₂ = 0`,
`μ = 1`, `ε = 0` (so `C = 0`): the mean gain is `−1`, not `tr((P₂−P₁)Cᵀ) = 0`. -/
theorem mean_needs_signal_in_class :
    ∫ _ω, (rss (1 : Matrix (Fin 1) (Fin 1) ℝ) (![1] + 0) -
      rss (0 : Matrix (Fin 1) (Fin 1) ℝ) (![1] + 0)) ∂(Measure.dirac ()) ≠
    trace ((0 - 1 : Matrix (Fin 1) (Fin 1) ℝ) * (0 : Matrix (Fin 1) (Fin 1) ℝ)ᵀ) := by
  simp [rss, dotProduct, mulVec]

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
