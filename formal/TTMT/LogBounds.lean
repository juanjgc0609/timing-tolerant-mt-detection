/-
Two elementary logarithm bounds used by the Laurent–Massart tail (`LaurentMassart.lean`):
- `neg_log_one_sub_le`: `−log(1−x) − x ≤ x² / (2(1−x))` for `0 ≤ x < 1`;
- `log_one_add_ge`: `u − u²/2 ≤ log(1+u)` for `0 ≤ u`.
Both by monotonicity: the difference vanishes at 0 and has a nonnegative derivative.
-/
import Mathlib.Analysis.SpecialFunctions.Log.Deriv
import Mathlib.Analysis.Calculus.Deriv.MeanValue

open Real Set

namespace TTMT

theorem neg_log_one_sub_le {x : ℝ} (hx0 : 0 ≤ x) (hx1 : x < 1) :
    -log (1 - x) - x ≤ x ^ 2 / (2 * (1 - x)) := by
  let h : ℝ → ℝ := fun y => y ^ 2 / (2 * (1 - y)) + y + log (1 - y)
  have hderiv : ∀ y, y < 1 → HasDerivAt h (y ^ 2 / (2 * (1 - y) ^ 2)) y := by
    intro y hy
    have h1y : (1 - y) ≠ 0 := by linarith
    have hA : HasDerivAt (fun y : ℝ => y ^ 2 / (2 * (1 - y)))
        ((2 * y * (2 * (1 - y)) - y ^ 2 * (2 * (-1))) / (2 * (1 - y)) ^ 2) y := by
      have hn : HasDerivAt (fun y : ℝ => y ^ 2) (2 * y) y := by
        simpa using hasDerivAt_pow 2 y
      have hd : HasDerivAt (fun y : ℝ => 2 * (1 - y)) (2 * (-1)) y := by
        simpa using ((hasDerivAt_id y).const_sub 1).const_mul 2
      exact hn.div hd (by positivity)
    have hL : HasDerivAt (fun y : ℝ => log (1 - y)) ((-1) / (1 - y)) y := by
      simpa using ((hasDerivAt_id y).const_sub 1).log h1y
    refine ((hA.add (hasDerivAt_id y)).add hL).congr_deriv ?_
    have h2 : (2 * (1 - y)) ^ 2 = 4 * (1 - y) ^ 2 := by ring
    rw [eq_comm, h2, div_add_one (by positivity), div_add_div _ _ (by positivity) h1y,
      div_eq_div_iff (by positivity) (by positivity)]
    ring
  have hmono : MonotoneOn h (Ico 0 1) := by
    refine monotoneOn_of_deriv_nonneg (convex_Ico 0 1) ?_ ?_ ?_
    · exact fun y hy => (hderiv y hy.2).continuousAt.continuousWithinAt
    · rw [interior_Ico]
      exact fun y hy => (hderiv y hy.2).differentiableAt.differentiableWithinAt
    · rw [interior_Ico]
      intro y hy
      rw [(hderiv y hy.2).deriv]
      positivity
  have := hmono ⟨le_rfl, by norm_num⟩ ⟨hx0, hx1⟩ hx0
  simp only [h, sub_zero, log_one, add_zero] at this
  norm_num at this
  linarith

theorem log_one_add_ge {u : ℝ} (hu : 0 ≤ u) : u - u ^ 2 / 2 ≤ log (1 + u) := by
  let h : ℝ → ℝ := fun y => log (1 + y) - y + y ^ 2 / 2
  have hderiv : ∀ y, 0 ≤ y → HasDerivAt h (y ^ 2 / (1 + y)) y := by
    intro y hy
    have h1y : (1 + y) ≠ 0 := by linarith
    have hL : HasDerivAt (fun y : ℝ => log (1 + y)) (1 / (1 + y)) y := by
      simpa using ((hasDerivAt_id y).const_add 1).log h1y
    have hq : HasDerivAt (fun y : ℝ => y ^ 2 / 2) (2 * y / 2) y := by
      simpa using (hasDerivAt_pow 2 y).div_const 2
    refine ((hL.sub (hasDerivAt_id y)).add hq).congr_deriv ?_
    rw [sub_add_eq_add_sub, div_add_div _ _ h1y two_ne_zero,
      div_sub_one (mul_ne_zero h1y two_ne_zero), div_eq_div_iff (mul_ne_zero h1y two_ne_zero) h1y]
    ring
  have hmono : MonotoneOn h (Ici 0) := by
    refine monotoneOn_of_deriv_nonneg (convex_Ici 0) ?_ ?_ ?_
    · exact fun y hy => (hderiv y hy).continuousAt.continuousWithinAt
    · rw [interior_Ici]
      exact fun y hy => (hderiv y (le_of_lt hy)).differentiableAt.differentiableWithinAt
    · rw [interior_Ici]
      intro y hy
      rw [(hderiv y (le_of_lt hy)).deriv]
      have : 0 < 1 + y := by linarith [mem_Ioi.mp hy]
      positivity
  have := hmono (mem_Ici.mpr le_rfl) (mem_Ici.mpr hu) hu
  simp only [h, add_zero, log_one, sub_zero] at this
  norm_num at this
  linarith

end TTMT
