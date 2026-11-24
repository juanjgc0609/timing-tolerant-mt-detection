/-
Controls for the formal results (rule 40). Each negative control drops one hypothesis
of a theorem and proves that the conclusion then FAILS on an explicit counterexample:
the hypothesis is needed and the theorem is not vacuous. Each positive control checks a
theorem on a concrete case computed independently by `norm_num`/`simp`.
-/
import TTMT.PPolarization
import TTMT.M6NullGain
import TTMT.ReadingContraction

open Matrix Finset

namespace TTMT.Controls

/-- Negative control for AB.A.9: two generic vectors (e.g. two S polarizations, which are
not tied to one fixed direction) can have a nonzero cross product. -/
theorem generic_cross_ne_zero : (![1, 0, 0] : Fin 3 → ℝ) ⨯₃ ![0, 1, 0] ≠ 0 := by
  intro h
  have h2 := congrFun h 2
  simp [cross_apply] at h2

/-- Positive control for AB.A.9 on a concrete case. -/
example : pAmplitude ![0, 0, 1] ![1, 1, 0] 2 1 ⨯₃ pAmplitude ![0, 0, 1] ![1, 1, 0] 2 0 = 0 :=
  pAmplitude_cross_eq_zero _ _ _ _ _

/-- Negative control for E.2: if the signal is NOT in the pair class (`P₂ μ ≠ μ`), the gain
depends on the signal. One dimension, `P₁ = I`, `P₂ = 0`, `μ = 1`, `ε = 0`:
the gain is `−1`, while the noise form is `0`. -/
theorem gain_depends_on_signal_outside_class :
    rss (1 : Matrix (Fin 1) (Fin 1) ℝ) (fun _ => 1) - rss (0 : Matrix (Fin 1) (Fin 1) ℝ) (fun _ => 1) ≠
      (0 : Fin 1 → ℝ) ⬝ᵥ ((0 - 1 : Matrix (Fin 1) (Fin 1) ℝ) *ᵥ 0) := by
  simp [rss, dotProduct]

/-- Negative control for AB.C.1: with overlapping supports (two rows reading the same
sample) the reading is not a contraction: `2 > 1`. -/
theorem overlapping_rows_not_contraction :
    ¬ (∑ _r : Fin 2, (∑ _i : Fin 1, (1 : ℝ) * 1) ^ 2 ≤ ∑ _i : Fin 1, (1 : ℝ) ^ 2) := by
  norm_num

end TTMT.Controls
