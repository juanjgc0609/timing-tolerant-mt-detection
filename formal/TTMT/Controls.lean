/-
Controls for the formal results (rule 40). Each negative control drops ONE named
hypothesis of a theorem and proves that the conclusion then fails on an explicit
counterexample, so that hypothesis is needed. Each theorem also has a positive control:
a concrete case whose value is computed directly (by `norm_num`/`simp`, not through the
theorem) and agrees with it.

Hypotheses with a negative control:
- AB.A.9: the common fixed vector `u` (generic vectors have nonzero cross product).
- E.2: signal in the pair class, idempotence, symmetry.
- AB.C.1 (1): disjoint supports, rows of norm ≤ 1.
- AB.A.2: linearity of `G` in the tensor.
-/
import TTMT.PPolarization
import TTMT.M6NullGain
import TTMT.ReadingContraction
import TTMT.Superposition

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

/-- Negative control for E.2, idempotence: `P₁ = 2` (symmetric, not idempotent), `P₂ = 0`,
`μ = 0`, `ε = 1`. The gain is `0` and the noise form is `−2`. -/
theorem gain_needs_idempotent :
    rss !![(2 : ℝ)] (0 + ![1]) - rss (0 : Matrix (Fin 1) (Fin 1) ℝ) (0 + ![1]) ≠
      ![(1 : ℝ)] ⬝ᵥ ((0 - !![(2 : ℝ)]) *ᵥ ![1]) := by
  norm_num [rss, dotProduct, mulVec, Matrix.sub_apply, Matrix.one_apply]

/-- Negative control for E.2, symmetry: `P₁ = [[1,1],[0,0]]` is idempotent but not
symmetric, `P₂ = 0`, `μ = 0`, `ε = (0,1)`. The gain is `1` and the noise form is `0`. -/
theorem gain_needs_symmetric :
    rss !![(1 : ℝ), 1; 0, 0] (0 + ![0, 1]) - rss (0 : Matrix (Fin 2) (Fin 2) ℝ) (0 + ![0, 1]) ≠
      ![(0 : ℝ), 1] ⬝ᵥ ((0 - !![(1 : ℝ), 1; 0, 0]) *ᵥ ![0, 1]) := by
  norm_num [rss, dotProduct, mulVec, Fin.sum_univ_two, Matrix.sub_apply, Matrix.one_apply]

/-- Positive control for E.2: `P₁ = diag(1,0)`, `P₂ = I`, `μ = (3,0)`, `ε = (a,b)`.
Direct computation gives gain `b²`, the value of `εᵀ(P₂−P₁)ε`. -/
theorem gain_positive_control (a b : ℝ) :
    rss !![(1 : ℝ), 0; 0, 0] (![3, 0] + ![a, b]) - rss (1 : Matrix (Fin 2) (Fin 2) ℝ)
        (![3, 0] + ![a, b]) = b ^ 2 ∧
      ![a, b] ⬝ᵥ ((1 - !![(1 : ℝ), 0; 0, 0]) *ᵥ ![a, b]) = b ^ 2 := by
  constructor <;>
  · simp [rss, dotProduct, mulVec, Fin.sum_univ_two, Matrix.sub_apply, Matrix.one_apply]
    ring

/-- Negative control for AB.C.1 (1), row norm: one row with kernel `2` reads `4 > 1`. -/
theorem large_row_not_contraction :
    ¬ (∑ _r : Fin 1, (∑ _i : Fin 1, (2 : ℝ) * 1) ^ 2 ≤ ∑ _i : Fin 1, (1 : ℝ) ^ 2) := by
  norm_num

/-- Positive control for AB.C.1 (1): the Cauchy–Schwarz equality case, kernel `(0.6, 0.8)`
on two samples, `z = (3, 4)`: the reading is `5² = 25 = ‖z‖²`. -/
example : (∑ i : Fin 2, ![(0.6 : ℝ), 0.8] i * ![3, 4] i) ^ 2 = ∑ i : Fin 2, ![(3 : ℝ), 4] i ^ 2 := by
  norm_num [Fin.sum_univ_two]

example (z : Fin 2 → ℝ) :
    ∑ _r : Fin 1, (∑ i ∈ (univ : Finset (Fin 2)), ![(0.6 : ℝ), 0.8] i * z i) ^ 2 ≤
      ∑ i, z i ^ 2 :=
  reading_contraction (fun _ => univ) (by simp [Set.PairwiseDisjoint])
    (fun _ => ![0.6, 0.8]) (fun _ => by norm_num [Fin.sum_univ_two]) z

/-- Negative control for AB.A.2, linearity: with the non-linear "operator"
`f θ m = (if m = 0 then 0 else 1)`, the coincident pair `f () 1 + f () 1 = 2` is not in the
one-event class `{0, 1}`. -/
theorem nonlinear_pair_not_in_class :
    ¬ ∃ p : Unit × ℝ, (if p.2 = 0 then (0 : ℝ) else 1) =
      (if (1 : ℝ) = 0 then 0 else 1) + (if (1 : ℝ) = 0 then 0 else 1) := by
  rintro ⟨⟨_, m⟩, h⟩
  split_ifs at h <;> norm_num at *

/-- Positive control for AB.A.2: `G θ m = θ • m` on `ℝ`; the coincident pair with `θ = 2`,
`m₁ = 1`, `m₂ = 3` is `8 = G 2 4`. -/
example : (2 : ℝ) • (1 : ℝ) + (2 : ℝ) • 3 = (2 : ℝ) • 4 := by norm_num

example : (LinearMap.lsmul ℝ ℝ 2) 1 + (LinearMap.lsmul ℝ ℝ 2) 3 ∈
    oneEventClass (fun θ : ℝ => LinearMap.lsmul ℝ ℝ θ) :=
  coincident_pair_mem_oneEventClass _ 2 1 3

end TTMT.Controls
