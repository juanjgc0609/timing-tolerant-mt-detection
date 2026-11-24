/-
Deterministic core of Proposition AB.C.1 (docs/proposal/Supplementary_Evidence_Appendix_AB.md,
§AB.C.1): the constant κ = 2.146 of AA.9.5(b) can be replaced by 1.

(1) `reading_contraction`: reading one sample per block and component is a contraction. If
    each row r of the reading is a kernel `k r` of norm ≤ 1 supported on a set `S r` of
    white-noise samples, and rows have pairwise disjoint supports (different components,
    different stations, or P and S of one station separated by more than the kernel
    support: the independence assumption of AA.9.5), then `∑_r ⟨k_r, z⟩² ≤ ‖z‖²`.
(2) `lipschitz_comp_one`: a composition of 1-Lipschitz maps is 1-Lipschitz;
    `lipschitz_sup`: the maximum over a finite set of readings of 1-Lipschitz functions is
    1-Lipschitz.

NOT formalized here: (a) that the radial map (cap, logarithmic) and `Q_hᵀ` are 1-Lipschitz
(Lemma AA.9.1); it enters as a hypothesis of `lipschitz_comp_one`; (b) Gaussian concentration
P(F ≥ E F + t) ≤ e^{-t²/2} and the bound E F ≤ √6 + √(2 n_b log L), which Mathlib lacks;
they are cited results.
-/
import Mathlib.Algebra.Order.BigOperators.Ring.Finset
import Mathlib.Topology.MetricSpace.Lipschitz

open Finset

namespace TTMT

/-- **AB.C.1 (1).** A reading with rows of norm ≤ 1 and disjoint supports is a contraction. -/
theorem reading_contraction {ι ρ : Type*} [Fintype ι] [Fintype ρ] [DecidableEq ι]
    (S : ρ → Finset ι) (hS : Set.PairwiseDisjoint (↑(univ : Finset ρ)) S)
    (k : ρ → ι → ℝ) (hk : ∀ r, ∑ i ∈ S r, k r i ^ 2 ≤ 1) (z : ι → ℝ) :
    ∑ r, (∑ i ∈ S r, k r i * z i) ^ 2 ≤ ∑ i, z i ^ 2 := by
  calc ∑ r, (∑ i ∈ S r, k r i * z i) ^ 2
      ≤ ∑ r, ∑ i ∈ S r, z i ^ 2 := by
        refine sum_le_sum fun r _ => ?_
        calc (∑ i ∈ S r, k r i * z i) ^ 2
            ≤ (∑ i ∈ S r, k r i ^ 2) * ∑ i ∈ S r, z i ^ 2 := sum_mul_sq_le_sq_mul_sq _ _ _
          _ ≤ 1 * ∑ i ∈ S r, z i ^ 2 :=
            mul_le_mul_of_nonneg_right (hk r) (sum_nonneg fun _ _ => sq_nonneg _)
          _ = ∑ i ∈ S r, z i ^ 2 := one_mul _
    _ = ∑ i ∈ univ.biUnion S, z i ^ 2 := (sum_biUnion hS).symm
    _ ≤ ∑ i, z i ^ 2 :=
        sum_le_sum_of_subset_of_nonneg (subset_univ _) fun _ _ _ => sq_nonneg _

/-- **AB.C.1 (2a).** The composition of three 1-Lipschitz maps (reading, radial map,
`Q_hᵀ` followed by the norm) is 1-Lipschitz. -/
theorem lipschitz_comp_one {α β γ δ : Type*} [PseudoEMetricSpace α] [PseudoEMetricSpace β]
    [PseudoEMetricSpace γ] [PseudoEMetricSpace δ] {A : α → β} {φ : β → γ} {g : γ → δ}
    (hA : LipschitzWith 1 A) (hφ : LipschitzWith 1 φ) (hg : LipschitzWith 1 g) :
    LipschitzWith 1 (g ∘ φ ∘ A) := by
  simpa using hg.comp (hφ.comp hA)

/-- **AB.C.1 (2b).** The maximum of a nonempty finite family of real 1-Lipschitz functions
is 1-Lipschitz. -/
theorem lipschitz_sup {α σ : Type*} [PseudoMetricSpace α] (s : Finset σ) (hs : s.Nonempty)
    (f : σ → α → ℝ) (hf : ∀ j ∈ s, LipschitzWith 1 (f j)) :
    LipschitzWith 1 fun x => s.sup' hs fun j => f j x := by
  refine LipschitzWith.of_le_add fun x y => ?_
  refine sup'_le hs _ fun j hj => ?_
  have h := (hf j hj).dist_le_mul x y
  have h1 : f j x ≤ f j y + dist x y := by
    have := (le_abs_self (f j x - f j y)).trans (by simpa [Real.dist_eq] using h)
    linarith
  have h2 := le_sup' (fun j => f j y) hj
  linarith

end TTMT
