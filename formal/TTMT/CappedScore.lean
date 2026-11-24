/-
The deterministic chain behind AB.C.1 and AA.9.1-AA.9.5, in Euclidean spaces.

Spaces.
- White noise `z ∈ EuclideanSpace ℝ ι`.
- Read vector `v ∈ R κ := PiLp 2 (fun _ : κ => EuclideanSpace ℝ (Fin 3))`: one block in `ℝ³`
  per (station, phase), with the Euclidean norm `‖v‖² = ∑_k ‖v_k‖²` of `ℝ^{3 n_b}`.
- `Q : F →L[ℝ] R κ` has orthonormal or zero columns (`F` is the 6-dim radiation space);
  `Qᵀ` is its adjoint `Q†`. We only use `‖Q y‖ ≤ ‖y‖`, which is proved from the columns in
  `norm_sum_smul_le` below.

Results.
- `readingMap_lipschitz`: the bridge. The one-sample-per-block reading with rows of norm
  ≤ 1 on disjoint supports is 1-Lipschitz from the noise to the read vector, now in the
  Euclidean metric (the review of PR #77 noted that `ι → ℝ` carries the sup metric).
- `norm_adjoint_le` (Lemma AA.9.1): `‖Q† a‖ ≤ ‖a‖` and `|‖Q†a‖ − ‖Q†b‖| ≤ ‖a − b‖`.
- `blockRadial_lipschitz`: the block-wise radial map is 1-Lipschitz on `R κ`.
- `tailFunction_lipschitz` (AB.C.1, deterministic part complete): for every finite set of
  readings, `F(z) = max_s ‖Q† φ(A_s z)‖` is 1-Lipschitz in the noise, for the identity,
  the cap and the logarithmic map. `read_le_tailFunction`: the score at any data-dependent
  reading is at most `F`.
- `influence_le` (AA.9.2), `cap_score_le` (AA.9.3, first claim), `cap_contamination`
  (AA.9.5(c)).

The probabilistic steps (Gaussian concentration of `F`, `E F ≤ √6 + √(2 n_b log L)`) are
not formalized here.
-/
import TTMT.ReadingContraction
import TTMT.RadialMap
import Mathlib.Analysis.InnerProductSpace.PiL2
import Mathlib.Analysis.InnerProductSpace.Adjoint

open Finset Real

namespace TTMT

/-- Read space: `n_b` blocks in `ℝ³` with the Euclidean norm of `ℝ^{3 n_b}`. -/
abbrev R (κ : Type*) [Fintype κ] := PiLp 2 (fun _ : κ => EuclideanSpace ℝ (Fin 3))

section Bridge

variable {ι ρ : Type*} [Fintype ι] [Fintype ρ] [DecidableEq ι]

/-- The reading map `z ↦ (∑_{i ∈ S r} k r i z_i)_r` between Euclidean spaces. -/
noncomputable def readingMap (S : ρ → Finset ι) (k : ρ → ι → ℝ) (z : EuclideanSpace ℝ ι) :
    EuclideanSpace ℝ ρ :=
  WithLp.toLp 2 fun r => ∑ i ∈ S r, k r i * z i

/-- **Bridge for AB.C.1.** The reading is 1-Lipschitz in the Euclidean metric. -/
theorem readingMap_lipschitz (S : ρ → Finset ι)
    (hS : Set.PairwiseDisjoint (↑(univ : Finset ρ)) S) (k : ρ → ι → ℝ)
    (hk : ∀ r, ∑ i ∈ S r, k r i ^ 2 ≤ 1) : LipschitzWith 1 (readingMap S k) := by
  refine LipschitzWith.of_dist_le_mul fun x y => ?_
  rw [NNReal.coe_one, one_mul]
  have hsq : dist (readingMap S k x) (readingMap S k y) ^ 2 ≤ dist x y ^ 2 := by
    rw [EuclideanSpace.dist_eq, EuclideanSpace.dist_eq, Real.sq_sqrt (by positivity),
      Real.sq_sqrt (by positivity)]
    simp only [readingMap, PiLp.toLp_apply, Real.dist_eq, sq_abs]
    have := reading_contraction S hS k hk (fun i => x i - y i)
    simpa [← Finset.sum_sub_distrib, mul_sub] using this
  exact (pow_le_pow_iff_left₀ dist_nonneg dist_nonneg two_ne_zero).mp hsq

end Bridge

section Projection

variable {F G : Type*} [NormedAddCommGroup F] [InnerProductSpace ℝ F] [CompleteSpace F]
  [NormedAddCommGroup G] [InnerProductSpace ℝ G] [CompleteSpace G]

/-- **Lemma AA.9.1 (contraction), second claim.** If `‖Q y‖ ≤ ‖y‖` for all `y`, then
`‖Q† a‖ ≤ ‖a‖`. -/
theorem norm_adjoint_le (Q : F →L[ℝ] G) (hQ : ∀ y, ‖Q y‖ ≤ ‖y‖) (a : G) :
    ‖ContinuousLinearMap.adjoint Q a‖ ≤ ‖a‖ := by
  set b := ContinuousLinearMap.adjoint Q a
  have h : ‖b‖ ^ 2 ≤ ‖b‖ * ‖a‖ := by
    calc ‖b‖ ^ 2 = @inner ℝ F _ b b := (real_inner_self_eq_norm_sq b).symm
      _ = @inner ℝ G _ a (Q b) := ContinuousLinearMap.adjoint_inner_left Q b a
      _ ≤ ‖a‖ * ‖Q b‖ := real_inner_le_norm a (Q b)
      _ ≤ ‖a‖ * ‖b‖ := mul_le_mul_of_nonneg_left (hQ b) (norm_nonneg a)
      _ = ‖b‖ * ‖a‖ := mul_comm _ _
  rcases eq_or_lt_of_le (norm_nonneg b) with h0 | hpos
  · rw [← h0]; exact norm_nonneg a
  · nlinarith

/-- **Lemma AA.9.1, first claim.** `|‖Q†a‖ − ‖Q†b‖| ≤ ‖a − b‖`. -/
theorem norm_adjoint_sub_le (Q : F →L[ℝ] G) (hQ : ∀ y, ‖Q y‖ ≤ ‖y‖) (a b : G) :
    |‖ContinuousLinearMap.adjoint Q a‖ - ‖ContinuousLinearMap.adjoint Q b‖| ≤ ‖a - b‖ := by
  calc |‖ContinuousLinearMap.adjoint Q a‖ - ‖ContinuousLinearMap.adjoint Q b‖|
      ≤ ‖ContinuousLinearMap.adjoint Q a - ContinuousLinearMap.adjoint Q b‖ :=
        abs_norm_sub_norm_le _ _
    _ = ‖ContinuousLinearMap.adjoint Q (a - b)‖ := by rw [map_sub]
    _ ≤ ‖a - b‖ := norm_adjoint_le Q hQ _

/-- Orthonormal-or-zero columns give `‖Q y‖ ≤ ‖y‖`: for pairwise orthogonal `q_j` of norm
`≤ 1`, `‖∑_j y_j q_j‖² = ∑_j y_j² ‖q_j‖² ≤ ∑_j y_j²`. -/
theorem norm_sum_smul_le {m : Type*} [Fintype m] (q : m → G)
    (horth : Pairwise fun i j => @inner ℝ G _ (q i) (q j) = 0) (hq : ∀ j, ‖q j‖ ≤ 1)
    (y : EuclideanSpace ℝ m) : ‖∑ j, y j • q j‖ ≤ ‖y‖ := by
  classical
  have hsq : ‖∑ j, y j • q j‖ ^ 2 ≤ ‖y‖ ^ 2 := by
    rw [EuclideanSpace.norm_eq, Real.sq_sqrt (by positivity), ← real_inner_self_eq_norm_sq,
      sum_inner]
    simp_rw [inner_sum, real_inner_smul_left, real_inner_smul_right]
    have : ∀ i, ∑ j, y i * (y j * @inner ℝ G _ (q i) (q j)) = y i ^ 2 * ‖q i‖ ^ 2 := by
      intro i
      rw [Finset.sum_eq_single i (fun j _ hji => by rw [horth (Ne.symm hji)]; ring)
        (fun h => absurd (mem_univ i) h), real_inner_self_eq_norm_sq]
      ring
    simp_rw [this, Real.norm_eq_abs, sq_abs]
    refine sum_le_sum fun i _ => ?_
    have : ‖q i‖ ^ 2 ≤ 1 := by
      have := hq i; nlinarith [norm_nonneg (q i)]
    nlinarith [sq_nonneg (y i)]
  exact (pow_le_pow_iff_left₀ (norm_nonneg _) (norm_nonneg _) two_ne_zero).mp hsq

end Projection

section Blocks

variable {κ : Type*} [Fintype κ]

/-- Block-wise radial map on the read space. -/
noncomputable def blockRadial (r : ℝ → ℝ) (v : R κ) : R κ :=
  WithLp.toLp 2 fun k => radial r (v k)

theorem blockRadial_lipschitz {r : ℝ → ℝ} (hr : GoodRadius r) :
    LipschitzWith 1 (blockRadial (κ := κ) r) := by
  refine LipschitzWith.of_dist_le_mul fun v w => ?_
  rw [NNReal.coe_one, one_mul]
  have hsq : dist (blockRadial r v) (blockRadial r w) ^ 2 ≤ dist v w ^ 2 := by
    rw [PiLp.dist_sq_eq_of_L2, PiLp.dist_sq_eq_of_L2]
    refine sum_le_sum fun k _ => ?_
    simp only [blockRadial, PiLp.toLp_apply, dist_eq_norm]
    exact pow_le_pow_left₀ (norm_nonneg _) (radial_lipschitz hr _ _) 2
  exact (pow_le_pow_iff_left₀ dist_nonneg dist_nonneg two_ne_zero).mp hsq

/-- Norm of a block-wise radial image: `‖φ(v)‖² = ∑_k r(‖v_k‖)²`. -/
theorem norm_blockRadial_sq {r : ℝ → ℝ} (hr : GoodRadius r) (v : R κ) :
    ‖blockRadial r v‖ ^ 2 = ∑ k, r ‖v k‖ ^ 2 := by
  rw [PiLp.norm_sq_eq_of_L2]
  simp only [blockRadial, PiLp.toLp_apply, norm_radial hr]

end Blocks

section Tail

variable {ι κ σ : Type*} [Fintype ι] [DecidableEq ι] [Fintype κ] [Fintype σ] [Nonempty σ]
  {F : Type*} [NormedAddCommGroup F] [InnerProductSpace ℝ F] [CompleteSpace F]

/-- The reading `A_s` with values in the read space (rows indexed by block and component). -/
noncomputable def readBlocks (S : κ × Fin 3 → Finset ι) (k : κ × Fin 3 → ι → ℝ)
    (z : EuclideanSpace ℝ ι) : R κ :=
  WithLp.toLp 2 fun b => WithLp.toLp 2 fun c => ∑ i ∈ S (b, c), k (b, c) i * z i

theorem readBlocks_lipschitz (S : κ × Fin 3 → Finset ι)
    (hS : Set.PairwiseDisjoint (↑(univ : Finset (κ × Fin 3))) S) (k : κ × Fin 3 → ι → ℝ)
    (hk : ∀ r, ∑ i ∈ S r, k r i ^ 2 ≤ 1) : LipschitzWith 1 (readBlocks S k) := by
  refine LipschitzWith.of_dist_le_mul fun x y => ?_
  rw [NNReal.coe_one, one_mul]
  have hsq : dist (readBlocks S k x) (readBlocks S k y) ^ 2 ≤ dist x y ^ 2 := by
    rw [PiLp.dist_sq_eq_of_L2, PiLp.dist_sq_eq_of_L2 x y]
    simp_rw [PiLp.dist_sq_eq_of_L2]
    simp only [readBlocks, PiLp.toLp_apply, Real.dist_eq, sq_abs]
    rw [← Fintype.sum_prod_type' (fun b c =>
      (∑ i ∈ S (b, c), k (b, c) i * x i - ∑ i ∈ S (b, c), k (b, c) i * y i) ^ 2)]
    have := reading_contraction S hS k hk (fun i => x i - y i)
    simpa [← Finset.sum_sub_distrib, mul_sub] using this
  exact (pow_le_pow_iff_left₀ dist_nonneg dist_nonneg two_ne_zero).mp hsq

/-- The function `F(z) = max_s ‖Q† φ(A_s z)‖` of the proof of AB.C.1. -/
noncomputable def tailFunction (S : σ → κ × Fin 3 → Finset ι) (k : σ → κ × Fin 3 → ι → ℝ)
    (r : ℝ → ℝ) (Q : F →L[ℝ] R κ) (z : EuclideanSpace ℝ ι) : ℝ :=
  univ.sup' univ_nonempty fun s =>
    ‖ContinuousLinearMap.adjoint Q (blockRadial r (readBlocks (S s) (k s) z))‖

/-- **AB.C.1, deterministic part, complete.** `F` is 1-Lipschitz in the white noise. -/
theorem tailFunction_lipschitz (S : σ → κ × Fin 3 → Finset ι)
    (hS : ∀ s, Set.PairwiseDisjoint (↑(univ : Finset (κ × Fin 3))) (S s))
    (k : σ → κ × Fin 3 → ι → ℝ) (hk : ∀ s r, ∑ i ∈ S s r, k s r i ^ 2 ≤ 1)
    {r : ℝ → ℝ} (hr : GoodRadius r) (Q : F →L[ℝ] R κ) (hQ : ∀ y, ‖Q y‖ ≤ ‖y‖) :
    LipschitzWith 1 (tailFunction S k r Q) := by
  refine lipschitz_sup univ univ_nonempty _ fun s _ => ?_
  have hQt : LipschitzWith 1 (ContinuousLinearMap.adjoint Q) := by
    refine LipschitzWith.of_dist_le_mul fun a b => ?_
    rw [NNReal.coe_one, one_mul, dist_eq_norm, dist_eq_norm, ← map_sub]
    exact norm_adjoint_le Q hQ _
  have h := lipschitzWith_one_norm.comp
    (hQt.comp ((blockRadial_lipschitz hr).comp (readBlocks_lipschitz (S s) (hS s) (k s) (hk s))))
  simp only [mul_one] at h
  exact h

/-- The score at any (data-dependent) reading `s⋆(z)` is at most `F(z)`. -/
theorem read_le_tailFunction (S : σ → κ × Fin 3 → Finset ι) (k : σ → κ × Fin 3 → ι → ℝ)
    (r : ℝ → ℝ) (Q : F →L[ℝ] R κ) (sel : EuclideanSpace ℝ ι → σ) (z : EuclideanSpace ℝ ι) :
    ‖ContinuousLinearMap.adjoint Q (blockRadial r (readBlocks (S (sel z)) (k (sel z)) z))‖ ≤
      tailFunction S k r Q z :=
  le_sup' (fun s => ‖ContinuousLinearMap.adjoint Q
    (blockRadial r (readBlocks (S s) (k s) z))‖) (mem_univ (sel z))

end Tail

section Influence

variable {κ : Type*} [Fintype κ] [DecidableEq κ]
  {F : Type*} [NormedAddCommGroup F] [InnerProductSpace ℝ F] [CompleteSpace F]

/-- **Several blocks (general form behind AA.9.2 and AA.9.5(c)).** If `v` and `v'` differ
only in the blocks of `D`, then
`|σ(v) − σ(v')| ≤ √(∑_{k∈D} (r‖v_k‖ + r‖v'_k‖)²)`. -/
theorem score_sub_le_blocks {r : ℝ → ℝ} (hr : GoodRadius r) (Q : F →L[ℝ] R κ)
    (hQ : ∀ y, ‖Q y‖ ≤ ‖y‖) (v v' : R κ) (D : Finset κ) (hD : ∀ k ∉ D, v k = v' k) :
    |‖ContinuousLinearMap.adjoint Q (blockRadial r v)‖ -
        ‖ContinuousLinearMap.adjoint Q (blockRadial r v')‖| ≤
      √(∑ k ∈ D, (r ‖v k‖ + r ‖v' k‖) ^ 2) := by
  refine (norm_adjoint_sub_le Q hQ _ _).trans ?_
  rw [← Real.sqrt_sq (norm_nonneg _)]
  refine Real.sqrt_le_sqrt ?_
  rw [PiLp.norm_sq_eq_of_L2, ← Finset.sum_subset (subset_univ D) ?_]
  · refine sum_le_sum fun k _ => ?_
    simp only [blockRadial, PiLp.sub_apply, PiLp.toLp_apply]
    have h1 := norm_sub_le (radial r (v k)) (radial r (v' k))
    rw [norm_radial hr, norm_radial hr] at h1
    exact pow_le_pow_left₀ (norm_nonneg _) h1 2
  · intro k _ hk
    simp [blockRadial, hD k hk]

/-- **AA.9.2 (influence of one block).** -/
theorem influence_le {r : ℝ → ℝ} (hr : GoodRadius r) (Q : F →L[ℝ] R κ)
    (hQ : ∀ y, ‖Q y‖ ≤ ‖y‖) (v v' : R κ) (b : κ) (hD : ∀ k, k ≠ b → v k = v' k) :
    |‖ContinuousLinearMap.adjoint Q (blockRadial r v)‖ -
        ‖ContinuousLinearMap.adjoint Q (blockRadial r v')‖| ≤ r ‖v b‖ + r ‖v' b‖ := by
  have h := score_sub_le_blocks hr Q hQ v v' {b} (fun k hk => hD k (by simpa using hk))
  rw [sum_singleton, Real.sqrt_sq (add_nonneg (hr.nonneg _ (norm_nonneg _))
    (hr.nonneg _ (norm_nonneg _)))] at h
  exact h

/-- **AA.9.5(c) for the cap.** With at most `m` contaminated blocks, the capped score moves
by at most `2c√m`. -/
theorem cap_contamination {c : ℝ} (hc : 0 ≤ c) (Q : F →L[ℝ] R κ) (hQ : ∀ y, ‖Q y‖ ≤ ‖y‖)
    (v v' : R κ) (D : Finset κ) (hD : ∀ k ∉ D, v k = v' k) :
    |‖ContinuousLinearMap.adjoint Q (blockRadial (fun t => min t c) v)‖ -
        ‖ContinuousLinearMap.adjoint Q (blockRadial (fun t => min t c) v')‖| ≤
      2 * c * √(D.card) := by
  refine (score_sub_le_blocks (goodRadius_cap hc) Q hQ v v' D hD).trans ?_
  have : ∑ k ∈ D, (min ‖v k‖ c + min ‖v' k‖ c) ^ 2 ≤ (2 * c) ^ 2 * D.card := by
    rw [mul_comm, ← nsmul_eq_mul, ← sum_const]
    refine sum_le_sum fun k _ => pow_le_pow_left₀ (add_nonneg (le_min (norm_nonneg _) hc)
      (le_min (norm_nonneg _) hc)) (by linarith [min_le_right ‖v k‖ c, min_le_right ‖v' k‖ c]) 2
  calc √(∑ k ∈ D, (min ‖v k‖ c + min ‖v' k‖ c) ^ 2) ≤ √((2 * c) ^ 2 * D.card) :=
        Real.sqrt_le_sqrt this
    _ = 2 * c * √(D.card) := by
        rw [Real.sqrt_mul (sq_nonneg _), Real.sqrt_sq (by linarith)]

/-- **AA.9.3, first claim (saturation bound).** `σ_{r_c}² ≤ n_b c²`, i.e. `s ≤ n_b c²/2`. -/
theorem cap_score_le {c : ℝ} (hc : 0 ≤ c) (Q : F →L[ℝ] R κ) (hQ : ∀ y, ‖Q y‖ ≤ ‖y‖)
    (v : R κ) :
    ‖ContinuousLinearMap.adjoint Q (blockRadial (fun t => min t c) v)‖ ^ 2 ≤
      Fintype.card κ * c ^ 2 := by
  calc ‖ContinuousLinearMap.adjoint Q (blockRadial (fun t => min t c) v)‖ ^ 2
      ≤ ‖blockRadial (fun t => min t c) v‖ ^ 2 :=
        pow_le_pow_left₀ (norm_nonneg _) (norm_adjoint_le Q hQ _) 2
    _ = ∑ k, min ‖v k‖ c ^ 2 := norm_blockRadial_sq (goodRadius_cap hc) v
    _ ≤ ∑ _k : κ, c ^ 2 := sum_le_sum fun k _ =>
        pow_le_pow_left₀ (le_min (norm_nonneg _) hc) (min_le_right _ _) 2
    _ = Fintype.card κ * c ^ 2 := by rw [sum_const, card_univ, nsmul_eq_mul]

end Influence

end TTMT
