/-
Orthonormal combinations of independent standard Gaussians are independent standard
Gaussians (`rotation_iid`), and a quadratic form in Gaussian noise is a weighted sum of
independent chi-square variables in the eigenbasis (`quadForm_eigen`). Together with
`M6NullTail.lean` this gives the full tail of Proposition E.2 for `ε = L X`, `C = L Lᵀ`
(`gain_null_tail`).
-/
import TTMT.M6NullTail
import TTMT.M6NullGain
import Mathlib.Probability.Distributions.Gaussian.HasGaussianLaw.Independence
import Mathlib.Probability.Moments.Covariance
import Mathlib.Analysis.Matrix.Spectrum

open MeasureTheory ProbabilityTheory Real Finset Matrix

namespace TTMT

variable {Ω ι : Type*} [MeasurableSpace Ω] {P : Measure Ω} [IsProbabilityMeasure P]
  [Fintype ι] [DecidableEq ι]

/-- **Rotation invariance.** If `X_i` are independent `N(0,1)` and the rows `u_j` are
orthonormal, then `Y_j = ∑_i u_j i X_i` are independent `N(0,1)`. -/
theorem rotation_iid (X : ι → Ω → ℝ) (hmeas : ∀ i, Measurable (X i))
    (hind : iIndepFun X P) (hlaw : ∀ i, P.map (X i) = gaussianReal 0 1) (u : ι → ι → ℝ)
    (horth : ∀ j k, ∑ i, u j i * u k i = if j = k then 1 else 0) :
    (∀ j, Measurable fun ω => ∑ i, u j i * X i ω) ∧
    iIndepFun (fun j ω => ∑ i, u j i * X i ω) P ∧
    ∀ j, P.map (fun ω => ∑ i, u j i * X i ω) = gaussianReal 0 1 := by
  set Y : ι → Ω → ℝ := fun j ω => ∑ i, u j i * X i ω
  have hYm : ∀ j, Measurable (Y j) := fun j =>
    Finset.measurable_sum _ fun i _ => measurable_const.mul (hmeas i)
  have hXg : ∀ i, HasGaussianLaw (X i) P := fun i =>
    (HasLaw.hasGaussianLaw (μ := gaussianReal 0 1) ⟨(hmeas i).aemeasurable, hlaw i⟩)
  have hvec : HasGaussianLaw (fun ω => (X · ω)) P := hind.hasGaussianLaw hXg
  let L : (ι → ℝ) →L[ℝ] (ι → ℝ) := LinearMap.toContinuousLinearMap (Matrix.mulVecLin u)
  have hYvec : HasGaussianLaw (fun ω => (Y · ω)) P := by
    have := hvec.map_fun L
    convert this using 2 with ω
    funext j
    show _ = Matrix.mulVec u _ j
    simp [Y, Matrix.mulVec, dotProduct]
  have hmem : ∀ i, MemLp (X i) 2 P := fun i => (hXg i).memLp_two
  have hvar : ∀ i, Var[X i; P] = 1 := by
    intro i
    rw [← variance_id_map (hmeas i).aemeasurable, hlaw i]
    simp
  have hmean : ∀ i, P[X i] = 0 := by
    intro i
    calc P[X i] = ∫ x, id x ∂(P.map (X i)) :=
          (integral_map (hmeas i).aemeasurable aestronglyMeasurable_id).symm
      _ = 0 := by rw [hlaw i]; simp
  have hcovX : ∀ i l, cov[X i, X l; P] = if i = l then 1 else 0 := by
    intro i l
    split_ifs with h
    · subst h; rw [covariance_self (hmeas i).aemeasurable, hvar]
    · exact (hind.indepFun h).covariance_eq_zero (hmem i) (hmem l)
  have hcov : ∀ j k, cov[Y j, Y k; P] = if j = k then 1 else 0 := by
    intro j k
    simp only [Y]
    rw [covariance_fun_sum_fun_sum (fun i => (hmem i).const_mul _)
      (fun i => (hmem i).const_mul _)]
    simp_rw [covariance_const_mul_left, covariance_const_mul_right, hcovX]
    rw [← horth j k]
    refine sum_congr rfl fun i _ => ?_
    rw [Finset.sum_eq_single i (fun l _ hl => by simp [Ne.symm hl])
      (fun h => absurd (mem_univ i) h)]
    simp
  refine ⟨hYm, ?_, ?_⟩
  · exact hYvec.iIndepFun_of_covariance_eq_zero fun j k hjk => by simp [hcov, hjk]
  · intro j
    have h := (hYvec.eval j).map_eq_gaussianReal
    have hm : P[Y j] = 0 := by
      simp only [Y]
      rw [integral_finsetSum _ fun i _ => ((hmem i).integrable one_le_two).const_mul _]
      simp [integral_const_mul, hmean]
    have hv : Var[Y j; P] = 1 := by
      rw [← covariance_self (hYm j).aemeasurable, hcov]; simp
    rw [h, hm, hv]
    simp

/-- The eigenvector rows of a symmetric matrix are orthonormal. -/
theorem eigen_rows_orthonormal (B : Matrix ι ι ℝ) (hB : B.IsHermitian) (j k : ι) :
    ∑ i, hB.eigenvectorBasis j i * hB.eigenvectorBasis k i = if j = k then 1 else 0 := by
  have h := orthonormal_iff_ite.mp hB.eigenvectorBasis.orthonormal j k
  rw [EuclideanSpace.inner_eq_star_dotProduct, star_trivial] at h
  rw [← h, dotProduct]
  exact sum_congr rfl fun i _ => mul_comm _ _

/-- **Spectral form of a quadratic form.** `xᵀ B x = ∑_j λ_j ⟨u_j, x⟩²` with `u_j` the
orthonormal eigenvectors of the symmetric matrix `B`. -/
theorem quadForm_eigen (B : Matrix ι ι ℝ) (hB : B.IsHermitian) (x : ι → ℝ) :
    x ⬝ᵥ (B *ᵥ x) =
      ∑ j, hB.eigenvalues j * (∑ i, hB.eigenvectorBasis j i * x i) ^ 2 := by
  set b := hB.eigenvectorBasis
  set c : ι → ℝ := fun j => ∑ i, b j i * x i
  have hc : ∀ j, inner ℝ (b j) (WithLp.toLp 2 x) = c j := by
    intro j
    rw [EuclideanSpace.inner_toLp_toLp, star_trivial, dotProduct]
    exact sum_congr rfl fun i _ => mul_comm _ _
  have hx : x = ∑ j, c j • ⇑(b j) := by
    have h := congrArg WithLp.ofLp (b.sum_repr' (WithLp.toLp 2 x))
    simp only [hc, WithLp.ofLp_sum, WithLp.ofLp_smul] at h
    exact h.symm
  have hBx : B *ᵥ x = ∑ j, (c j * hB.eigenvalues j) • ⇑(b j) := by
    conv_lhs => rw [hx]
    rw [Matrix.mulVec_sum]
    refine sum_congr rfl fun j _ => ?_
    rw [Matrix.mulVec_smul, hB.mulVec_eigenvectorBasis, smul_smul]
  rw [hBx, dotProduct_sum]
  refine sum_congr rfl fun j _ => ?_
  rw [dotProduct_smul, smul_eq_mul]
  have : x ⬝ᵥ ⇑(b j) = c j := by
    simp only [c, dotProduct]; exact sum_congr rfl fun i _ => mul_comm _ _
  rw [this]; ring

/-- **Proposition E.2, tail (full).** Noise `ε = L X` with `X` i.i.d. `N(0,1)`, so `ε` has
covariance `C = L Lᵀ`; `A` is symmetric (for E.2, `A = P₂ − P₁`) and `λ_j` are the eigenvalues
of `B = Lᵀ A L`. Then the null gain `εᵀ A ε` exceeds `tr(A C) + 2√(v⁺τ) + 2aτ + 2√(v⁻τ)`
with probability at most `2e^{-τ}`, where `v^± = ∑ (λ_j)_±²` and `a ≥ max_j (λ_j)_+`. -/
theorem gain_null_tail (X : ι → Ω → ℝ) (hmeas : ∀ i, Measurable (X i))
    (hind : iIndepFun X P) (hlaw : ∀ i, P.map (X i) = gaussianReal 0 1)
    (L A : Matrix ι ι ℝ) (hB : (Lᵀ * A * L).IsHermitian)
    {a : ℝ} (ha : ∀ j, max (hB.eigenvalues j) 0 ≤ a)
    (hvp : 0 < ∑ j, max (hB.eigenvalues j) 0 ^ 2)
    (hvn : 0 < ∑ j, max (-hB.eigenvalues j) 0 ^ 2) {τ : ℝ} (hτ : 0 < τ) :
    P.real {ω | (A * (L * Lᵀ)).trace + 2 * √((∑ j, max (hB.eigenvalues j) 0 ^ 2) * τ) +
        2 * a * τ + 2 * √((∑ j, max (-hB.eigenvalues j) 0 ^ 2) * τ) ≤
        (L *ᵥ fun i => X i ω) ⬝ᵥ (A *ᵥ (L *ᵥ fun i => X i ω))} ≤ 2 * exp (-τ) := by
  set lam := hB.eigenvalues
  set u : ι → ι → ℝ := fun j i => hB.eigenvectorBasis j i
  obtain ⟨hYm, hYind, hYlaw⟩ :=
    rotation_iid X hmeas hind hlaw u (eigen_rows_orthonormal _ hB)
  have hT : (A * (L * Lᵀ)).trace = ∑ j, lam j := by
    rw [← Matrix.mul_assoc, Matrix.trace_mul_comm, ← Matrix.mul_assoc,
      hB.trace_eq_sum_eigenvalues]
    simp [lam]
  have hQ : ∀ x : ι → ℝ, (L *ᵥ x) ⬝ᵥ (A *ᵥ (L *ᵥ x)) =
      ∑ j, lam j * (∑ i, u j i * x i) ^ 2 := by
    intro x
    rw [← quadForm_eigen _ hB x, ← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec,
      Matrix.dotProduct_mulVec x Lᵀ, Matrix.vecMul_transpose]
  have h := gain_tail_eigen (fun j ω => ∑ i, u j i * X i ω) hYm hYind hYlaw lam ha hvp hvn hτ
  convert h using 2
  ext ω
  simp only [Set.mem_ofPred_eq]
  rw [hT, hQ]

omit [DecidableEq ι] in
/-- `Lᵀ A L` is symmetric whenever `A` is; this supplies the hypothesis `hB` above. -/
theorem sandwich_isHermitian (L A : Matrix ι ι ℝ) (hA : Aᵀ = A) : (Lᵀ * A * L).IsHermitian := by
  unfold Matrix.IsHermitian
  rw [conjTranspose_eq_transpose_of_trivial, transpose_mul, transpose_mul, transpose_transpose,
    hA, Matrix.mul_assoc]

/-- **Proposition E.2 in its original form.** For orthogonal projectors `P₁, P₂` that both
contain the signal `μ`, and noise `ε = L X`, the RSS gain `RSS₁ − RSS₂` obeys the tail bound
of `gain_null_tail` with `A = P₂ − P₁`; the threshold's centre is `tr((P₂ − P₁) C)`. -/
theorem gain_null_tail_rss (X : ι → Ω → ℝ) (hmeas : ∀ i, Measurable (X i))
    (hind : iIndepFun X P) (hlaw : ∀ i, P.map (X i) = gaussianReal 0 1)
    (P₁ P₂ L : Matrix ι ι ℝ) (hs₁ : P₁ᵀ = P₁) (hi₁ : P₁ * P₁ = P₁)
    (hs₂ : P₂ᵀ = P₂) (hi₂ : P₂ * P₂ = P₂) (μ : ι → ℝ)
    (hμ₁ : P₁ *ᵥ μ = μ) (hμ₂ : P₂ *ᵥ μ = μ) (hB : (Lᵀ * (P₂ - P₁) * L).IsHermitian)
    {a : ℝ} (ha : ∀ j, max (hB.eigenvalues j) 0 ≤ a)
    (hvp : 0 < ∑ j, max (hB.eigenvalues j) 0 ^ 2)
    (hvn : 0 < ∑ j, max (-hB.eigenvalues j) 0 ^ 2) {τ : ℝ} (hτ : 0 < τ) :
    P.real {ω | ((P₂ - P₁) * (L * Lᵀ)).trace +
        2 * √((∑ j, max (hB.eigenvalues j) 0 ^ 2) * τ) + 2 * a * τ +
        2 * √((∑ j, max (-hB.eigenvalues j) 0 ^ 2) * τ) ≤
        rss P₁ (μ + L *ᵥ fun i => X i ω) - rss P₂ (μ + L *ᵥ fun i => X i ω)} ≤
      2 * exp (-τ) := by
  have h := gain_null_tail X hmeas hind hlaw L (P₂ - P₁) hB ha hvp hvn hτ
  convert h using 2
  ext ω
  simp only [Set.mem_ofPred_eq]
  rw [gain_null_eq P₁ P₂ hs₁ hi₁ hs₂ hi₂ μ _ hμ₁ hμ₂]

end TTMT
