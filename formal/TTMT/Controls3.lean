/-
Controls for the probabilistic results (rule 40). Negative controls drop ONE named hypothesis
and prove the conclusion fails on an explicit configuration (Gaussian masses are bounded below
by `gaussian_Icc_ge`). Positive controls instantiate the theorem on a concrete configuration
(one or two independent `N(0,1)`) to show its hypotheses can all hold at once.

Hypotheses with a negative control:
- `neg_log_one_sub_le`: `0 ≤ x` (`x = −3`), `x < 1` (`x = 1.1`). `log_one_add_ge`: `0 ≤ u`.
- `laurent_massart_upper`: independence (two copies of one `N(0,1)`), `λ ≥ 0`, `λ ≤ a`,
  `v > 0`, the `N(0,1)` law (`X ≡ 3`).
- `laurent_massart_lower`: independence (64 copies), `λ ≥ 0`, `v > 0`, the law (`X ≡ 0`).
- `rotation_iid`: orthonormal rows (`u = 2`), independence, the law (`X ≡ 0`).
- `gain_tail_eigen`: independence (`lam5` on `X5`), `max(λ_j,0) ≤ a`, the law, and
  `v⁺ > 0`, `v⁻ > 0` together.

Hypotheses that are not needed for the conclusion to be true (kept for the proof), shown here:
- `τ > 0` in every tail bound: `tail_trivial_of_nonpos` (for `τ ≤ 0`, `e^{−τ} ≥ 1`).
- Measurability: implied by the law hypothesis, `law_forces_aemeasurable`.
- `v⁺ > 0` alone or `v⁻ > 0` alone in `gain_tail_eigen` (the other tail carries the bound).

Exempt, with the reason: `quadForm_eigen`, `eigen_rows_orthonormal`, `gain_null_tail` (`B`
symmetric is a typing dependency of the eigenvalues; the other hypotheses are passed to
`gain_tail_eigen` and `rotation_iid`, which have controls); `gain_null_tail_rss` (adds the
hypotheses of `gain_null_eq`, controlled in `Controls.lean`); `sandwich_isHermitian` and the
moment-generating-function lemmas (intermediate steps of the controlled theorems).
-/
import TTMT.GaussianRotation
import TTMT.LogBounds

open MeasureTheory ProbabilityTheory Real Finset Matrix Set

namespace TTMT.Controls3

/-! ### Two tools for the negative controls -/

/-- A Gaussian interval has mass at least `length × density at the far end` (`0 ≤ a ≤ b`). -/
theorem gaussian_Icc_ge {a b : ℝ} (ha : 0 ≤ a) (hab : a ≤ b) :
    (b - a) * gaussianPDFReal 0 1 b ≤ (gaussianReal 0 1).real (Icc a b) := by
  rw [measureReal_def, gaussianReal_apply_eq_integral 0 one_ne_zero,
    ENNReal.toReal_ofReal (setIntegral_nonneg measurableSet_Icc
      fun x _ => gaussianPDFReal_nonneg _ _ _)]
  calc (b - a) * gaussianPDFReal 0 1 b = ∫ _ in Icc a b, gaussianPDFReal 0 1 b := by
        rw [setIntegral_const, volume_real_Icc_of_le hab, smul_eq_mul]
    _ ≤ ∫ x in Icc a b, gaussianPDFReal 0 1 x := by
        refine setIntegral_mono_on (integrableOn_const (by simp))
          (integrable_gaussianPDFReal 0 1).integrableOn measurableSet_Icc fun x hx => ?_
        simp only [gaussianPDFReal, NNReal.coe_one, mul_one, sub_zero]
        gcongr
        all_goals linarith [hx.1, hx.2]

/-- Numeric comparison used by every tail control: if `2.9 < (b − a) e^{τ − b²/2}` then
`e^{−τ} < (b − a) φ(b)`. -/
theorem exp_lt_mass {a b τ : ℝ} (key : 2.9 < (b - a) * exp (τ - b ^ 2 / 2)) :
    exp (-τ) < (b - a) * gaussianPDFReal 0 1 b := by
  have hs : √(2 * π) < 2.9 := by
    rw [Real.sqrt_lt' (by norm_num)]; nlinarith [pi_le_four]
  have hs0 : 0 < √(2 * π) := Real.sqrt_pos.mpr (by positivity)
  simp only [gaussianPDFReal, NNReal.coe_one, mul_one, sub_zero]
  rw [← mul_assoc, mul_comm (b - a), mul_assoc, inv_mul_eq_div, lt_div_iff₀ hs0]
  have he : exp (τ - b ^ 2 / 2) = exp τ * exp (-b ^ 2 / 2) := by
    rw [← exp_add]; ring_nf
  rw [he] at key
  have hτ : 0 < exp (-τ) := exp_pos _
  have h1 : exp (-τ) * exp τ = 1 := by rw [← exp_add]; simp
  calc exp (-τ) * √(2 * π) < exp (-τ) * 2.9 := by gcongr
    _ < exp (-τ) * ((b - a) * (exp τ * exp (-b ^ 2 / 2))) := by gcongr
    _ = (b - a) * exp (-b ^ 2 / 2) := by
        rw [show exp (-τ) * ((b - a) * (exp τ * exp (-b ^ 2 / 2))) =
          (exp (-τ) * exp τ) * ((b - a) * exp (-b ^ 2 / 2)) by ring, h1, one_mul]

/-- Under `N(0,1)` on `ℝ`, the identity has law `N(0,1)` (used by every instance below). -/
theorem map_id_gaussian : (gaussianReal 0 1).map (fun ω : ℝ => ω) = gaussianReal 0 1 :=
  Measure.map_id

/-! ### Log bounds -/

/-- Positive control: `neg_log_one_sub_le` at `x = 1/2`. -/
example : -log (1 - 1 / 2 : ℝ) - 1 / 2 ≤ (1 / 2 : ℝ) ^ 2 / (2 * (1 - 1 / 2)) :=
  neg_log_one_sub_le (by norm_num) (by norm_num)

/-- Negative control, `0 ≤ x` dropped: at `x = −3`, `−log 4 + 3 > 9/8` (as `log 4 < 1.875`). -/
theorem neg_log_needs_nonneg :
    ¬ (-log (1 - (-3)) - (-3) ≤ (-3 : ℝ) ^ 2 / (2 * (1 - (-3)))) := by
  have h : log 4 < 1.875 := by
    rw [Real.log_lt_iff_lt_exp (by norm_num)]
    have := quadratic_le_exp_of_nonneg (x := 1.875) (by norm_num)
    linarith
  norm_num; linarith

/-- Negative control, `x < 1` dropped: at `x = 1.1` the left side is `> −1.1` (`log 0.1 < 0`)
and the right side is `−6.05`. -/
theorem neg_log_needs_lt_one :
    ¬ (-log (1 - 1.1) - 1.1 ≤ (1.1 : ℝ) ^ 2 / (2 * (1 - 1.1))) := by
  have h : log (1 - 1.1 : ℝ) < 0 := by
    rw [show (1 - 1.1 : ℝ) = -0.1 by norm_num, Real.log_neg_eq_log]
    exact Real.log_neg (by norm_num) (by norm_num)
  have : (1.1 : ℝ) ^ 2 / (2 * (1 - 1.1)) = -6.05 := by norm_num
  rw [this]; linarith

/-- Positive control: `log_one_add_ge` at `u = 1`. -/
example : (1 : ℝ) - 1 ^ 2 / 2 ≤ log (1 + 1) := log_one_add_ge zero_le_one

/-- Negative control, `0 ≤ u` dropped: at `u = −0.99`, `log 0.01 = 2 log 0.1 ≤ −1.8` while
`u − u²/2 = −1.48005`. -/
theorem log_add_needs_nonneg : ¬ ((-0.99 : ℝ) - (-0.99) ^ 2 / 2 ≤ log (1 + (-0.99))) := by
  have h1 : log (0.1 : ℝ) ≤ 0.1 - 1 := Real.log_le_sub_one_of_pos (by norm_num)
  have h2 : log (1 + (-0.99) : ℝ) = 2 * log 0.1 := by
    rw [show (1 + (-0.99) : ℝ) = 0.1 ^ 2 by norm_num, Real.log_pow]; norm_num
  rw [h2]; norm_num; linarith

/-! ### Laurent–Massart, upper tail -/

/-- Positive control: every hypothesis of `laurent_massart_upper` holds for one `N(0,1)`
variable (`Ω = ℝ`, `X = id`, `λ = 1`, `a = 1`, `τ = 1`). -/
example : (gaussianReal 0 1).real {ω : ℝ | 2 * √((∑ _j : Unit, (1 : ℝ) ^ 2) * 1) + 2 * 1 * 1 ≤
    ∑ _j : Unit, 1 * (ω ^ 2 - 1)} ≤ exp (-1) :=
  laurent_massart_upper (fun _ ω => ω) (fun _ => measurable_id) iIndepFun.of_subsingleton
    (fun _ => map_id_gaussian) (fun _ => 1) (fun _ => zero_le_one) (fun _ => le_rfl)
    (by simp) one_pos

/-- Negative control, independence dropped: two copies of one `N(0,1)` (`λ = (1, 1)`, `a = 1`,
`τ = 32`). The event is `Z² ≥ 41`, which contains `[6.5, 7.5]`, of mass `> e^{−32}`. -/
theorem lm_upper_needs_indep :
    ¬ ((gaussianReal 0 1).real {ω : ℝ | 2 * √((∑ _j : Fin 2, (1 : ℝ) ^ 2) * 32) + 2 * 1 * 32 ≤
      ∑ _j : Fin 2, 1 * (ω ^ 2 - 1)} ≤ exp (-32)) := by
  have hsq : √((∑ _j : Fin 2, (1 : ℝ) ^ 2) * 32) = 8 := by
    rw [show (∑ _j : Fin 2, (1 : ℝ) ^ 2) * 32 = 8 ^ 2 by norm_num]
    exact Real.sqrt_sq (by norm_num)
  rw [not_le]
  refine (exp_lt_mass (a := 6.5) (b := 7.5) ?_).trans_le
    ((gaussian_Icc_ge (by norm_num) (by norm_num)).trans (measureReal_mono ?_))
  · have := add_one_le_exp (32 - 7.5 ^ 2 / 2); norm_num at this ⊢; linarith
  · intro x hx
    rw [Set.mem_ofPred_eq, hsq]
    simp only [Fin.sum_univ_two]
    nlinarith [hx.1]

/-- Negative control, `λ ≥ 0` dropped: `λ = −1`, `a = −1`, `τ = 4`. The event is `Z² ≤ 5`,
which contains `[0, 2]`, of mass `> e^{−4}`. -/
theorem lm_upper_needs_nonneg :
    ¬ ((gaussianReal 0 1).real {ω : ℝ | 2 * √((∑ _j : Unit, (-1 : ℝ) ^ 2) * 4) + 2 * (-1) * 4 ≤
      ∑ _j : Unit, (-1) * (ω ^ 2 - 1)} ≤ exp (-4)) := by
  have hsq : √((∑ _j : Unit, (-1 : ℝ) ^ 2) * 4) = 2 := by
    rw [show (∑ _j : Unit, (-1 : ℝ) ^ 2) * 4 = 2 ^ 2 by simp; norm_num]
    exact Real.sqrt_sq (by norm_num)
  rw [not_le]
  refine (exp_lt_mass (a := 0) (b := 2) ?_).trans_le
    ((gaussian_Icc_ge le_rfl (by norm_num)).trans (measureReal_mono ?_))
  · have := add_one_le_exp (4 - 2 ^ 2 / 2); norm_num at this ⊢; linarith
  · intro x hx
    rw [Set.mem_ofPred_eq, hsq]
    simp only [Finset.univ_unique, Finset.sum_singleton]
    nlinarith [hx.1, hx.2]

/-- Negative control, `λ ≤ a` dropped: `λ = 1`, `a = −100`, `τ = 1`. The threshold is `−198`
and `Z² − 1 ≥ −1`, so the event is everything, of mass `1 > e^{−1}`. -/
theorem lm_upper_needs_le_a :
    ¬ ((gaussianReal 0 1).real {ω : ℝ | 2 * √((∑ _j : Unit, (1 : ℝ) ^ 2) * 1) + 2 * (-100) * 1 ≤
      ∑ _j : Unit, 1 * (ω ^ 2 - 1)} ≤ exp (-1)) := by
  have : {ω : ℝ | 2 * √((∑ _j : Unit, (1 : ℝ) ^ 2) * 1) + 2 * (-100) * 1 ≤
      ∑ _j : Unit, 1 * (ω ^ 2 - 1)} = univ := by
    refine eq_univ_of_forall fun x => ?_
    simp only [Set.mem_ofPred_eq, Finset.univ_unique, Finset.sum_singleton]
    norm_num; nlinarith [sq_nonneg x]
  rw [this, probReal_univ, not_le]
  exact (exp_lt_exp.mpr (by norm_num)).trans_eq exp_zero

/-- Negative control, `v > 0` dropped: `λ = 0`, `a = 0`, `τ = 1`; the event `0 ≤ 0` is
everything. -/
theorem lm_upper_needs_v_pos :
    ¬ ((gaussianReal 0 1).real {ω : ℝ | 2 * √((∑ _j : Unit, (0 : ℝ) ^ 2) * 1) + 2 * 0 * 1 ≤
      ∑ _j : Unit, 0 * (ω ^ 2 - 1)} ≤ exp (-1)) := by
  have : {ω : ℝ | 2 * √((∑ _j : Unit, (0 : ℝ) ^ 2) * 1) + 2 * 0 * 1 ≤
      ∑ _j : Unit, 0 * (ω ^ 2 - 1)} = univ := eq_univ_of_forall fun x => by simp
  rw [this, probReal_univ, not_le]
  exact (exp_lt_exp.mpr (by norm_num)).trans_eq exp_zero

/-- Negative control, the `N(0,1)` law dropped: `X ≡ 3` (measurable, trivially independent),
`λ = 1`, `a = 1`, `τ = 1`. The event `4 ≤ 8` is everything. -/
theorem lm_upper_needs_law :
    ¬ ((gaussianReal 0 1).real {_ω : ℝ | 2 * √((∑ _j : Unit, (1 : ℝ) ^ 2) * 1) + 2 * 1 * 1 ≤
      ∑ _j : Unit, 1 * ((3 : ℝ) ^ 2 - 1)} ≤ exp (-1)) := by
  have : {_ω : ℝ | 2 * √((∑ _j : Unit, (1 : ℝ) ^ 2) * 1) + 2 * 1 * 1 ≤
      ∑ _j : Unit, 1 * ((3 : ℝ) ^ 2 - 1)} = univ := eq_univ_of_forall fun x => by norm_num
  rw [this, probReal_univ, not_le]
  exact (exp_lt_exp.mpr (by norm_num)).trans_eq exp_zero

/-! ### Laurent–Massart, lower tail -/

/-- Positive control: every hypothesis of `laurent_massart_lower` holds for one `N(0,1)`. -/
example : (gaussianReal 0 1).real {ω : ℝ | ∑ _j : Unit, 1 * (ω ^ 2 - 1) ≤
    -(2 * √((∑ _j : Unit, (1 : ℝ) ^ 2) * 1))} ≤ exp (-1) :=
  laurent_massart_lower (fun _ ω => ω) (fun _ => measurable_id) iIndepFun.of_subsingleton
    (fun _ => map_id_gaussian) (fun _ => 1) (fun _ => zero_le_one) (by simp) one_pos

/-- Negative control, independence dropped: 64 copies of one `N(0,1)`, `λ = 1`, `τ = 9`. The
event is `64 (Z² − 1) ≤ −48`, i.e. `Z² ≤ 1/4`, which contains `[0, 1/2]`. -/
theorem lm_lower_needs_indep :
    ¬ ((gaussianReal 0 1).real {ω : ℝ | ∑ _j : Fin 64, 1 * (ω ^ 2 - 1) ≤
      -(2 * √((∑ _j : Fin 64, (1 : ℝ) ^ 2) * 9))} ≤ exp (-9)) := by
  have hsq : √((∑ _j : Fin 64, (1 : ℝ) ^ 2) * 9) = 24 := by
    rw [show (∑ _j : Fin 64, (1 : ℝ) ^ 2) * 9 = 24 ^ 2 by simp; norm_num]
    exact Real.sqrt_sq (by norm_num)
  rw [not_le]
  refine (exp_lt_mass (a := 0) (b := 1 / 2) ?_).trans_le
    ((gaussian_Icc_ge le_rfl (by norm_num)).trans (measureReal_mono ?_))
  · have := add_one_le_exp (9 - (1 / 2) ^ 2 / 2); norm_num at this ⊢; linarith
  · intro x hx
    rw [Set.mem_ofPred_eq, hsq]
    simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
    norm_num
    nlinarith [hx.1, hx.2]

/-- Negative control, `λ ≥ 0` dropped: `λ = −1`, `τ = 16`. The event is `Z² ≥ 9`, which
contains `[3, 4]`, of mass `> e^{−16}`. -/
theorem lm_lower_needs_nonneg :
    ¬ ((gaussianReal 0 1).real {ω : ℝ | ∑ _j : Unit, (-1) * (ω ^ 2 - 1) ≤
      -(2 * √((∑ _j : Unit, (-1 : ℝ) ^ 2) * 16))} ≤ exp (-16)) := by
  have hsq : √((∑ _j : Unit, (-1 : ℝ) ^ 2) * 16) = 4 := by
    rw [show (∑ _j : Unit, (-1 : ℝ) ^ 2) * 16 = 4 ^ 2 by simp; norm_num]
    exact Real.sqrt_sq (by norm_num)
  rw [not_le]
  refine (exp_lt_mass (a := 3) (b := 4) ?_).trans_le
    ((gaussian_Icc_ge (by norm_num) (by norm_num)).trans (measureReal_mono ?_))
  · have := add_one_le_exp (16 - 4 ^ 2 / 2); norm_num at this ⊢; linarith
  · intro x hx
    rw [Set.mem_ofPred_eq, hsq]
    simp only [Finset.univ_unique, Finset.sum_singleton]
    nlinarith [hx.1, hx.2]

/-- Negative control, `v > 0` dropped: `λ = 0`; the event `0 ≤ 0` is everything. -/
theorem lm_lower_needs_v_pos :
    ¬ ((gaussianReal 0 1).real {ω : ℝ | ∑ _j : Unit, 0 * (ω ^ 2 - 1) ≤
      -(2 * √((∑ _j : Unit, (0 : ℝ) ^ 2) * 1))} ≤ exp (-1)) := by
  have : {ω : ℝ | ∑ _j : Unit, 0 * (ω ^ 2 - 1) ≤
      -(2 * √((∑ _j : Unit, (0 : ℝ) ^ 2) * 1))} = univ := eq_univ_of_forall fun x => by simp
  rw [this, probReal_univ, not_le]
  exact (exp_lt_exp.mpr (by norm_num)).trans_eq exp_zero

/-- Negative control, the `N(0,1)` law dropped: `X ≡ 0`, `λ = 1`, `τ = 1/4`; the event
`−1 ≤ −1` is everything. -/
theorem lm_lower_needs_law :
    ¬ ((gaussianReal 0 1).real {_ω : ℝ | ∑ _j : Unit, 1 * ((0 : ℝ) ^ 2 - 1) ≤
      -(2 * √((∑ _j : Unit, (1 : ℝ) ^ 2) * (1 / 4)))} ≤ exp (-(1 / 4))) := by
  have hsq : √((∑ _j : Unit, (1 : ℝ) ^ 2) * (1 / 4)) = 1 / 2 := by
    rw [show (∑ _j : Unit, (1 : ℝ) ^ 2) * (1 / 4) = (1 / 2) ^ 2 by simp; norm_num]
    exact Real.sqrt_sq (by norm_num)
  have : {_ω : ℝ | ∑ _j : Unit, 1 * ((0 : ℝ) ^ 2 - 1) ≤
      -(2 * √((∑ _j : Unit, (1 : ℝ) ^ 2) * (1 / 4)))} = univ :=
    eq_univ_of_forall fun x => by rw [Set.mem_ofPred_eq, hsq]; norm_num
  rw [this, probReal_univ, not_le]
  exact (exp_lt_exp.mpr (by norm_num)).trans_eq exp_zero

/-- `τ > 0` is not needed for the bounds to be true: for `τ ≤ 0`, `e^{−τ} ≥ 1` bounds any
probability. It is kept because the proof divides by `√τ`. -/
theorem tail_trivial_of_nonpos {Ω : Type*} [MeasurableSpace Ω] (P : Measure Ω)
    [IsProbabilityMeasure P] (S : Set Ω) {τ : ℝ} (hτ : τ ≤ 0) : P.real S ≤ exp (-τ) :=
  (measureReal_le_one).trans (Real.one_le_exp (by linarith))

/-- Measurability is implied by the law hypothesis: the image of a probability measure under a
function that is not a.e. measurable is a Dirac mass, never `N(0,1)` (which has no atoms). -/
theorem law_forces_aemeasurable {Ω : Type*} [MeasurableSpace Ω] (P : Measure Ω)
    [IsProbabilityMeasure P] (f : Ω → ℝ) (h : P.map f = gaussianReal 0 1) : AEMeasurable f P := by
  by_contra hf
  rw [Measure.map_of_not_aemeasurable_of_ne_zero hf (IsProbabilityMeasure.ne_zero P)] at h
  have h1 : gaussianReal 0 1 {Classical.ofNonempty} = 0 := by
    rw [gaussianReal_apply 0 one_ne_zero]
    exact setLIntegral_measure_zero _ _ (Real.volume_singleton)
  have := congrArg (fun μ : Measure ℝ => μ {Classical.ofNonempty}) h
  simp only [Measure.dirac_apply_of_mem (Set.mem_singleton _), h1] at this
  exact one_ne_zero this

/-! ### Rotation invariance -/

/-- `c · Z` with `Z ~ N(0,1)` is `N(0,1)` only if `c² = 1` (its variance is `c²`). -/
theorem scaled_ne_std {c : ℝ} (hc : c ^ 2 ≠ 1) :
    (gaussianReal 0 1).map (fun ω : ℝ => c * ω) ≠ gaussianReal 0 1 := by
  intro h
  have hv := variance_id_map (μ := gaussianReal 0 1) (X := fun ω : ℝ => c * ω)
    (measurable_const.mul measurable_id).aemeasurable
  rw [h, variance_id_gaussianReal, variance_const_mul, variance_fun_id_gaussianReal] at hv
  exact hc (by simpa using hv.symm)

/-- Two independent `N(0,1)` coordinates on `ℝ²`. -/
noncomputable abbrev P2 : Measure (Fin 2 → ℝ) := Measure.pi fun _ => gaussianReal 0 1

theorem P2_law (i : Fin 2) : P2.map (fun ω => ω i) = gaussianReal 0 1 := by
  rw [show (fun ω : Fin 2 → ℝ => ω i) = Function.eval i from rfl, Measure.pi_map_eval]
  simp

theorem P2_indep : iIndepFun (fun (i : Fin 2) (ω : Fin 2 → ℝ) => ω i) P2 :=
  iIndepFun_pi (X := fun _ => id) fun _ => aemeasurable_id

/-- A rotation with rational entries `(3/5, 4/5; −4/5, 3/5)`. -/
noncomputable def rot : Fin 2 → Fin 2 → ℝ := ![![3 / 5, 4 / 5], ![-4 / 5, 3 / 5]]

theorem rot_orth : ∀ j k, ∑ i, rot j i * rot k i = if j = k then 1 else 0 := by
  intro j k; fin_cases j <;> fin_cases k <;> simp [rot, Fin.sum_univ_two] <;> norm_num

/-- Positive control: `rotation_iid` applied to two independent `N(0,1)` and `rot`. -/
example : iIndepFun (fun j ω => ∑ i, rot j i * (fun (i : Fin 2) (ω : Fin 2 → ℝ) => ω i) i ω) P2 :=
  (rotation_iid _ (fun i => measurable_pi_apply i) P2_indep P2_law rot rot_orth).2.1

/-- Negative control, orthonormality dropped: `u = 2` (one coordinate) gives `2Z ~ N(0,4)`. -/
theorem rotation_needs_orth :
    (gaussianReal 0 1).map (fun ω : ℝ => ∑ _i : Unit, (2 : ℝ) * ω) ≠ gaussianReal 0 1 := by
  simpa using scaled_ne_std (c := 2) (by norm_num)

/-- Negative control, independence dropped: `X₀ = X₁ = Z`; row 0 of `rot` gives `(7/5) Z`. -/
theorem rotation_needs_indep :
    (gaussianReal 0 1).map (fun ω : ℝ => ∑ i, rot 0 i * ω) ≠ gaussianReal 0 1 := by
  have : (fun ω : ℝ => ∑ i, rot 0 i * ω) = fun ω => 7 / 5 * ω := by
    funext ω; simp [rot, Fin.sum_univ_two]; ring
  rw [this]; exact scaled_ne_std (by norm_num)

/-- Negative control, the `N(0,1)` law dropped: `X ≡ 0` gives `Y ≡ 0 = 0 · Z`. -/
theorem rotation_needs_law :
    (gaussianReal 0 1).map (fun _ω : ℝ => ∑ i, rot 0 i * 0) ≠ gaussianReal 0 1 := by
  have : (fun _ω : ℝ => ∑ i, rot 0 i * 0) = fun ω => 0 * ω := by funext ω; simp
  rw [this]; exact scaled_ne_std (by norm_num)

/-- Positive control: `quadForm_eigen` for `B = 1` (`B` symmetric is a typing dependency:
the eigenvalues are only defined for symmetric `B`, so it has no negative control). -/
example (x : Fin 2 → ℝ) := quadForm_eigen (1 : Matrix (Fin 2) (Fin 2) ℝ) isHermitian_one x

/-! ### Null-gain tail -/

/-- Weights `λ = (1, −1)`. -/
noncomputable def lam2 : Fin 2 → ℝ := ![1, -1]

/-- Positive control: `gain_tail_eigen` with `λ = (1, −1)`, `a = 1`, `τ = 1`. -/
example : P2.real {ω | ∑ j, lam2 j + 2 * √((∑ j, max (lam2 j) 0 ^ 2) * 1) +
    2 * 1 * 1 + 2 * √((∑ j, max (-lam2 j) 0 ^ 2) * 1) ≤
    ∑ j, lam2 j * ω j ^ 2} ≤ 2 * exp (-1) :=
  gain_tail_eigen (fun (j : Fin 2) (ω : Fin 2 → ℝ) => ω j) (fun i => measurable_pi_apply i)
    P2_indep P2_law lam2 (fun j => by fin_cases j <;> simp [lam2])
    (by simp [lam2, Fin.sum_univ_two]) (by simp [lam2, Fin.sum_univ_two]) one_pos

/-- Negative control, `v⁺ > 0` and `v⁻ > 0` dropped together: `λ = 0`, `a = 0`, `τ = 1`; the
event `0 ≤ 0` is everything and `1 > 2e^{−1}`. (Each alone is not needed for truth: if one
side is empty, the other Laurent–Massart bound alone gives `≤ e^{−τ}`.) -/
theorem tail_needs_some_weight :
    ¬ (P2.real {ω | ∑ _j : Fin 2, (0 : ℝ) + 2 * √((∑ _j : Fin 2, max (0 : ℝ) 0 ^ 2) * 1) +
      2 * 0 * 1 + 2 * √((∑ _j : Fin 2, max (-0 : ℝ) 0 ^ 2) * 1) ≤
      ∑ j : Fin 2, (0 : ℝ) * ω j ^ 2} ≤ 2 * exp (-1)) := by
  have : {ω : Fin 2 → ℝ | ∑ _j : Fin 2, (0 : ℝ) + 2 * √((∑ _j : Fin 2, max (0 : ℝ) 0 ^ 2) * 1) +
      2 * 0 * 1 + 2 * √((∑ _j : Fin 2, max (-0 : ℝ) 0 ^ 2) * 1) ≤
      ∑ j : Fin 2, (0 : ℝ) * ω j ^ 2} = univ := eq_univ_of_forall fun x => by simp
  rw [this, probReal_univ, not_le]
  have h := add_one_le_exp 1
  have : exp (-1) * exp 1 = 1 := by rw [← exp_add]; simp
  nlinarith [exp_pos (-1), quadratic_le_exp_of_nonneg (x := 1) zero_le_one]

/-- Mass of the unit box `[0,1]²` under `P2` exceeds `e^{−10}`. -/
theorem P2_box_mass : exp (-10) < P2.real (univ.pi fun _ => Icc 0 1) := by
  have hm : exp (-5) < (gaussianReal 0 1).real (Icc 0 1) :=
    (exp_lt_mass (a := 0) (b := 1) (by
      have := add_one_le_exp (5 - 1 ^ 2 / 2); norm_num at this ⊢; linarith)).trans_le
      (gaussian_Icc_ge le_rfl zero_le_one)
  have hP : P2.real (univ.pi fun _ => Icc 0 1) = (gaussianReal 0 1).real (Icc 0 1) ^ 2 := by
    simp only [measureReal_def, Measure.pi_pi, Fin.prod_univ_two, ENNReal.toReal_mul]; ring
  rw [hP, show exp (-10) = exp (-5) ^ 2 by rw [← exp_nat_mul]; norm_num]
  gcongr

/-- Negative control, `max(λ_j, 0) ≤ a` dropped: `λ = (1, −1)`, `a = −100`, `τ = 12`. The
threshold is `4√12 − 2400 < −1 ≤ ω₀² − ω₁²` on `[0,1]²`, of mass `> e^{−11} > 2e^{−12}`. -/
theorem tail_needs_le_a :
    ¬ (P2.real {ω | ∑ j, lam2 j + 2 * √((∑ j, max (lam2 j) 0 ^ 2) * 12) +
      2 * (-100) * 12 + 2 * √((∑ j, max (-lam2 j) 0 ^ 2) * 12) ≤
      ∑ j, lam2 j * ω j ^ 2} ≤ 2 * exp (-12)) := by
  have hm : exp (-11) < (gaussianReal 0 1).real (Icc 0 1) ^ 2 := by
    have hm : exp (-5.5) < (gaussianReal 0 1).real (Icc 0 1) :=
      (exp_lt_mass (a := 0) (b := 1) (by
        have := add_one_le_exp (5.5 - 1 ^ 2 / 2); norm_num at this ⊢; linarith)).trans_le
        (gaussian_Icc_ge le_rfl zero_le_one)
    rw [show exp (-11) = exp (-5.5) ^ 2 by rw [← exp_nat_mul]; norm_num]
    gcongr
  have hP : P2.real (univ.pi fun _ => Icc 0 1) = (gaussianReal 0 1).real (Icc 0 1) ^ 2 := by
    simp only [measureReal_def, Measure.pi_pi, Fin.prod_univ_two, ENNReal.toReal_mul]; ring
  have h2 : 2 * exp (-12) ≤ exp (-11) := by
    have h := add_one_le_exp 1
    rw [show (-11 : ℝ) = 1 + -12 by norm_num, exp_add]
    nlinarith [exp_pos (-12)]
  have hs : √12 < 4 := by rw [Real.sqrt_lt' (by norm_num)]; norm_num
  rw [not_le]
  refine h2.trans_lt (hm.trans_le ?_)
  rw [← hP]
  refine measureReal_mono fun x hx => ?_
  rw [Set.mem_univ_pi] at hx
  have h0 := hx 0; have h1 := hx 1
  simp only [Set.mem_ofPred_eq, Fin.sum_univ_two, lam2, Matrix.cons_val_zero,
    Matrix.cons_val_one, Matrix.cons_val_fin_one]
  norm_num
  nlinarith [h0.1, h0.2, h1.1, h1.2]

/-- Negative control, the `N(0,1)` law dropped: `X₀ ≡ 3`, `X₁ ≡ 0`, `λ = (1, −1)`, `a = 1`,
`τ = 1`. The event `6 ≤ 9` is everything and `1 > 2e^{−1}`. -/
theorem tail_needs_law :
    ¬ (P2.real {_ω | ∑ j, lam2 j + 2 * √((∑ j, max (lam2 j) 0 ^ 2) * 1) +
      2 * 1 * 1 + 2 * √((∑ j, max (-lam2 j) 0 ^ 2) * 1) ≤
      ∑ j, lam2 j * ![(3 : ℝ), 0] j ^ 2} ≤ 2 * exp (-1)) := by
  have : {_ω : Fin 2 → ℝ | ∑ j, lam2 j +
      2 * √((∑ j, max (lam2 j) 0 ^ 2) * 1) +
      2 * 1 * 1 + 2 * √((∑ j, max (-lam2 j) 0 ^ 2) * 1) ≤
      ∑ j, lam2 j * ![(3 : ℝ), 0] j ^ 2} = univ :=
    eq_univ_of_forall fun x => by simp [lam2, Fin.sum_univ_two]; norm_num
  rw [this, probReal_univ, not_le]
  nlinarith [exp_pos (-1), quadratic_le_exp_of_nonneg (x := 1) zero_le_one,
    show exp (-1) * exp 1 = 1 by rw [← exp_add]; simp]

/-- Five weights `(1, 1, 1, 1, −1)` on the coordinates `(ω₀, ω₀, ω₀, ω₀, ω₁)`: each is
`N(0,1)` and measurable, but the first four are the same variable. -/
noncomputable def lam5 : Fin 5 → ℝ := fun j => if j = 4 then -1 else 1

def X5 : Fin 5 → (Fin 2 → ℝ) → ℝ := fun j ω => if j = 4 then ω 1 else ω 0

/-- Negative control, independence dropped: weights `lam5` on `X5`, `a = 1`, `τ = 25`. The
threshold is `3 + 20 + 50 + 10 = 83` and the gain is `4ω₀² − ω₁²`, which is `≥ 87` on
`[4.7, 5.7] × [0, 1]`; that box has mass `> e^{−19} e^{−5} ≥ 2e^{−25}`. -/
theorem tail_needs_indep :
    ¬ (P2.real {ω | ∑ j, lam5 j + 2 * √((∑ j, max (lam5 j) 0 ^ 2) * 25) + 2 * 1 * 25 +
      2 * √((∑ j, max (-lam5 j) 0 ^ 2) * 25) ≤ ∑ j, lam5 j * X5 j ω ^ 2} ≤ 2 * exp (-25)) := by
  have e1 : ∑ j, lam5 j = 3 := by simp [lam5, Fin.sum_univ_five]; norm_num
  have e2 : √((∑ j, max (lam5 j) 0 ^ 2) * 25) = 10 := by
    rw [show (∑ j, max (lam5 j) 0 ^ 2) * 25 = 10 ^ 2 by simp [lam5, Fin.sum_univ_five]; norm_num]
    exact Real.sqrt_sq (by norm_num)
  have e3 : √((∑ j, max (-lam5 j) 0 ^ 2) * 25) = 5 := by
    rw [show (∑ j, max (-lam5 j) 0 ^ 2) * 25 = 5 ^ 2 by simp [lam5, Fin.sum_univ_five]; norm_num]
    exact Real.sqrt_sq (by norm_num)
  have e4 : ∀ ω : Fin 2 → ℝ, ∑ j, lam5 j * X5 j ω ^ 2 = 4 * ω 0 ^ 2 - ω 1 ^ 2 := by
    intro ω; simp [lam5, X5, Fin.sum_univ_five]; ring
  have m0 : exp (-19) < (gaussianReal 0 1).real (Icc 4.7 5.7) :=
    (exp_lt_mass (a := 4.7) (b := 5.7) (by
      have := add_one_le_exp (19 - 5.7 ^ 2 / 2); norm_num at this ⊢; linarith)).trans_le
      (gaussian_Icc_ge (by norm_num) (by norm_num))
  have m1 : exp (-5) < (gaussianReal 0 1).real (Icc 0 1) :=
    (exp_lt_mass (a := 0) (b := 1) (by
      have := add_one_le_exp (5 - 1 ^ 2 / 2); norm_num at this ⊢; linarith)).trans_le
      (gaussian_Icc_ge le_rfl zero_le_one)
  have hP : P2.real (univ.pi ![Icc (4.7 : ℝ) 5.7, Icc 0 1]) =
      (gaussianReal 0 1).real (Icc 4.7 5.7) * (gaussianReal 0 1).real (Icc 0 1) := by
    simp only [measureReal_def, Measure.pi_pi, Fin.prod_univ_two, ENNReal.toReal_mul]
    rfl
  have h2 : 2 * exp (-25) ≤ exp (-19) * exp (-5) := by
    have h := add_one_le_exp 1
    rw [← exp_add, show (-19 + -5 : ℝ) = 1 + -25 by norm_num, exp_add]
    nlinarith [exp_pos (-25)]
  rw [not_le]
  refine h2.trans_lt ((mul_lt_mul'' m0 m1 (exp_pos _).le (exp_pos _).le).trans_le ?_)
  rw [← hP]
  refine measureReal_mono fun x hx => ?_
  rw [Set.mem_univ_pi] at hx
  have h0 : x 0 ∈ Icc (4.7 : ℝ) 5.7 := hx 0
  have h1 : x 1 ∈ Icc (0 : ℝ) 1 := hx 1
  rw [Set.mem_ofPred_eq, e1, e2, e3, e4]
  nlinarith [h0.1, h0.2, h1.1, h1.2]

/-- Positive control for `gain_null_tail`: `L = 1`, `A = diag(1, −1)` on `ℝ²`, `τ = 1`. All
hypotheses are discharged; `det A = −1 = λ₀ λ₁` forces one eigenvalue of each sign. -/
example : True := by
  set A : Matrix (Fin 2) (Fin 2) ℝ := diagonal ![1, -1]
  have hA : Aᵀ = A := diagonal_transpose _
  have hB := sandwich_isHermitian 1 A hA
  have hdet : hB.eigenvalues 0 * hB.eigenvalues 1 = -1 := by
    have hd : (1ᵀ * A * 1).det = -1 := by simp [A, det_diagonal, Fin.prod_univ_two]
    have h := hB.det_eq_prod_eigenvalues
    rw [hd, Fin.prod_univ_two] at h
    exact h.symm
  have hp : 0 < ∑ j, max (hB.eigenvalues j) 0 ^ 2 := by
    rw [Fin.sum_univ_two]
    rcases lt_or_ge 0 (hB.eigenvalues 0) with h | h
    · have := max_eq_left h.le; nlinarith [sq_nonneg (max (hB.eigenvalues 1) 0)]
    · have h1 : 0 < hB.eigenvalues 1 := by nlinarith
      have := max_eq_left h1.le; nlinarith [sq_nonneg (max (hB.eigenvalues 0) 0)]
  have hn : 0 < ∑ j, max (-hB.eigenvalues j) 0 ^ 2 := by
    rw [Fin.sum_univ_two]
    rcases lt_or_ge (hB.eigenvalues 0) 0 with h | h
    · have := max_eq_left (neg_pos.mpr h).le; nlinarith [sq_nonneg (max (-hB.eigenvalues 1) 0)]
    · have h1 : hB.eigenvalues 1 < 0 := by nlinarith
      have := max_eq_left (neg_pos.mpr h1).le
      nlinarith [sq_nonneg (max (-hB.eigenvalues 0) 0)]
  have ha : ∀ j, max (hB.eigenvalues j) 0 ≤
      max (hB.eigenvalues 0) 0 + max (hB.eigenvalues 1) 0 := by
    intro j; fin_cases j
    · simp
    · simp
  have := gain_null_tail (fun (j : Fin 2) (ω : Fin 2 → ℝ) => ω j) (fun i => measurable_pi_apply i) P2_indep P2_law
    1 A hB ha hp hn one_pos
  trivial

end TTMT.Controls3
