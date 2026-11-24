/-
Radial maps of Appendix AA.9 (docs/proposal/Supplementary_Evidence_Appendices_G_to_AA.md,
§AA.9 and the proof of AA.9.5(b)) are 1-Lipschitz.

A radial map is `φ(x) = r(‖x‖) x / ‖x‖` (and `φ(0) = 0`). `radial_lipschitz` proves, in any
real inner product space, that `φ` is 1-Lipschitz whenever `r ≥ 0`, `r(t) ≤ t` and `r` is
1-Lipschitz on `[0, ∞)`. The written proof states the constant `max(sup r', sup r(t)/t)`;
the formal statement uses the two conditions it needs, which hold for the three maps used:
- the identity `r(t) = t`;
- the cap `r_c(t) = min(t, c)`, `c ≥ 0`;
- the logarithmic map `ℓ_c(t) = c √(log(1 + t²/c²))`, `c > 0`, via the mean value theorem
  on `(0, ∞)` and `ℓ_c'(t) ≤ 1` (from `log(1+u) ≥ u/(1+u)`), as in the written proof.
-/
import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Deriv
import Mathlib.Analysis.SpecialFunctions.Sqrt
import Mathlib.Analysis.Calculus.Deriv.MeanValue

open Real

namespace TTMT

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

/-- The radial map with radius function `r`: `x ↦ (r ‖x‖ / ‖x‖) • x` (zero at zero). -/
noncomputable def radial (r : ℝ → ℝ) (x : E) : E := (r ‖x‖ / ‖x‖) • x

/-- Conditions on the radius function under which the radial map is 1-Lipschitz. -/
structure GoodRadius (r : ℝ → ℝ) : Prop where
  nonneg : ∀ t, 0 ≤ t → 0 ≤ r t
  le_self : ∀ t, 0 ≤ t → r t ≤ t
  lip : ∀ a b, 0 ≤ a → 0 ≤ b → |r a - r b| ≤ |a - b|

theorem norm_radial {r : ℝ → ℝ} (hr : GoodRadius r) (x : E) : ‖radial r x‖ = r ‖x‖ := by
  unfold radial
  rcases eq_or_ne x 0 with rfl | hx
  · simp only [norm_zero, smul_zero]
    have := hr.le_self 0 le_rfl; have := hr.nonneg 0 le_rfl; linarith
  · have hn : 0 < ‖x‖ := norm_pos_iff.mpr hx
    rw [norm_smul, Real.norm_eq_abs, abs_of_nonneg (div_nonneg (hr.nonneg _ hn.le) hn.le),
      div_mul_cancel₀ _ hn.ne']

/-- **Radial maps are 1-Lipschitz** (the step "the block-wise radial map is 1-Lipschitz" of
AA.9.5(b) and AB.C.1). -/
theorem radial_lipschitz {r : ℝ → ℝ} (hr : GoodRadius r) (x y : E) :
    ‖radial r x - radial r y‖ ≤ ‖x - y‖ := by
  have h0 : radial r (0 : E) = 0 := by simp [radial]
  rcases eq_or_ne x 0 with rfl | hx
  · rw [h0, zero_sub, norm_neg, zero_sub, norm_neg, norm_radial hr]
    exact hr.le_self _ (norm_nonneg _)
  rcases eq_or_ne y 0 with rfl | hy
  · rw [h0, sub_zero, sub_zero, norm_radial hr]
    exact hr.le_self _ (norm_nonneg _)
  have ha : 0 < ‖x‖ := norm_pos_iff.mpr hx
  have hb : 0 < ‖y‖ := norm_pos_iff.mpr hy
  set a := ‖x‖ with ha_def
  set b := ‖y‖ with hb_def
  set s := r a / a
  set t := r b / b
  set p := @inner ℝ E _ x y
  have hsa : s * a = r a := div_mul_cancel₀ _ ha.ne'
  have htb : t * b = r b := div_mul_cancel₀ _ hb.ne'
  have hs0 : 0 ≤ s := div_nonneg (hr.nonneg _ ha.le) ha.le
  have ht0 : 0 ≤ t := div_nonneg (hr.nonneg _ hb.le) hb.le
  have hs1 : s ≤ 1 := (div_le_one ha).mpr (hr.le_self _ ha.le)
  have ht1 : t ≤ 1 := (div_le_one hb).mpr (hr.le_self _ hb.le)
  have hp : p ≤ a * b := real_inner_le_norm x y
  have hlip : (r a - r b) ^ 2 ≤ (a - b) ^ 2 := by
    have := hr.lip a b ha.le hb.le
    nlinarith [abs_nonneg (r a - r b), sq_abs (r a - r b), sq_abs (a - b)]
  have key : ‖s • x - t • y‖ ^ 2 ≤ ‖x - y‖ ^ 2 := by
    rw [norm_sub_sq_real, norm_sub_sq_real, norm_smul, norm_smul, real_inner_smul_left,
      real_inner_smul_right, Real.norm_eq_abs, Real.norm_eq_abs, abs_of_nonneg hs0,
      abs_of_nonneg ht0]
    rw [← ha_def, ← hb_def, hsa, htb]
    have hst : 0 ≤ 1 - s * t := by nlinarith
    have h1 := mul_nonneg hst (sub_nonneg.mpr hp)
    have hprod : r a * r b = s * t * (a * b) := by rw [← hsa, ← htb]; ring
    nlinarith [h1, hprod, hlip]
  exact (pow_le_pow_iff_left₀ (norm_nonneg _) (norm_nonneg _) two_ne_zero).mp key

/-- The identity radius. -/
theorem goodRadius_id : GoodRadius id :=
  ⟨fun _ h => h, fun _ _ => le_rfl, fun _ _ _ _ => le_rfl⟩

/-- The cap `r_c(t) = min(t, c)`. -/
theorem goodRadius_cap {c : ℝ} (hc : 0 ≤ c) : GoodRadius fun t => min t c :=
  ⟨fun _ h => le_min h hc, fun _ _ => min_le_left _ _, fun a b _ _ => by
    simpa using abs_min_sub_min_le_max a c b c⟩

/-- The logarithmic radius `ℓ_c(t) = c √(log(1 + t²/c²))`. -/
noncomputable def logRadius (c : ℝ) (t : ℝ) : ℝ := c * √(log (1 + t ^ 2 / c ^ 2))

theorem logRadius_le_self {c : ℝ} (hc : 0 < c) {t : ℝ} (ht : 0 ≤ t) : logRadius c t ≤ t := by
  unfold logRadius
  have hu : 0 ≤ t ^ 2 / c ^ 2 := by positivity
  have hlog : log (1 + t ^ 2 / c ^ 2) ≤ (t / c) ^ 2 := by
    have := log_le_sub_one_of_pos (by linarith : (0 : ℝ) < 1 + t ^ 2 / c ^ 2)
    rw [div_pow]; linarith
  have : √(log (1 + t ^ 2 / c ^ 2)) ≤ t / c := by
    calc √(log (1 + t ^ 2 / c ^ 2)) ≤ √((t / c) ^ 2) := Real.sqrt_le_sqrt hlog
      _ = t / c := Real.sqrt_sq (div_nonneg ht hc.le)
  calc c * √(log (1 + t ^ 2 / c ^ 2)) ≤ c * (t / c) := mul_le_mul_of_nonneg_left this hc.le
    _ = t := mul_div_cancel₀ t hc.ne'

theorem logRadius_mono {c : ℝ} (hc : 0 < c) {a b : ℝ} (ha : 0 ≤ a) (hab : a ≤ b) :
    logRadius c a ≤ logRadius c b := by
  unfold logRadius
  refine mul_le_mul_of_nonneg_left (Real.sqrt_le_sqrt ?_) hc.le
  refine log_le_log (by positivity) ?_
  have : a ^ 2 ≤ b ^ 2 := pow_le_pow_left₀ ha hab 2
  have hc2 : 0 < c ^ 2 := by positivity
  linarith [div_le_div_of_nonneg_right this hc2.le]

theorem hasDerivAt_logRadius {c : ℝ} (hc : 0 < c) {t : ℝ} (ht : 0 < t) :
    HasDerivAt (logRadius c)
      (c * ((2 * t / c ^ 2) / (1 + t ^ 2 / c ^ 2) / (2 * √(log (1 + t ^ 2 / c ^ 2))))) t := by
  have hpos : 0 < 1 + t ^ 2 / c ^ 2 := by positivity
  have hu : 0 < t ^ 2 / c ^ 2 := by positivity
  have hL : 0 < log (1 + t ^ 2 / c ^ 2) := log_pos (by linarith)
  have h1 : HasDerivAt (fun t : ℝ => 1 + t ^ 2 / c ^ 2) (2 * t / c ^ 2) t := by
    have := ((hasDerivAt_pow 2 t).div_const (c ^ 2)).const_add 1
    convert this using 1; push_cast; ring
  have h2 := h1.log hpos.ne'
  have h3 := h2.sqrt hL.ne'
  exact h3.const_mul c

theorem deriv_logRadius_le_one {c : ℝ} (hc : 0 < c) {t : ℝ} (ht : 0 < t) :
    deriv (logRadius c) t ≤ 1 := by
  rw [(hasDerivAt_logRadius hc ht).deriv]
  set u := t ^ 2 / c ^ 2 with hu
  have hu0 : 0 < u := by positivity
  have hL : 0 < log (1 + u) := log_pos (by linarith)
  have hlog : u / (1 + u) ≤ log (1 + u) := by
    have := one_sub_inv_le_log_of_pos (by linarith : (0 : ℝ) < 1 + u)
    have : 1 - (1 + u)⁻¹ = u / (1 + u) := by field_simp; ring
    linarith
  set L := log (1 + u)
  have hsL : 0 < √L := Real.sqrt_pos.mpr hL
  have hsq : √L ^ 2 = L := Real.sq_sqrt hL.le
  have hcu : c ^ 2 * u = t ^ 2 := by rw [hu]; field_simp
  -- the derivative equals t / (c (1+u) √L)
  have hform : c * ((2 * t / c ^ 2) / (1 + u) / (2 * √L)) = t / (c * (1 + u) * √L) := by
    field_simp
  rw [hform, div_le_one (by positivity)]
  -- t ≤ c (1+u) √L  ⇐  t² ≤ c² (1+u)² L
  have hrhs : 0 ≤ c * (1 + u) * √L := by positivity
  have : t ^ 2 ≤ (c * (1 + u) * √L) ^ 2 := by
    rw [mul_pow, mul_pow, hsq]
    have h1 : u ≤ (1 + u) ^ 2 * L := by
      have : u / (1 + u) * (1 + u) ^ 2 ≤ L * (1 + u) ^ 2 :=
        mul_le_mul_of_nonneg_right hlog (by positivity)
      have e : u / (1 + u) * (1 + u) ^ 2 = u * (1 + u) := by field_simp
      nlinarith
    nlinarith [mul_le_mul_of_nonneg_left h1 (sq_nonneg c)]
  exact (pow_le_pow_iff_left₀ ht.le hrhs two_ne_zero).mp this

theorem goodRadius_log {c : ℝ} (hc : 0 < c) : GoodRadius (logRadius c) := by
  refine ⟨fun t _ => mul_nonneg hc.le (Real.sqrt_nonneg _), fun t ht => logRadius_le_self hc ht,
    ?_⟩
  -- 1-Lipschitz on [0, ∞): mean value theorem on the interior (0, ∞)
  have hcont : Continuous (logRadius c) := by
    unfold logRadius
    refine continuous_const.mul (Continuous.sqrt (Continuous.log (by fun_prop) fun x => ?_))
    have : 0 ≤ x ^ 2 / c ^ 2 := by positivity
    linarith
  have hmvt : ∀ a b, 0 ≤ a → a ≤ b → logRadius c b - logRadius c a ≤ 1 * (b - a) := by
    intro a b ha hab
    refine (convex_Ici (0 : ℝ)).image_sub_le_mul_sub_of_deriv_le hcont.continuousOn ?_ ?_ a ha
      b (ha.trans hab) hab
    · rw [interior_Ici]
      exact fun t ht => (hasDerivAt_logRadius hc ht).differentiableAt.differentiableWithinAt
    · rw [interior_Ici]; exact fun t ht => deriv_logRadius_le_one hc ht
  intro a b ha hb
  rcases le_total a b with hab | hab
  · have h1 := hmvt a b ha hab
    have h2 := logRadius_mono hc ha hab
    rw [abs_of_nonpos (by linarith), abs_of_nonpos (by linarith)]; linarith
  · have h1 := hmvt b a hb hab
    have h2 := logRadius_mono hc hb hab
    rw [abs_of_nonneg (by linarith), abs_of_nonneg (by linarith)]; linarith

end TTMT
