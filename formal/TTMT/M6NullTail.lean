/-
Tail of the M6 null gain (Proposition E.2, second claim; teoria-fase5/teoria-m6-prueba.md),
in the eigen-coordinates of the gain.

`gain_tail_eigen`: let `X_j` be independent standard Gaussians and `λ_j` real weights of
either sign, `v⁺ = ∑ (λ_j⁺)²`, `v⁻ = ∑ (λ_j⁻)²` (both positive), `a ≥ λ_j` for all `j`, and
`m₀ = ∑ λ_j`. Then for every `τ > 0`
  `P(∑ λ_j X_j² ≥ m₀ + 2√(v⁺τ) + 2aτ + 2√(v⁻τ)) ≤ 2e^{−τ}`,
which is the bound of E.2 (upper Laurent–Massart tail of the positive part, lower tail of the
negative part, and a union bound).

The step from `εᵀ(P₂−P₁)ε` with `ε ~ N(0, C)` to the eigen-coordinates is in
`GaussianRotation.lean`.
-/
import TTMT.LaurentMassart

open MeasureTheory ProbabilityTheory Real Finset

namespace TTMT

variable {Ω ι : Type*} [MeasurableSpace Ω] {P : Measure Ω} [IsProbabilityMeasure P] [Fintype ι]

theorem gain_tail_eigen (X : ι → Ω → ℝ) (hmeas : ∀ j, Measurable (X j))
    (hind : iIndepFun X P) (hlaw : ∀ j, P.map (X j) = gaussianReal 0 1) (lam : ι → ℝ)
    {a : ℝ} (ha : ∀ j, max (lam j) 0 ≤ a)
    (hvp : 0 < ∑ j, max (lam j) 0 ^ 2) (hvn : 0 < ∑ j, max (-lam j) 0 ^ 2)
    {τ : ℝ} (hτ : 0 < τ) :
    P.real {ω | ∑ j, lam j + 2 * √((∑ j, max (lam j) 0 ^ 2) * τ) + 2 * a * τ +
        2 * √((∑ j, max (-lam j) 0 ^ 2) * τ) ≤ ∑ j, lam j * X j ω ^ 2} ≤ 2 * exp (-τ) := by
  set lp : ι → ℝ := fun j => max (lam j) 0
  set ln : ι → ℝ := fun j => max (-lam j) 0
  have hsplit : ∀ j, lam j = lp j - ln j := by
    intro j; simp only [lp, ln]
    rcases le_total (lam j) 0 with h | h
    · rw [max_eq_right h, max_eq_left (by linarith)]; ring
    · rw [max_eq_left h, max_eq_right (by linarith)]; ring
  have hU := laurent_massart_upper X hmeas hind hlaw lp (fun j => le_max_right _ _) ha hvp hτ
  have hL := laurent_massart_lower X hmeas hind hlaw ln (fun j => le_max_right _ _) hvn hτ
  have hsub : {ω | ∑ j, lam j + 2 * √((∑ j, lp j ^ 2) * τ) + 2 * a * τ +
        2 * √((∑ j, ln j ^ 2) * τ) ≤ ∑ j, lam j * X j ω ^ 2} ⊆
      {ω | 2 * √((∑ j, lp j ^ 2) * τ) + 2 * a * τ ≤ ∑ j, lp j * (X j ω ^ 2 - 1)} ∪
      {ω | ∑ j, ln j * (X j ω ^ 2 - 1) ≤ -(2 * √((∑ j, ln j ^ 2) * τ))} := by
    intro ω hω
    simp only [Set.mem_ofPred_eq, Set.mem_union] at hω ⊢
    by_contra h
    push Not at h
    obtain ⟨h1, h2⟩ := h
    have e : ∑ j, lam j * X j ω ^ 2 =
        ∑ j, lam j + (∑ j, lp j * (X j ω ^ 2 - 1) - ∑ j, ln j * (X j ω ^ 2 - 1)) := by
      rw [← sum_sub_distrib, ← sum_add_distrib]
      exact sum_congr rfl fun j _ => by rw [hsplit j]; ring
    linarith
  calc P.real _ ≤ P.real ({ω | 2 * √((∑ j, lp j ^ 2) * τ) + 2 * a * τ ≤
          ∑ j, lp j * (X j ω ^ 2 - 1)} ∪
        {ω | ∑ j, ln j * (X j ω ^ 2 - 1) ≤ -(2 * √((∑ j, ln j ^ 2) * τ))}) :=
        measureReal_mono hsub
    _ ≤ _ := measureReal_union_le _ _
    _ ≤ exp (-τ) + exp (-τ) := add_le_add hU hL
    _ = 2 * exp (-τ) := by ring

end TTMT
