/-
Proposition AB.A.2 (superposition; docs/proposal/Supplementary_Evidence_Appendix_AB.md, §AB.A.2):
if the two events coincide (θ₁ = θ₂), the pair signal G(θ₁)m₁ + G(θ₁)m₂ = G(θ₁)(m₁ + m₂)
lies in the one-event class, so the margin δ = inf_{μ ∈ M₁} ‖μ_ξ − μ‖ is zero and no test
tells "one" from "two".

Formal hypothesis: G(θ) is linear in the tensor m for each fixed θ (the model of §AB.A.0).
Nothing else is used.
-/
import Mathlib.Basic.Real.Basic
import Mathlib.Algebra.Module.LinearMap.Defs

namespace TTMT

variable {Θ V W : Type*} [AddCommMonoid V] [AddCommMonoid W] [Module ℝ V] [Module ℝ W]

/-- One-event class: the signals `G θ m`. -/
def oneEventClass (G : Θ → V →ₗ[ℝ] W) : Set W := Set.range fun p : Θ × V => G p.1 p.2

/-- **Proposition AB.A.2.** With θ₁ = θ₂ the pair signal is in the one-event class
(zero margin against that class). -/
theorem coincident_pair_mem_oneEventClass (G : Θ → V →ₗ[ℝ] W) (θ : Θ) (m₁ m₂ : V) :
    G θ m₁ + G θ m₂ ∈ oneEventClass G :=
  ⟨(θ, m₁ + m₂), (G θ).map_add m₁ m₂⟩

end TTMT
