# Formal proofs (Lean 4 + Mathlib)

Machine-checked versions of the deterministic cores of results in Appendix AB and in the
theory of the M6 test. Lean `v4.35.0`, Mathlib tag `v4.35.0` (see `lean-toolchain`,
`lakefile.toml`). CI: `.github/workflows/lean.yml` builds the library and fails on any
`sorry`, `admit` or `axiom`; `formal/AxiomCheck.lean` prints the axioms each theorem uses
(only `propext`, `Classical.choice`, `Quot.sound`, the standard ones of Mathlib).

Build locally: `cd formal && lake exe cache get && lake build`.

## Status labels (rule 31)

A Lean proof is strong verification of the **formal statement**. It does not replace the
independent review nor the director's review of the results that support the contribution,
and it covers only what the formal statement says. Each row lists what the formal statement
leaves out or changes with respect to the written one.

| Result | Lean file / theorem | What is formalized | Left out / changed |
|:--|:--|:--|:--|
| AB.A.9 (P carries no first-order term; F3) | `PPolarization.lean`: `pAmplitude_cross_eq_zero`, `kappa_P_eq_zero` | For any fixed vector `u`, any `n_src`, `c` and any two 3x3 tensors, the cross product of the two P amplitudes is 0, so `κ_iP = 0` | The far-field form of the operator is an input (definition), not derived from elastodynamics |
| E.2, deterministic part (null gain of M6) | `M6NullGain.lean`: `gain_null_eq`, `gain_null_pivotal` | With `P₁`, `P₂` orthogonal projectors and the signal in both classes, `RSS₁(μ+ε) − RSS₂(μ+ε) = εᵀ(P₂−P₁)ε` for every ε: the null gain does not depend on the event | Fixed times only (the search inflation of E.2 is not covered, as in E.2). Mean `tr(C P₂) − tr(C P₁)` and the Laurent–Massart tail not formalized |
| AB.C.1, deterministic part (κ → 1) | `ReadingContraction.lean`: `reading_contraction`, `lipschitz_comp_one`, `lipschitz_sup` | The one-sample-per-block reading with rows of norm ≤ 1 and disjoint supports is a contraction; composition of 1-Lipschitz maps and the finite maximum of 1-Lipschitz functions are 1-Lipschitz | The 1-Lipschitz property of the radial map and of `Q_hᵀ` (Lemma AA.9.1) enters as a hypothesis; Gaussian concentration and `E F ≤ √6 + √(2 n_b log L)` are cited, not formalized (not in Mathlib) |
| AB.A.2 (coincident events are not identifiable) | `Superposition.lean`: `coincident_pair_mem_oneEventClass` | If `G(θ)` is linear in the tensor, the pair signal with `θ₁ = θ₂` lies in the one-event class | Margin 0 follows from membership; the testing consequence (Theorem AB.A.1) not formalized |

Not yet formalized (candidates, in priority order): AB.A.3(a) station decomposition of the
timing-tolerant margin; AB.C.3 (bounded influence implies saturation); AB.C.2 (conditional
tail; needs rotation-invariant laws on the sphere); AB.A.1 (Le Cam bound; needs total
variation between Gaussians).
