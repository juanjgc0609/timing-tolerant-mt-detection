# Formal proofs (Lean 4 + Mathlib)

Machine-checked versions of the deterministic cores of results in Appendix AB and in the
theory of the M6 test. Lean `v4.35.0`, Mathlib tag `v4.35.0` (see `lean-toolchain`,
`lakefile.toml`). CI: `.github/workflows/lean.yml` builds the library and fails on any
`sorry`, `admit` or `axiom`; `formal/AxiomCheck.lean` prints the axioms each theorem uses, and
CI fails unless every line lists exactly `propext`, `Classical.choice`, `Quot.sound` (the
standard ones of Mathlib).

Build locally: `cd formal && lake exe cache get && lake build`.

## Status labels (rule 31)

A Lean proof is strong verification of the **formal statement**. It does not replace the
independent review nor the director's review of the results that support the contribution,
and it covers only what the formal statement says. Each row lists what the formal statement
leaves out or changes with respect to the written one.

| Result | Lean file / theorem | What is formalized | Left out / changed |
|:--|:--|:--|:--|
| AB.A.9 (P carries no first-order term; F3) | `PPolarization.lean`: `pAmplitude_cross_eq_zero`, `kappa_P_eq_zero` | For any fixed vector `u`, any `n_src`, `c` and any two 3x3 tensors, the cross product of the two P amplitudes is 0, so `κ_iP = 0` | The far-field form of the operator is an input (definition), not derived from elastodynamics |
| E.2, null gain at fixed times | `M6NullGain.lean`: `gain_null_eq`, `gain_null_pivotal` | With `P₁`, `P₂` orthogonal projectors and the signal in both classes, `RSS₁(μ+ε) − RSS₂(μ+ε) = εᵀ(P₂−P₁)ε` for every ε | Fixed times only (the search inflation is not covered, as in E.2) |
| E.2, mean | `M6NullMean.lean`: `integral_quadForm`, `integral_gain_null`, `trace_mul_sum_vecMulVec`, `trace_eq_lambda_mul_card` | `E[gain] = tr((P₂−P₁)Cᵀ)` for any noise law with second moments `C`; `tr(C P) = ∑_k q_kᵀCq_k`, and `= λ_r dim P` when each basis vector gives `λ_r` | That each column gives exactly `λ_r = rᵀCr/rᵀr` (stationary noise, separated pulses away from the edges) is a hypothesis. The Laurent–Massart tail is not formalized |
| AB.C.1, deterministic part, complete | `CappedScore.lean`: `readingMap_lipschitz`, `readBlocks_lipschitz`, `blockRadial_lipschitz`, `tailFunction_lipschitz`, `read_le_tailFunction`; `ReadingContraction.lean`; `RadialMap.lean` | In Euclidean spaces: the reading (rows of norm ≤ 1, disjoint supports) is 1-Lipschitz from the white noise; the block-wise radial map is 1-Lipschitz for the identity, the cap and the logarithmic map `ℓ_c`; `Q†` is a contraction; so `F(z) = max_s ‖Q† φ(A_s z)‖` is 1-Lipschitz, and the score at any data-dependent reading is `≤ F` | The radial lemma uses `0 ≤ r(t) ≤ t` and `r` 1-Lipschitz instead of the written constant `max(sup r', sup r(t)/t)` (equivalent for the three maps used). Gaussian concentration `P(F ≥ EF + t) ≤ e^{-t²/2}` and `EF ≤ √6 + √(2 n_b log L)` are cited, not formalized (not in Mathlib) |
| AA.9.1 (contraction) | `CappedScore.lean`: `norm_adjoint_le`, `norm_adjoint_sub_le`, `norm_sum_smul_le` | `‖Q†a‖ ≤ ‖a‖` and `|‖Q†a‖ − ‖Q†b‖| ≤ ‖a−b‖` whenever `‖Qy‖ ≤ ‖y‖`, which holds for pairwise orthogonal columns of norm ≤ 1 (orthonormal or zero) | — |
| AA.9.2 (influence of one block) | `influence_le`, `score_sub_le_blocks` | `|σ(v) − σ(v')| ≤ r(‖v_b‖) + r(‖v'_b‖)`; for several blocks, `√(∑_{k∈D}(r‖v_k‖ + r‖v'_k‖)²)` | Corollaries (a)-(c) of AA.9.2 (the cap bound `2c`, the logarithmic bound with `M`, unboundedness for the identity) are not formalized; (a) is the special case `m = 1` of `cap_contamination` |
| AA.9.3 (saturation bound, first claim) | `cap_score_le` | `σ_{r_c}² ≤ n_b c²` | The amplitude-free formula under full saturation is not formalized |
| AA.9.5(c) (contamination) | `cap_contamination` | With at most `m` contaminated blocks the capped score moves by at most `2c√m` | — |
| AB.A.2 (coincident events are not identifiable) | `Superposition.lean`: `coincident_pair_mem_oneEventClass`, `coincident_pair_margin_zero` | If `G(θ)` is linear in the tensor, the pair signal with `θ₁ = θ₂` lies in the one-event class and its distance to that class is 0 | The testing consequence (Theorem AB.A.1) is not formalized |

Controls (rule 40): `Controls.lean` and `Controls2.lean`. Every result in the table has a
positive control (a concrete computation, or an instance with a nonzero `Q` showing the
hypotheses can hold together) and a negative control for each hypothesis listed in the file
headers. Exempt from negative controls, with the reason: `tailFunction_lipschitz` (each of its
hypotheses already has one in the lemma it feeds), the Lipschitz composition and maximum
lemmas, and the trace identities (no hypothesis to drop).

Not yet formalized:
- Probabilistic steps: the Laurent–Massart tail of E.2 (chi-square sums), Gaussian
  concentration for Lipschitz functions (Tsirelson–Ibragimov–Sudakov) and the expected
  maximum of sub-Gaussian variables used in AB.C.1 / AA.9.5(b). Mathlib v4.35 has Gaussian
  laws and sub-Gaussian MGFs but not these inequalities.
- Other results: AB.A.3(a) station decomposition of the timing-tolerant margin; AB.C.3
  (bounded influence implies saturation); AB.C.2 (conditional tail; rotation-invariant laws
  on the sphere); AB.A.1 (Le Cam bound; total variation between Gaussians).
