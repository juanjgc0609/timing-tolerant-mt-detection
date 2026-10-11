# Independent audit 2 of the doctoral proposal

**Object:** `/mnt/user-data/outputs/docs/Doctoral_Proposal.md` (hereafter P; 8011 lines, about 96.000 words).
**Evidence consulted:** `/home/claude/work/` (PREREGISTRO_1j_1k.txt with its timestamps and those of the files, task2a/2b/2c_report.txt, analyze_2a/2b/2c_out.txt, rft_report.txt, station_corrections_report.txt, qm_report.txt, gen_realistic_report.txt, certified_fw_report.txt, log2b_attempt1.txt, task2b_selfcheck.txt, task2b_parts/) and `/mnt/user-data/outputs/docs/artifact_manifest_and_error_log.md`.
**Independence:** `edits/independent_audit_1.md` and `audit_status.md` were not read.
**Reading:** the Abstract, Sections 1 to 4, 7.5, 13 to 20 and 23, and Appendices Y.11 to Y.20, Z and AA were read in full; the rest (5, 6, 8 to 12, 21, 22, G to X, Y.1 to Y.10) was reviewed through headings and passages.
**Fixed rubric:** Novelty and positioning 20 %; Theoretical rigour 15 %; Experimental design 20 %; Empirical evidence 15 %; Internal coherence and writing 10 %; Feasibility 10 %; Reproducibility 10 %. Total = weighted sum x 10.

---

## 1. Verdict

**Score: 51 / 100. Simulated decision for a doctoral proposal: major revision before approval.**

P is an honest and exceptionally self-critical document. It records negative results, withdraws claims and declares deviations. Even so, it does not yet meet the standard of a defensible Q1 proposal, for three structural reasons:

1. **All the evidence supporting the single hypothesis G1 is synthetic, from a single generator and outside the declared regime.** The band is 1.4 to 2.8 Hz versus 0.1 to 2 Hz, the domain is 30 x 30 km versus 120 x 120 x 40 km, and the medium is homogeneous or a 1D gradient.
2. **The measured advantage of MT-FW rests on comparators that operate at a structural disadvantage.** The associators are starved of picks because an STA/LTA recovers only 10 to 17 % of the arrivals in white noise. QuakeMigrate uses STA/LTA onsets, while MT-FW uses a matched filter close to the true wavelet. BOMP and mechanism-corrected stacking were read at exact times on a coarse grid, without the timing tolerance that MT-FW has. In all cases the operating point was chosen by oracle on the evaluation windows.
3. **The convex framework that gives the thesis its name (group BLASSO) performs below the greedy detector that is deployed.** The theory of Appendix AA covers the method that loses, and the method that wins (MT-FW) has only partial guarantees (AA.4(e)-(h)).

Added to this is a traceability problem detected in this audit: **the preregistration stamp of Task 2B (14:30 UTC) is later than the files that contain its results** (section 6, I-1). The document is also about three to four times longer than is reasonable for a proposal, and repeats the same status paragraph in at least seven places.

The dominant virtue is the methodological discipline on synthetic data: rules written before running, negative controls that break the geometry, an ablation of the operator, a hierarchical bootstrap over 30 networks, Holm, an independent verifier and pinning of the BLAS kernel. The dominant sin is that this discipline was applied to a designed scenario, without real data, and that the title, Section 11, Section 13 and the budget remain organised around an earlier thesis (constraint generation and differentiable representation), not around G1.

---

## 2. Evaluation by criteria

| Criterion | Weight | Score (0-10) | Weighted | Justification |
|:--|:-:|:-:|:-:|:--|
| Novelty and positioning | 20 % | 5.0 | 1.00 | The positioning is exemplary: it acknowledges as direct antecedents [129] (BOMP with moment-tensor blocks), GRiD MT [132], [123], mechanism-corrected stacking [121], [122], [133], subspace theory [118], [137], EVT thresholds [125] and the conformal outlier test [139] (L439-470, L566-603, L2017-2046). The residual novelty is modest and the author describes it as such. It comprises the off-the-grid formulation with amplitudes in $\mathbb R^6$, which is not implemented; the identity certificate = GLRT, which P calls a "modest observation" (L576-577); a certified search that does not prune (App. Z); and the conformal calibration, which is a transfer. The method that wins is, in its algorithmic core, "a grid relative" of [129] (L2022-2023) with peak search and cap. SOTA comparators declared mandatory are missing: MALMI [110], the stacking of [112], REAL, GENIE, SIGVISA and the NETVISA-type azimuth baseline (L2682-2730). WEAK BASELINE in the part on close antecedents (see D2). |
| Theoretical rigour | 15 % | 6.0 | 0.90 | The proofs of App. AA that were checked are correct: existence by Banach-Alaoglu (AA.1), GLRT identity (AA.2), $O(1/k)$ rate of Frank-Wolfe with exact oracle (AA.4(b)), a posteriori bound (AA.4(d)), validity of the cell bound and exact B&B (AA.3, Cor. AA.1), finite-sample validity of the conformal threshold (AA.5(c)) and exact support recovery on the grid (AA.6(a)), whose constants were checked. The Proposition / Conjecture / Approximation labelling is rigorous and AA.8 lists what is not covered. Four limitations weigh. (i) The results are standard adaptations of [113], [124], [139] and [47]. (ii) Recovery on the continuum is Conjecture AA.1. (iii) The deployed algorithm has no rate and no monotone descent when the cap is active. (iv) The theory does not cover model error (AA.8), which is precisely what the evidence measures. The observability certificate is leading-order for a reduced model. |
| Experimental design | 20 % | 5.0 | 1.00 | In favour: decision rules written before execution, MT-FW-NODE control, ablation of the operator (2B), close comparators (2C), hierarchical bootstrap with 30 networks, Holm and a deployable rule 2A-e. Against, there are five problems. (a) Operating point by oracle on the evaluation for all methods (L2294-2296), and the only deployable rule (2A-e) fails in its strict reading. (b) The comparators are at a structural disadvantage: STA/LTA that recovers 9.6 to 16.8 % of the arrivals (L6455-6457), BOMP and MTSTACK without timing tolerance (L7100-7104), QuakeMigrate with STA/LTA versus a near-oracle matched filter (L5967-5968). (c) A single generator, white noise, double-couple sources inside the detector subspace and two seeds. (d) There is no ablation of the peak search within MT-FW, although it is the element to which the advantage is attributed (L7087-7089). (e) Sections 11 and 13 (Tier 1 to 3, constraint-attribution experiment) do not contain the design of G1 (section 6, I-5 and I-6). |
| Empirical evidence | 15 % | 4.0 | 0.60 | The figures are internally solid: 16 of 16 verified figures match the files (section 7). But the evidence is entirely synthetic at the core, as the Abstract itself admits (L53). The only real data (T.1) are unfavourable to the polarization channel, which passes the gate at 1 of 5 stations. The advantages over close antecedents in the 4-station cells with strong clutter do not survive Holm (Table Y.41). The 2B ablation passes at the minimum with weak clutter (59 % of bootstrap replicates). Against GaMMA with optimistic picks, MT-FW loses at 0.1 in 4 of 6 weak cells. There are claims of real-world scope ("national network", Section 2) without support: *Undemonstrated speculative assertion*. |
| Internal coherence and writing | 10 % | 4.0 | 0.40 | The technical English is correct and the terminology precise. However, the length (about 96.000 words) is disproportionate for a proposal, and the same summary of 2A to 2C is repeated in L53-64, L259-311, L838, L1123-1185, L2439-2473, L2929-2953, L3102-3139 and in Y.18a and Y.20a. The title and the presentation announce a thesis different from G1. Nine text-text or text-evidence inconsistencies were detected (section 6), including a dangling reference (L145), a self-imposed taxonomy that the text itself violates (L774-776 versus L817 and L825) and the contradiction between "every claim is bound to a declared operating regime" (L34) and "the evidence ... does not lie inside this regime" (L241). |
| Feasibility | 10 % | 5.0 | 0.50 | The real test of G1 (decimation of the RSNC) is cheap in compute and well specified (L2788-2814). But semester 1 (L2895) accumulates too much: real test of G1, orientation audit, calibration of C-IV (tasks 1c to 1n), REAL, GENIE, azimuth, second generator, compute and domains; the gate requires the real result of G1 at month 6 (L2909). There are 12 specific objectives, 20 hypotheses and 7 publication units. The budget of 1670 to 2870 GPU-h (Table 5) funds conditional lines (3D elastodynamics, learned representation) and not the experiment that decides the thesis. There is an unresolved risk: MT-FW depends on the coherence of sign and polarization across stations, and T.1 shows that P polarization is reliable at 1 of 5 stations. The "per-station ladder" mechanism (L1129-1130) is untested. |
| Reproducibility | 10 % | 7.0 | 0.70 | This is a real strength: code per task, seeds, bit-for-bit control of the regeneration of windows, pinning of the OpenBLAS kernel (E31), an independent verifier with its own scripts (`verif_2a/`), an artifact manifest and a register of 33 errors. Deductions: the preregistration is a local, mutable text file (mtime 20:47 UTC on 26-09, later than all the results), without hash or external deposit (OSF or AsPredicted), and its 2B stamp is incompatible with the files (I-1). There is no public repository or DOI until month 42-48 (L2902). Bit-for-bit reproduction depends on the BLAS kernel of the machine. |
| **Total** | **100 %** | | **5.10** | **51 / 100** |

**Score a reasonable committee would give today: around 60 / 100, conditional approval.** A doctoral committee rewards the volume and honesty of the preliminary work, the correct theory and the registration discipline, and penalises the absence of real data in a proposal less than a Q1 reviewer does. Even so, it would require four things: reduce the document to 40 to 60 pages centred on G1; align title, objectives, Section 13 and budget with G1; present a minimal real pilot before candidacy; and explain the 2B timestamp.

---

## 3. Ten main weaknesses, by severity, with feasible corrections

**D1 (critical). Exclusively synthetic evidence, from one generator and outside the declared regime.**
G1 rests on Y.10 to Y.20. All those tasks use Ricker wavelets of 1.4 to 2.8 Hz, sampling at 16 Hz, white noise, a 30 x 30 km domain, a homogeneous or 1D-gradient medium and double-couple sources inside the detector subspace (L241-248, L6540-6545). In that noise PhaseNet obtains an F1 of 0.025 on P (L5792), which indicates a scenario very far from real data. P acknowledges this, but the Abstract and Section 20 present the advantage as the central result.
*Correction:* before candidacy, run a minimal real pilot: one dense subnetwork of the RSNC, one month, three decimated subnetworks of 4 and 6 stations, MT-FW+SC versus QuakeMigrate and GaMMA with PhaseNet picks, conformal threshold on noise windows and a recall versus false events per day endpoint above $M_c$. On synthetic data, add real noise (Task 1p) and non-double-couple sources. The cost is weeks, not months: P itself says it does not require a waveform-grade velocity model (L2813).

**D2 (critical). Comparators at a structural disadvantage; the advantage does not isolate a merit of the method.**
- (a) GaMMA and PyOcto are "starved": only 20 % of the events receive 3 picks (L6456-6457), because the STA/LTA was calibrated in a domain of isolated pulses in white noise.
- (b) QuakeMigrate receives STA/LTA onsets, while MT-FW correlates with a 2 Hz Ricker (L5967-5968).
- (c) BOMP and MTSTACK read the matched filter at the exact predicted time on a 3 km and 0.1 s grid. P itself shows that with sources on the nodes BOMP outperforms MT-FW (0.97 versus 0.76, L7054-7056) and that "the element that matters on these data is the peak search" (L7087).

The comparison with antecedents therefore measures the timing tolerance and not the method. Label: WEAK BASELINE.
*Correction:* (i) rerun of 2C with BOMP+SC, MTSTACK+SC and GL-FW+SC equipped with the same peak search of ±0.25 s, or with per-station shifts as nuisance parameters; (ii) ablation of MT-FW without peak search and MT-FW without cap; (iii) QuakeMigrate fed with the same matched filter as onset function, or MALMI with learned probabilities; (iv) associators with a picker that reaches at least 40 % of arrivals, for example an STA/LTA on the output of the matched filter. Everything runs in hours on the same 1200 windows.

**D3 (high). Operating points chosen by oracle on the evaluation; the only deployable rule fails.**
In 2A, 2B and 2C each method chooses its threshold, and the associators also their configuration, on the evaluation windows (L2294-2296). Rule 2A-e, with thresholds frozen on another seed, fails in its strict reading: the realized rate is 0.22 to 0.23 versus a limit of 0.20 (L6429-6438). The primary endpoint of P is "recall at a fixed realized false-event rate", but the rate is only "fixed" ex post.
*Correction:* redefine the primary synthetic endpoint with a threshold frozen by conformal calibration on noise windows of the same regime, which P already validates in Y.17, and report the realized rate with its interval. Rewrite the strict or per-cell reading of 2A-e before the next test.

**D4 (high). The convex framework that titles the thesis performs below the deployed greedy method; the theory covers the method that loses.**
GL-FW+SC, solved with a relative gap of up to $1.4\times10^{-8}$, stays between 0.06 and 0.61 below MT-FW in 23 of 24 cases (L1976-1980, Table Y.39). MT-FW has only AA.4(e)-(h): no rate and no guaranteed descent with an active cap. That sliding Frank-Wolfe will close the gap is an *Undemonstrated speculative assertion* (L1989-1995).
*Correction:* (i) reformulate the methodological contribution as "greedy detector with peak search robust to grid and timing errors, with its partial theory", and keep the BLASSO as framework; (ii) implement in 4 to 6 weeks a prototype of 1D or 2D sliding FW (in epicentre and origin time) on the windows of 2A, with debiasing, as a feasibility experiment prior to T2; (iii) prove, for MT-FW with cap, at least the condition $\alpha^*$ of AA.4(g), which already guarantees descent, and implement it.

**D5 (high). Traceability of the preregistration.**
Section 2B of the preregistration declares "14:30 UTC" and "before computing it" (L6564-6566, PREREGISTRO L325). However, `task2b_selfcheck.txt` is from 14:17 UTC, `log2b_attempt1.txt` (with the per-network recall of MT-FW-FREE+SC) from 14:20, `task2b_parts/` from 14:19 to 14:25 and `analyze_2b_out.txt` from 14:26:43 UTC. Either the stamp is wrong or the registration came later. In both cases the "preregistered" status of 2B is not verifiable. Moreover, the registration file is mutable and unique, and 2C was registered after seeing 2A and 2B on the same windows (declared in L6764-6766).
*Correction:* reclassify 2B as "rule written alongside the execution, not verifiable as prior" until evidence is provided (for example, an edit history); from now on, deposit each registration on OSF with a SHA-256 hash of the rules file and of the frozen code before execution.

**D6 (high). Structural misalignment between the declared thesis (G1) and the apparatus of the document.**
- The title and subtitle (L2-3, L10-13) announce "Information Sufficiency, Feasibility", "Constraint Generation" and "Differentiable Wavefield Representation"; the last is later declared "not a contribution of the thesis" (L339, L1239, L3031).
- Section 11 sets the "central comparison" against ABL and SIGVISA (L2308-2312), but ABL was retired as the reference for G1 (L1133-1135).
- In Table 6, G1 appears in Tier 2 of "constraint attribution" (L2565), whose design (L2597-2600) does not include MT-FW, QuakeMigrate or BOMP.
- The Tier 1 pilot "does not test G1" (L2552), but L790 calls it "the confirmatory programme on real data".
- The budget (L2347-2373) funds the conditional lines.

*Correction:* rewrite title, Abstract, 4.5, 11, 13, 12 and 19 around G1. Families A to F, the representation and optimal transport move to a "secondary programme" appendix. Target: 40 to 60 pages.

**D7 (medium-high). Fragile statistics in the secondary results that support the attribution.**
- 2B passes with the minimum of 2 of 3 cells with weak clutter and survives in 59 to 62 % of replicates (L6681-6687); with five networks per cell.
- In the 4-station cell with strong clutter, MT-FW does not outperform BOMP or MTSTACK after Holm, and on seed 7331 not QM-3C+SC either (p Holm 0.044 and 0.098; Table Y.40, analyze_2c_out L177-179).
- Only 2 seeds. The replication unit (network) is respected only in the comparison with QuakeMigrate (L2630-2631).
- The history 1j to 2A amounts to about twelve tasks in four days, with reformulations after negative results: the garden of forking paths is not corrected at programme level.

*Correction:* repeat 2B and 2A-c on 30 networks with hierarchical bootstrap; declare the complete family of tests of Section Y and a global correction (BH over 2A to 2C); add a third block of seeds never used in any previous task.

**D8 (medium). The physical channel that MT-FW exploits is in doubt on real data.**
MT-FW requires coherence of sign and relative amplitude across stations (L1914-1916). T.1 shows that P polarization passes the gate at 1 of 5 stations of the RSNC, and at the only one examined in depth it is at chance level (L3170-3172). 2B further indicates that with weak clutter much of the advantage over QuakeMigrate does not need the operator (Table Y.36). That MT-FW retains an advantage with a per-station ladder (L1129-1130) is an *Undemonstrated speculative assertion*.
*Correction:* a registered synthetic experiment with 1 of 4 stations reliable in polarization and the rest with random or free direction, to measure whether MT-FW with the ladder maintains G1; on real data, report MT-FW-FREE+SC as a mandatory arm of the decimation.

**D9 (medium). Oversizing of the scope and of the text.**
There are 20 hypotheses, 12 objectives, 7 publication units, 5 constraint families, learned representation, optimal transport, SIGVISA and p4est, in a document of about 96.000 words with massive redundancy. A committee cannot assess what is being defended without reconstructing it. The risk of dispersion is real: the first semester is already overloaded (L2895, L2909).
*Correction:* reduce to G1 plus two secondary hypotheses (sufficiency ladder and observability certificate) and move everything else to "conditional future work"; condense Appendices G to Y into a 10-page appendix with a traceability table (task, rule, result, file) and the rest as supplementary material.

**D10 (medium-low). Inaccuracies in the Abstract and in the summaries.**
The Abstract states "the better of block orthogonal matching pursuit and mechanism-corrected stacking by 0.08 to 0.52, in all twelve cells" (L55-57) without saying that this holds only at cap 0.25. At cap 0.1 the difference is +0.026 with an interval that includes zero (Table Y.38). It also omits that the lower bound of 0.08 corresponds to a cell whose intervals include zero and that does not survive Holm. It also omits the Holm exception against QuakeMigrate. The figure "17 to 22 times for the empirical ones" (L2804, L2840, L6160) does not match Table Y.28 (18.8 to 22.2).
*Correction:* specify the cap in each range; include "not significant after Holm in the strong-clutter 4-station cells"; correct to "19 to 22".

---

## 4. Real strengths

1. **Honest and up-to-date positioning.** It identifies and runs the closest antecedents ([129], [122]), acknowledges that its statistic is equivalent, up to normalization, to the variance reduction of GRiD MT, and formulates the novelty as "absence in targeted searches, not priority" (L581, L2037). This is rare in doctoral proposals and protects against OVERCLAIMING.
2. **Culture of registration and of negatives.** There are 33 documented errors with symptom, cause and detection; four positives withdrawn after audit (E24 to E26 and the ties effect); controls that break the geometry (MT-FW-NODE); an ablation of the operator; and rules with an interpretation written in advance ("if in neither, the claim ... is withdrawn", L6312-6314).
3. **False-alarm calibration well resolved.** Y.17 refutes with data the EEC prediction of $\chi^2_6$ (it overestimates $P(M>u)$ by between 20 and 60 times), shows that strong clutter breaks every threshold calibrated on Gaussian noise and adopts the conformal threshold with a finite-sample guarantee (AA.5(c), correct proof). It is a clean and transferable methodological contribution.
4. **Exemplary computational verification.** There is a certified oracle bit for bit against the dense loop, exact regeneration of 1200 windows, detection and pinning of the BLAS kernel, an independent verifier with FROC re-derivation and positive controls for GaMMA and PyOcto (0.88 and 0.69).
5. **Correct and well-delimited theory.** AA.6(a), recovery on the grid with explicit constants, and AA.6(c)-(d), curvature of the certificate equal to the efficient Fisher information with the mechanism as nuisance parameter, legitimately connect the observability certificate with the non-degenerate source condition.
6. **Statistics of 2C appropriate to the replication unit.** Hierarchical network-window bootstrap over 30 networks, bootstrap p-values with Holm over the family of 36 comparisons and declaration of 9 deviations.

---

## 5. Undemonstrated speculative assertions (selection)

- L1989-1995: that sliding Frank-Wolfe with per-station shifts and debiasing will make the convex method competitive. *Undemonstrated speculative assertion.*
- L2812-2814: that the real test of G1 "needs no waveform-grade velocity model" for MT-FW, whose advantage on synthetic data depends on the radiation and polarization geometry. *Undemonstrated speculative assertion.*
- L1129-1130: applying the physics "at the level of the fidelity ladder each station supports" without any experiment that validates it. *Undemonstrated speculative assertion.*
- L3186-3187: "a route to forming events where current pipelines cannot". *Undemonstrated speculative assertion.*
- L436-437: the regime crossover with SIGVISA as a contribution, without a plan for access to its code or an estimated cost. *Undemonstrated speculative assertion.*
- L378-379 (Section 2): determining epicentres "where current pipelines require three or more detecting stations"; in X.4 every method emits false events in 60 to 100 % of the realizations at 3 stations. *Undemonstrated speculative assertion.*

---

## 6. Text-evidence and text-text inconsistencies (with line numbers of P)

| ID | Lines of P | Inconsistency | Evidence |
|:--|:--|:--|:--|
| I-1 | L795-796; L6564-6566 | "each registered before it was computed"; "written ... at 14:30 UTC ... before any quantity of Task 2B was computed". | `task2b_selfcheck.txt` 14:17:01 UTC; `log2b_attempt1.txt` 14:20:05 UTC, with per-network recall of MT-FW-FREE+SC; `task2b_parts/*` 14:19 to 14:25 UTC; `analyze_2b_out.txt` 14:26:43 UTC; `task2b_report.txt` 14:30:47 UTC. All outputs are earlier than the registration stamp. |
| I-2 | L55-57 (Abstract) | Range against BOMP and MTSTACK "0.08 to 0.52 ... in all twelve cells" without cap and without the statistical exception. | Table Y.38 (L6906-6919): only at cap 0.25; at cap 0.1, +0.026 [-0.074, +0.156] in the strong 4-station cell (seed 7331); at 0.25, +0.077 with intervals that include zero; not significant after Holm (Table Y.41; analyze_2c_out L169, L178-179). |
| I-3 | L34 versus L241-248 | "Every claim is bound to a declared operating regime (Section 1.5)" versus "The evidence gathered so far does not lie inside this regime". | Direct internal contradiction. |
| I-4 | L774-776 versus L817, L825, L830, L841-842 | "Five labels are used throughout and no others". | The master table uses "preliminary support", "open", "preliminary null" and "conditional on distance spread". |
| I-5 | L789-791 versus L2550-2553 and L2788-2814 | L790: "The confirmatory programme on real data is the Tier-1 pilot of Section 13.0". | Section 13.0 defines a pilot of synthetic realizations (30 stations, $f_{\max}=0.8$ Hz, "one medium") that "does not test G1"; the real test of G1 is the decimation of the RSNC (Section 16, E1). |
| I-6 | L2308-2312 (Section 11); L2565 (Table 6); L2597-2600 (13.4) | The "central comparison" is against ABL and SIGVISA; G1 appears in Tier 2 of constraint attribution, whose conditions do not include the comparators of G1. | L1133-1135: ABL "is retired as the reference for G1"; L1131-1133: the comparators of G1 are published associators and migration-based locators. |
| I-7 | L2-3, L10-13 (title) versus L339, L1239, L3031 | The title highlights "Differentiable Wavefield Representation" and "Constraint Generation". | The representation "is an option, not a contribution of the thesis"; the single hypothesis is G1 (MT detector). |
| I-8 | L2804, L2840, L6160 | "17 to 22 times for the empirical ones". | Table Y.28 (L6147-6150) and `rft_report.txt` L447: GEV 21.86, Gumbel 18.80, conformal 22.18 (Bonferroni 19.40). The 17.22 corresponds to EEC-time, which is not empirical. Correct: "19 to 22". |
| I-9 | L145 | "Section 11.1 separates the two". | Section 11 has no subsection 11.1 (L2284-2315). Dangling reference. |
| I-10 | L2622 versus L2650 | "For $E_{\mathrm{det}}$ the replication unit is the network" versus "The independent replication unit is the complete realization". | Ambiguity of unit within the same section. Minor. |
| I-11 | L1129-1133 (G1: $\ge0.05$ at 0.1 and at 0.25 versus the "best-tuned published associator") versus L838 and L1135 ("supported") | With optimistic picks, GaMMA-DY+SC outperforms MT-FW+SC at 0.1 in 4 of 6 weak cells (by up to 0.11), and in single-source windows in 5 of 6. G1, as formulated, does not hold in that arm, although P declares it as a qualification. | analyze_2a_out Table 2 (thin), column "- GaMMA-DY+SC" at 0.1: -0.02, -0.01, -0.08, +0.04, +0.04, -0.11. |

No numerical discrepancies were found between P and the result files in the verified figures (section 7), except I-8, which is already in `rft_report.txt` itself.

---

## 7. Verification of figures and references

**Figures (16 verified; all match, except for the label of I-8):**

| # | Figure in P | Lines | File | Result |
|:-:|:--|:--|:--|:--|
| 1 | One-sided 8 st. (2024), weak: MT-FW+SC 0.55/0.65; QM-3C+SC 0.35/0.38; GaMMA-DY+SC 0.63/0.63 | L6322 | analyze_2a_out L13 | matches |
| 2 | 2A-a: +0.09 to +0.49 (+0.14 to +0.49 and +0.25 to +0.46 weak; +0.09 to +0.37 and +0.17 to +0.41 strong) | L6389-6391 | analyze_2a_out Tables 3 | matches |
| 3 | 2A-e: realized rate 0.22 (weak) and 0.23 (strong) versus limit 0.20; "up to 2.3 times" | L6429-6431 | analyze_2a_out L56, L132 | matches |
| 4 | STA/LTA recovers 9.6 to 16.8 % of the arrivals | L6455 | analyze_2a_out Table 5 (0.096 to 0.168) | matches |
| 5 | SC residual: 0.12 to 0.21 s (weak), 0.20 to 0.28 s (strong), statics 0.26 to 0.34 s | L6416-6417 | analyze_2a_out Table 5 | matches |
| 6 | Optimistic clutter picks: 27.25 to 32.60 per window | L6402 | analyze_2a_out L148-153 | matches |
| 7 | 2B: rule in 59 and 62 % of replicates (weak), 98 and 97 % (strong) | L6683-6684 | analyze_2b_out (thin) | matches |
| 8 | 2C-c: 4 st. (7331) weak, MT-FW+SC 0.466, QM-3C+SC 0.160 | L6950 | analyze_2c_out L107 | matches |
| 9 | Holm: strong, 4 st. (7331), +0.099 [+0.029, +0.181], p Holm 0.044 | L7006 | analyze_2c_out L177 | matches |
| 10 | GL-FW versus MT-FW: -0.06 to -0.61 in 23 of 24; the only "equal" -0.013 | L6991-6994 | analyze_2c_out L92 | matches |
| 11 | Conformal on Gaussian noise: 0.007/0.037/0.095/0.239 | L1960, L6130 | rft_report L64 | matches |
| 12 | Strong clutter, calibrated on Gaussian: 0.60/0.78/0.90/0.96 | L1965, L6159 | rft_report L75 | matches |
| 13 | K=20: combined rms 0.088 s; with catalogue errors 0.173 s | L6015-6016 | station_corrections_report L35-36 | matches |
| 14 | MT-FW - QM-3C: +0.28 [+0.10, +0.43] at 0.25 | L5943 | qm_report L180 | matches |
| 15 | Picked fraction: degrade_y 0.420, STA/LTA 0.143, PhaseNet 0.041 | L5818-5820 | gen_realistic_report L133-135 | matches |
| 16 | Certified oracle: 215.982 of 225.444 evaluations; 13.1 s versus 55.9 s | L7240-7242, L7260 | certified_fw_report L14, L102 | matches |

Also checked by hand were Table 3 of solver cost (h, $N_{\mathrm{dof}}$, steps and PFLOP are consistent with $h=V_{S,\min}/(pf_{\max})$ and work proportional to $f^4$), the sample size (25), where $n=100$ with $\sigma_d=0.178$ implies $\Delta\approx0.05$, and the constants of the proof of AA.6(a).

**References (6 checked; 3 by web search):**
- [109] Winder, Bacon, Smith, Hudson, White, "QuakeMigrate ...", *Seismica* 5(1), 2026, doi 10.26443/seismica.v5i1.1854. Verified on the web: published on 7 May 2026, vol. 5, no. 1, matching authors.
- [16] Puente Huerta, Sippl, Munchmeyer, McBrearty, "Benchmarking seismic phase associators: Insights from synthetic scenarios", *Seismica*. Verified on the web (article 1559, arXiv 2501.03621).
- [125] Hirano, Kawakata, Doi, "A matched-filter technique with an objective threshold", *Sci. Rep.* 12, 22090 (2022), doi 10.1038/s41598-022-25839-2. Verified on the web.
- [129] Vera Rodriguez, Sacchi, Gu, *GJI* 188(3):1188-1202, 2012; [47] Duval and Peyre, *FoCM* 15:1315-1355, 2015; [116] Denoyelle, Duval, Peyre, Soubies, *Inverse Problems* 36:014001, 2020. They match the known bibliographic record. The use P makes of them is correct; in particular, it acknowledges that Theorem 3 of [116] is scalar and 1D (L7687-7692).

Web sources: [Seismica, QuakeMigrate](https://seismica.library.mcgill.ca/article/view/1854); [Seismica, Benchmarking seismic phase associators](https://seismica.library.mcgill.ca/article/view/1559); [arXiv 2501.03621](https://arxiv.org/abs/2501.03621); [Scientific Reports, matched filter with objective threshold](https://www.nature.com/articles/s41598-022-25839-2).

---

## 8. What would raise the score by 5 points (from 51 to about 56)

The shortest and most feasible route in 4 to 6 weeks, without real data, adds +5 points (+6 if it also improves novelty):

1. **Comparators with equal timing tolerance and ablation of the peak search** (D2): rerun of 2C with BOMP, MTSTACK and GL-FW with ±0.25 s per station and phase, plus MT-FW without peak search and without cap, on the same 1200 windows and the 30 networks. Experimental design goes from 5 to 6 (+2.0 points). If MT-FW retains the advantage, novelty also rises from 5 to 5.5 (+1.0).
2. **Deployable endpoint** (D3): the conformal threshold frozen on noise windows of the same regime as the primary reading, with the realized rate reported. It is included in the +2.0 of design.
3. **Restructuring of the document** (D6, D9, D10 and I-2 to I-11): title, Abstract, 4.5, 11, 12, 13 and 19 aligned with G1, a maximum of 60 pages and the rest moved to supplementary material. Coherence goes from 4 to 6 (+2.0).
4. **External preregistration with hash and clarification of I-1** (D5): reproducibility goes from 7 to 8 (+1.0).

With a minimal real pilot of decimation of the RSNC (D1), evidence would go from 4 to 6 (+3.0) and feasibility from 5 to 6 (+1.0), and the score would reach about 60 to 61 on this rubric.
