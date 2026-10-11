# Independent audit 4 — Doctoral proposal «Timing-Tolerant Moment-Tensor Detection of Overlapping Earthquakes»

Date of the review: 2026-10-02. External reviewer without access to the production process or to previous audits.
Object: `Doctoral_Proposal.md` (4102 lines, read in full) and
`Supplementary_Evidence_Appendices_G_to_AA.md` (Y.18–Y.26, Z and AA read; the rest, through headings and summaries).
Evidence checked: registration and hash chain, outputs `analyze_2a…2f`, `analyze_2i_primary_out.txt`,
verifications 2E/2F/2I, `explore_contrib/REPORT*.txt`, error register (E1–E39).
Line numbers refer to the main document unless marked «Suppl.».

---

## 1. Verdict

**Score with the fixed rubric: 57,8 / 100. Simulated decision: major revision (approvable as a doctoral proposal
only with a reduction of scope and a reformulation of the title and of the central thesis).**

The document is, in its registration discipline, far above the standard of a doctoral proposal: the figures
cited match the result files in all the checks performed, the hash chain is
intact and the negative results are reported without embellishment. Its main defect is one of substance and not of
honesty: after seven preregistered tasks, what remains standing is (i) a synthetic advantage of a family
of detectors whose decisive ingredient is a peak-search heuristic of ±0,25 s that the antecedent [129]
equals as soon as it is granted it, (ii) a gain of 0,02–0,07 in recall from bounding the influence of each station,
a technique that the text itself declares standard, and (iii) a theory (Appendix AA) that covers a convex problem that is not
the detector that performs. All the evidence for the primary hypothesis is synthetic, from a single generator family,
outside the declared operating regime, with the operating point chosen on the evaluation windows, and the only
measurement on real data (reliable polarization at 1 of 5 stations of the RSNC) points against the physical assumption of the
operator that the detector uses. The framework of «five constraints generated from elastodynamics», which occupies the
bulk of the text, takes part in none of the confirmatory tests.

---

## 2. Table by criteria

| Criterion | Weight | Score (0–10) | Contribution | Justification |
|:--|:-:|:-:|:-:|:--|
| Novelty and positioning | 20 % | 5,5 | 1,10 | Exhaustive and truthful positioning against [129], [122], GRiD MT, NETVISA/SIGVISA. But the novelty that survives is avowedly empirical: per-station clipping and Huber fitting are standard (l. 77–78, 716–718, 2276–2280), the peak search is a transferable heuristic, G1-M was refuted. The bibliographic search protocol has not been executed (l. 3689–3692): every claim of absence of an antecedent is an *Undemonstrated speculative assertion* until it is executed. |
| Theoretical rigour | 15 % | 6,0 | 0,90 | The proofs of AA.1–AA.6 are correct and their hypotheses are declared, but they are known results (group-lasso optimality, Jaggi's rate, conformal validity) applied to the problem. The theory covers (16a), not the deployed statistic: AA.8 admits it. The «robust certificate» has only Lemma AA.3 (trivial influence bound) and AA.4(g) shows that with an active cap the residual does not even decrease. Recovery on the continuum is Conjecture AA.1. The Proposition of 6.1 has a remainder «stated, not bounded» (l. 1632–1634). |
| Experimental design | 20 % | 7,0 | 1,40 | Real strength: rules written before running, bit-for-bit controls, ablations (2B, 2D), geometry-breaking controls, independent adversarial verification (2E, 2F, 2I), hierarchical bootstrap by network, superiority / non-inferiority / equivalence rules. Weaknesses: operating point chosen on the evaluation in all tasks; a single generator; regime different from the declared one; STA/LTA picker that leaves GaMMA and PyOcto without picks; five networks per cell in 2A, 2B and 2D; equivalence margins (0,05) larger than the claimed effect (0,03). |
| Empirical evidence | 15 % | 4,5 | 0,675 | Exclusively synthetic for G1, G1-M and G1-R. G1, as formulated (≥0,05 at 0,1 **and** 0,25), does not hold against GaMMA with optimistic picks in 4 of 6 cells of weak clutter. Effect of G1-R: 0,043–0,074, which falls to 0,019–0,045 when origin time is required; absolute recall in dense windows 0,22–0,29 (0,07–0,10 with origin). There is no deployable operating point. The only real data are unfavourable. |
| Internal coherence and writing | 10 % | 5,5 | 0,55 | Remarkable text-file numerical consistency (no discrepancy in the checks). But: 49 000 words plus 80 000 of supplement; section 1.6 is ~3100 words in two paragraphs; the abstract has 299 words; the same inventory of results is repeated in the abstract, 1.6, 4.3, 7.5, 13.0c, 17, 20 and 23; the title and the thesis of 1.3 no longer describe what the evidence supports; sections 5–12 inherit a programme (PINN, optimal transport, Gaussian beams) that the text itself downgrades to a conditional option. |
| Feasibility | 10 % | 5,0 | 0,50 | The decisive test (G1' by decimation of the RSNC) is cheap and well prioritized (month 6). But the plan keeps twelve objectives, twenty-five hypotheses, two tracks, an elastodynamic solver, SGT, learned representation and seven publication units for a single doctoral student; E1 (months 1–6) accumulates decimation, orientation audit, REAL, GENIE, azimuth baseline and second generator. The budget is a «planning model» without confirmed allocation. |
| Reproducibility | 10 % | 6,5 | 0,65 | Verifiable hash chain, code hashes before running, script manifest, register of 39 errors, bit-for-bit regeneration, BLAS kernel pinned and documented. Against: nothing is public (no repository, no DOI, no third-party seal); the hashes live in the same container as the registration; the results depend on the OpenBLAS kernel (E31); the environment is not packaged; the registration of 2A–2C has no contemporaneous hash. |
| **Total** |  |  | **5,775** | **57,8 / 100** |

**Score a reasonable doctoral committee would give today: ≈ 70 / 100** (conditional approval). A committee rewards
honesty, the volume of prior work and statistical maturity more than this rubric does, and penalises
synthetic-only evidence less at the proposal stage; it would nonetheless require cutting the scope, changing the title
and fixing the success criterion on real data.

---

## 3. Checks performed

**Hash chain.** SHA-256 of successive prefixes of `PREREGISTRO_1j_1k.txt` (cuts at end of line) against
`PREREGISTRO_hashes.txt`: the 12 hashes are reproduced, with monotone prefixes ending at lines 431, 436, 447,
497, 501, 508, 542, 550, 557, 588, 596 and 604 (end of file, 56 925 bytes, current hash `c2745a94…`). The hash of the
2D registration ends at the end of the section «Task 2D» (l. 431) and that of 2I at the end of «Task 2I» (l. 588). Limit:
lines 1–431, which contain the registrations of 2A, 2B and 2C, have no hash earlier than that of 26-Sep 21:13, so
their content at the time of writing is not verifiable (the text declares this, l. 938–958). The hash file
is local and mutable; the chain proves consistency, not priority.

**Numbers against files (14 checked, 14 match).**

| Claim in the text | File | Result |
|:--|:--|:--|
| 2A, one-sided 8 st. (2024): MT-FW+SC 0,55/0,65 | `analyze_2a_out.txt` l. 13 | matches |
| 2A-a, 6 st. (7331): +0,49 [+0,25, +0,65] at 0,1 | ibid. l. 37 | matches |
| 2A-e, one-sided (7331): realized rate 0,22 with budget 0,1 | ibid. l. 56 | matches |
| 2D-a: range −0,155 to +0,091 | `analyze_2d_out.txt` l. 106, 110 | matches |
| 2E-a: +0,039 [−0,025, +0,122], p = 0,135 | `analyze_2e_out.txt` l. 38 | matches (0,1349) |
| 2F-a: +0,072 [+0,048, +0,092]; +0,074; +0,051; +0,043 | `analyze_2f_out.txt` l. 26–30 | matches |
| 2F-b: +0,060 [+0,046, +0,076]; NOCAP−BOMP −0,001 | ibid. rule 2F-b | matches |
| 2F-c: p against the lower margin 0,0170 and 0,0405 | ibid. l. 49–50 | matches |
| 2I-a: +0,058 [+0,039, +0,076]; +0,053 [+0,037, +0,070] | `analyze_2i_primary_out.txt` | matches |
| 2I-b: −0,008 (bound −0,0234); −0,002 (bound −0,0184) | ibid. | matches |
| 2I-c: +0,0099 [−0,0100, +0,0272]; +0,0073 [−0,0132, +0,0235] | ibid. | matches |
| 2I with origin ≤0,5 s: +0,031 and +0,019; weak standard 6262 −0,032 | `verif_2i/…` l. 96, 106 | matches |
| Verif. 2E: −0,014 with 10× iterations; ceilings 0,08–0,24 | `verif_2e/…` l. 110, 157 | matches |
| C8: HUBER-N − BOMP-PK +0,064 [+0,028, +0,105]; slack fraction 0,21 / 0,14 | `REPORT_C8.txt` l. 44–45, 193 | matches |

The secondary arms of 2I appear as «NOT RUN (parts missing)» in the primary output, consistent with the
text («still running»); the process `task2i_run.py` was still active during the review and was not touched.

**References (4 checked online plus a check from memory).** [144] Pham, GJI 238(2):783, DOI resolves to
`academic.oup.com/gji/article/238/2/783`: correct. [109] QuakeMigrate, Seismica, article 1854: exists with that
title. [16] Puente Huerta et al., Seismica, article 1559 (preprint arXiv:2501.03621): correct. [9] *Nature
Communications* 2026 version of HARPA: **not located** in the search; only the preprint arXiv:2307.07572 is on record.
The volume, article number and DOI must be verified. From memory, [47], [116], [124], [129], [139], [142], [143], [43]
and [85] (174 events, 150 + 24) are correct. [19] and [23] are 2026 preprints that I could not check.

---

## 4. Ten weaknesses, by severity, with concrete correction

1. **All the evidence for the primary hypothesis is synthetic, from a single generator and outside the declared regime.**
   The tests use 1,4–2,8 Hz, 30×30 km, 4–8 stations, a 1-D gradient medium, white noise and double-couple sources
   inside the detector subspace; the primary regime is 0,1–2 Hz and 120×120×40 km (l. 215–239). The «fresh
   seeds» are replicates of the same generator, not independent evidence.
   *Correction:* before defending the proposal, (a) run the second generator already promised (another synthesis code,
   e.g. 1-D Green's functions by reflectivity, real noise from the RSNC, non-double-couple sources) and (b) deliver the
   first result of G1' on the RSNC. In the meantime, rewrite table 1.5 so that the declared regime is the
   one that is tested, or declare the primary regime as an objective without evidence.

2. **G1 does not hold as formulated and there is no deployable operating point.** The rule requires ≥0,05 at 0,1
   **and** at 0,25; GaMMA with optimistic picks is ahead at 0,1 in 4 of 6 cells of weak clutter (l. 271–273). The
   operating point of all methods is chosen on the evaluation; the frozen threshold exceeds the budget 2×
   (2A-e) and the conformal one leaves 1,2–2,2 false events per window with events, 3–8 times QuakeMigrate (l. 304–309).
   *Correction:* register and run a task with a frozen threshold on windows **with** events from disjoint networks
   (calibration of the 2A-e type with ≥30 networks), reporting recall at the frozen threshold with the realized rate; and declare G1 in
   the abstract as «not met at 0,1 against an associator with good picks», not as a subordinate caveat.

3. **The advantage over the state of the art rests on a transferable heuristic, and the comparison with associators is
   weak.** 2D showed that [129] with the same peak search equals MT-FW; the «family» advantage over GaMMA and
   PyOcto is obtained with an STA/LTA that hits 10–17 % of the arrivals (l. 2992–2994), and QuakeMigrate was calibrated
   only in weak clutter and in two of three cells (Suppl. l. 2885–2887). MALMI, [112], REAL and GENIE have not been run.
   *Correction:* repeat 2A with a picker that does not leave the associators without picks (PhaseNet retrained or fine-tuned on
   the generator, or real picks in G1'), QuakeMigrate calibrated per regime and cell, and MALMI; present the conclusion
   as «pick-free detection with timing tolerance and deflation», not as a method contribution of its own.

4. **The term «certificate» and the title claim more than has been demonstrated.** The title promises «A Robust Subspace
   Certificate» (l. 2, 10), but the text itself states that the statistic with peak search and cap «is not that
   certificate» (l. 400, 1292) and AA.8 treats it as a deviation without theory. There is no optimization problem of which the
   Huber or the clipped statistic is the dual certificate, nor a false-alarm guarantee, nor a recovery guarantee.
   *Correction:* call the object «bounded-influence subspace statistic» until T2 formulates the problem
   (e.g. BLASSO with a per-block Huber loss) and proves its optimality condition; retitle the proposal around
   the preregistered benchmark and the test on a decimated network.

5. **Small effect size and misaligned statistical margins.** G1-R gains 0,043–0,074 (0,019–0,045 with origin
   ≤0,5 s; at 6262 the point criterion of 0,03 would fail) with absolute recall ≈0,24. The equivalence margin on
   isolated events (±0,05) is larger than the superiority threshold (0,03) and of the order of the gain: a cost of
   up to 0,05 on isolated events, which are the majority in operation, is not excluded; in fact with weak clutter at
   6262 Huber loses 0,032 with an interval that excludes zero, and in 2F the strong clutter gives −0,026 and −0,032. The registered hit
   rule ignores the origin time, and ≈43 % of the location hits have an error >0,5 s.
   *Correction:* set as the future primary rule the hit with origin ≤0,5 s; equivalence with margin ≤0,02 and per
   regime; add a net-utility criterion weighted by the realistic proportion of dense and isolated windows.

6. **The only real data contradict an assumption of the operator and the plan does not address it for MT-FW.** P polarization
   passes the gate at 1 of 5 stations of the RSNC (l. 377–379, 1836–1841). The radiation operator imposes P
   along the ray and transverse S at all stations, and 2B indicates that with strong clutter the advantage depends on that
   operator. The risk table treats this only for C-IV (l. 3132), not for the primary detector.
   *Correction:* add an explicit risk «radiation operator invalid at unreliable stations» with its
   response: a variant of MT-FW with a free block (MT-FW-FREE) or weighted by $\hat r_s$ at stations that fail,
   registered as an arm of G1', and a stopping criterion if less than a declared fraction of the network passes the orientation
   audit.

7. **Disproportionate scope and incoherence between the framework and what is tested.** The thesis of 1.3 («every piece of
   physics is a constraint») and sections 5–12 (five constraints, optimal transport, beams, learned residual,
   interior-point) take part in no preregistered task; 708 of 1670 GPU-h are conditional on F1 and 605
   fund lines that «cannot decide the thesis» (l. 2628–2634). There are 25 hypotheses and 12 objectives.
   *Correction:* split into a core (G1/G1', G1-R, benchmark, T2) and an annex of optional lines without budget or
   committed milestones; reduce to ≤6 hypotheses with a threshold; remove from the body sections 7.3, 10 and 12.2–12.3 or
   move them to an appendix.

8. **The theory is of the framework, not of the method, and the open mathematical problem is not formulated.** The BLASSO on the grid
   stays 0,06–0,61 below the greedy method, off the grid it does not win, and what performs (peak + deflation + bounding) only
   has AA.4(e)–(h). Three quarters of the measured slack remain «unexplained» without a registered hypothesis about its cause.
   *Correction:* formulate in T2 a provable and falsifiable statement about the deployed statistic (e.g. control of
   false alarm under $H_0$ of the Huber statistic with maxima over ±w; a bias bound for the deflation with cap) and a
   registered hypothesis about the slack (global placement of atoms at 5–10 km) with its experiment.

9. **Non-external preregistration and incomplete results presented as support.** The hashes are stored alongside the
   registration; «sent to the user» is not a third party (the text admits it, l. 963–965); 2A–2C lack a hash; the
   2I registration was confirmed after a selection among 11 candidates and with the run interrupted twice; the
   level and dispersion arms —those that would say whether «the dependence on a constant is removed»— are not reported, and in the
   pilot Huber was worse than the cap at dispersion 0,8 and the gain disappears at c = 8.
   *Correction:* seal the current registration today with RFC 3161 or OSF and a public signed commit; do not claim «removes the
   dependence on a tuned constant» (l. 411–412, 3245) until the arms are reported; say «replaces a calibrated
   constant by a rule, with a useful interval of levels not yet measured on evaluation seeds».

10. **Length and writing prevent evaluation.** 130 000 words in total; abstract of 299 words; 1.6 is a
    paragraph of more than 3000 words; the same results are listed at least seven times; the date of the document
    (September 2026) is earlier than results it contains (2 October).
    *Correction:* body ≤12 000 words; a single table «claim – evidence – status – limit» in 1.6; abstract
    ≤250 words; move Y.1–Y.17 to a cited technical report.

---

## 5. Strengths

- **Text-file fidelity.** Fourteen figures checked, fourteen matches, including interval bounds and
  p-values; rounding differences declared (e.g. +0,0494 printed as +0,05, Suppl. l. 3067).
- **Negatives reported and with consequence.** G1-M was registered with a refutation rule and declared refuted; 2D
  withdrew the distinction of MT-FW from [129]; four previous positives were withdrawn after audit; 39 design
  errors documented.
- **Adversarial verification that narrows the scope.** The verifications of 2E, 2F and 2I reproduced FROC and
  bootstrap with their own code and added limits (hit rule with origin, a ranking mechanism and not a
  detection mechanism, 9× cost) that the text incorporates.
- **Well-designed breaking controls.** Subspace control from another node (2A-b), identity operator (2B),
  ablation of peak and cap (2D), oracle ceilings that show that the benchmark is not saturated (Y.22).
- **Above-average statistics.** Replication unit = network, hierarchical bootstrap with the operating point
  re-chosen per replicate, Holm, TOST and non-inferiority with margins fixed beforehand; sample size of 2I planned
  with an effect reduced by the selection.
- **Closest antecedents confronted and run**, not only cited ([129], [122], real QuakeMigrate 1.2.2).
- **Cheap and early decisive test**: G1' depends neither on the solver, nor on the learned representation, nor on a waveform-grade
  velocity model, and it has a gate at month 6.

---

## 6. Inconsistencies (with line number)

| # | Lines | Inconsistency |
|:-:|:--|:--|
| 1 | 5 versus Suppl. 4934, and l. 345–371 | The date of the document is «September 2026», but it contains results of Task 2I, registered and run on 2 October 2026. |
| 2 | 2, 10, 13–14, 404 versus 400, 1292 | The title and G1-R call the clipped or Huber statistic a «certificate»; the text says that the statistic with peak and cap «is not that certificate» and AA.8 treats it as a deviation without theory. |
| 3 | 897 versus 3409 and 3518–3519 | The row «Abandoned as claims» lists six entries (it includes G1-M); l. 3409 says «five claims abandoned in Section 4.3» and l. 3518–3519 «two refuted … three further». |
| 4 | 1022 versus 1005, 1232 | «E1 to E4 record what was ruled out», but E5 exists in the master table and in Family E. |
| 5 | 3197 versus 2876–2877 | Table 8 assigns G1 to package E4 (constraint attribution); section 13.4 states that this experiment «does not test G1». |
| 6 | 3437–3438 versus 1011, 1296–1302 | «In three preregistered synthetic tasks» versus «four preregistered tasks (Y.18 to Y.21)» for the same evidence of G1 (depending on whether 2D is counted or not); the count changes from section to section. |
| 7 | 2787, 2908–2911, 3663 versus 2921–2922, 2959, 3651 | Heterogeneous multiplicity correction: Holm in tasks 2C–2F and in G1'; Benjamini-Hochberg registered as the procedure of the confirmatory family. |
| 8 | 411–412, 3245, 3467–3468 versus 370–371, 3138; Suppl. 4866–4870, 5214–5217 | It is claimed that 2I «removes the dependence … on a tuned constant» while the level and dispersion arms are not reported and the pilot shows dependence on the level (null gain at c = 8) and worse performance at dispersion 0,8. |
| 9 | 168–171, 3488–3489 versus 173–183, 2165–2172 | The one-sentence thesis holds constrained optimization to be what solves the problem; the following text acknowledges that what passed the tests is a greedy search with a peak heuristic and not the convex solution or the constraints C-I to C-V. |
| 10 | 3626, 1285–1286 versus 3628, 1330 | G1 requires ≥0,05 at 0,1 and 0,25; G1-R, presented as the central methodological contribution, requires only ≥0,03 at 0,25, with an equivalence margin of 0,05 larger than the effect threshold itself. |
| 11 | 2816–2818, 3643 versus 2828–2834 | The Tier-1 pilot on which the «frozen» margins depend is run at 30 stations and 0,8 Hz, a configuration where the text declares that no effect is expected; the thresholds of G1 are frozen «after the Tier-1 pilot» (l. 807) but the pilot «does not test G1». |
| 12 | 3111–3112 versus 2253–2255, 3132 | G1' «needs no waveform-grade velocity model», but the radiation operator uses P and S amplitudes and polarities at 1–3 Hz and the polarization fails at 4 of 5 stations; there is no registered arm that addresses this. |
| 13 | 55–79 | Abstract of 299 words versus the usual limit of 250; the «mathematical» and «structural» contributions of 20.0 do not appear in it despite l. 3230 («the same four named in the Abstract»). |
| 14 | App. F [9] | The *Nature Communications* 2026 version of HARPA could not be located; verify it or cite only the preprint. |

No numerical discrepancy was found between the text and the result files.

---

## 7. What would raise the score

**+5 points (≈63), achievable in weeks, without new field data:**
1. Retitle and rewrite the abstract, 1.3 and 23 around what the evidence supports: preregistered benchmark,
   timing-tolerant family, design result on bounded influence; remove «certificate» from the title.
2. Cut the body to ≤12 000 words with a single evidence table, and relegate the unused sections 5–12.
3. Seal the registration externally (RFC 3161 or OSF) and publish code, pinned environment (container with declared BLAS
   kernel) and registration in a repository with a DOI.
4. Report the secondary arms of 2I as they were registered and adjust the sentence about the constant.
5. Correct inconsistencies 1, 3, 4, 5, 6 and 10 and execute the bibliographic protocol of Appendix E with counts.

**+15 points (≈73), requires new experimental work:**
1. **First result of G1' on the decimated RSNC**, preregistered with an external seal, with a frozen threshold and false events per
   day on continuous data, against QuakeMigrate, MALMI and associators with PhaseNet; any sign will do.
2. **Second independent generator** (another synthesis code, real noise, sources outside the double couple) for G1 and
   G1-R, with a hit rule that requires origin ≤0,5 s and per-regime equivalence with margin ≤0,02.
3. **Deployable operating point**: frozen calibration on windows with events in disjoint networks, with a realized
   rate within a declared factor of the budget.
4. **A theoretical result on the deployed detector** (false-alarm control or selection guarantee of the
   bounded-influence statistic with peak search), so that the mathematical contribution ceases to be the
   theory of a method that is not used.
5. Formal cut of the programme to a three-year core with ≤6 hypotheses and four publication units.
