# Independent audit 3 (Q1 reviewer)

**Object.** Doctoral proposal `Doctoral_Proposal.md` (hereafter M, 3719 lines, ~42 000 words, read in full) and supplementary volume `Supplementary_Evidence_Appendices_G_to_AA.md` (hereafter S; Y.11 to Y.21, Z and AA read in full; G, H, T.1, W and X reviewed). Evidence checked in `/home/claude/work` (preregistration, hashes, analysis outputs of Tasks 2A to 2D, component reports) and error register `artifact_manifest_and_error_log.md` (hereafter R). No previous audit was read.

**Fixed rubric.** Novelty and positioning 20 %; theoretical rigour 15 %; experimental design 20 %; empirical evidence 15 %; internal coherence and writing 10 %; feasibility 10 %; reproducibility 10 %. Total = weighted sum × 10.

---

## 1. Verdict

**Rubric score: 51/100. Simulated decision: Major Revision.** A reasonable committee today: **62/100**, that is, approval conditional on a substantial rewrite and a reduction of scope.

**Critical verdict.** The prior work is unusually honest and well instrumented. The registration hashes verify, the reproductions are bit for bit and the negative results are reported. But that same work dismantled the methodological contribution that the title and hypothesis G1 still announce. The detector that performs (MT-FW) is, with the same peak search, indistinguishable from the block orthogonal search of Vera Rodriguez et al. [129]: rule 2D-a fails and all the hierarchical intervals contain zero. The certified convex problem stays between 0.11 and 0.30 below the detector in 11 of 12 cases. Moreover, no operating point is deployable, all the evidence is synthetic and lies outside the declared regime, and the only distinctive contribution that remains, the off-the-grid group BLASSO with per-station time shifts, is not implemented. The proposal is today a benchmark programme with a mathematical promise, not a method with evidence. The main sin is one of positioning: the title ("with a Convex Certificate") and the formulation of G1 ("the certificate of a group BLASSO") describe an object that, when measured, is not the one that works.

---

## 2. Evaluation by criteria

| Criterion | Weight | Score (0-10) | Weighted | Summary justification |
|:--|:-:|:-:|:-:|:--|
| Novelty and positioning | 20 % | 4.0 | 0.80 | The block MT dictionary with deflation is [129]. The statistic at zero is the variance reduction of GRiD MT [132], [123]. The conformal threshold is a transfer of [139], the Frank-Wolfe rate is [124] and the BLASSO/SFW machinery is [113]-[117]. MT-FW does not outperform BOMP-PK (2D-a fails). What remains distinctive (per-station shifts as off-the-grid nuisance parameters) is not implemented, and moreover has uncited antecedents: the MT inversion with per-station time shifts of the CAP method (Zhu and Helmberger, 1996) and its differentiable version (Phạm, GJI 238(2):783, 2024). The preregistered benchmark against GaMMA, PyOcto and QuakeMigrate is indeed a contribution, but a synthetic one. |
| Theoretical rigour | 15 % | 5.5 | 0.83 | Appendix AA is correct as far as it was verified (constants of AA.4(b), quantiles of AA.5(a), chain of inequalities of AA.6(a); see §7). Its results are mostly standard transfers and cover the idealized convex problem, not the deployed detector: AA.8 and M1918-1920 admit that the implemented statistic is neither the certificate nor a $\chi^2$ field. Conjecture AA.1 (exact recovery on the continuum, amplitudes in $\mathbb R^6$) is not proven. The remainder $O(V/(f_{\max}d_{\min}))$ of (8b) is not bounded. Families C-II to C-V lack their own theory beyond their statement. |
| Experimental design | 20 % | 6.0 | 1.20 | Strengths: mismatched generator, station corrections from a catalogue with error, QuakeMigrate as published, controls that break the geometry, ablations within the same solver, hierarchical bootstrap and Holm. Weaknesses: operating point chosen by oracle on the evaluation in 2A-2C and in part (i) of 2D; five networks per cell and two seeds; rules decided on point estimates ("2 of 3 cells"); designs of 2B-2D informed by previous results and audits; pick-based comparators throttled by the picker (10-17 % of arrivals picked); QuakeMigrate calibrated only in weak clutter. The confirmatory plan G1' is sketched but is incoherent across sections (see §6). |
| Empirical evidence | 15 % | 4.5 | 0.68 | Everything is synthetic: one generator family, white noise, 1D medium, double-couple sources inside the detector subspace, sampling at 16 Hz, a 30×30 km domain, outside the regime of M1.5. The advantage over QuakeMigrate is robust over 30 networks, but it belongs to the timing-tolerant family, not to the method. There is no deployable operating point: with conformal thresholds, 1.2-2.2 false events per window with events. On real data only the polarization gate was measured, and it passes at 1 of 5 stations. |
| Internal coherence and writing | 10 % | 4.0 | 0.40 | ~105 000 words between M and S, with the same evidence restated at least six times (Abstract, 1.6, 4.3, 7.5, 13.0c, 20.0, 23). There are structural contradictions (what "the core" is, reciprocity and unique continuation as constraint generators, rows of Appendix D), obsolete text in S and R, and a misrepresented citation ([16]). The prose is dense and precise locally, but unreadable for a committee as a document. |
| Feasibility | 10 % | 5.0 | 0.50 | G1' is cheap in compute (140 GPU-h). By contrast, the 48-month plan for a single doctoral student includes an elastodynamic solver, SGT, Gaussian beams, a learned residual, optimal transport, SBI, single station, network design, two real domains and eight theoretical packages. The milestone of month 6 (RSNC decimation with MT-FW and BOMP-PK, orientation audit, REAL, GENIE, azimuth baseline and second generator) is overloaded. The risk of transferring the MT detector to the RSNC (polarization reliable at 1 of 5 stations) has no mitigation within the operator. |
| Reproducibility | 10 % | 7.0 | 0.70 | Verified: the hash chain of the 2D preregistration (prefixes of 431 and 436 lines and the complete file), the hashes of the 2D code and the declared edit of `analyze_2d.py`. There are bit-for-bit controls and an OpenBLAS kernel pinned per network. Limitations: 2B and 2C without a verifiable seal; an undeclared edit of `analyze_2b.py` after its first output; manifest paths (`results/task2a/`...) that do not match the actual location; no public repository or environment file cited in M18. |
| **Total** | | | **5.11 → 51/100** | |

**Score of a reasonable committee today: 62/100.** A doctoral committee in applied mathematics would reward the integrity, the preregistration discipline and the volume of pilots. It would penalise the length, the lack of focus and the fact that the announced methodological contribution is no longer distinctive. The foreseeable outcome is an approval conditional on cutting the document and the scope and on reformulating G1 and the title.

---

## 3. Ten weaknesses by severity, with concrete corrections

**D1 (critical). The announced methodological contribution is not distinctive and the distinctive one does not yet exist.** With the same peak search, BOMP-PK equals MT-FW: differences of -0.155 to +0.091 at 0.25, all hierarchical intervals containing zero, and the rule holds in 11 % and 0 % of the replicates (S3815-3828; `analyze_2d_out.txt` lines 83 and 146). Without a cap, MT-FW and BOMP-PK give the same recall in almost all cells (S3839-3844; Table Y.48, identical rows). The certified group BLASSO stays below with and without peak search (S3846-3853). Even so, the title (M2, M10) promises "a Convex Certificate", and G1 (M629-631, M737-740, M1110-1121) defines the detector as "the certificate of a group BLASSO", when AA.8 (S4870-4891) and M1918-1920 acknowledge that the statistic with peak search and cap is not that certificate. *Undemonstrated speculative assertion:* that the BLASSO with off-the-grid shifts "turns the heuristic peak search into a principled tolerance and outperforms it" (M173-176). **Correction:** (i) retitle without "Convex Certificate", for example "Timing-Tolerant Moment-Tensor Detection in Sparse Networks: a Preregistered Benchmark and an Off-the-Grid Formulation"; (ii) reformulate G1 on the timing-tolerant family with deflation and turn T2 into the primary methodological hypothesis H-M, with a preregistered and externally sealed rule against BOMP-PK on at least 30 networks per cell; (iii) position T2 against CAP (Zhu and Helmberger, 1996), Phạm (GJI 2024) and the joint MT-location inversion by stacking of P and S diffractions (PS-JMDS, Geophysics 86(6), 2021) before claiming novelty.

**D2 (critical). Exclusively synthetic evidence, from one generator family and outside the declared regime.** M218-225 and S2941-2946 acknowledge it: wavelets of 1.4-2.8 Hz, sampling at 16 Hz, a 30×30 km domain at 5-13 km depth and a 1D medium, versus the primary regime of 0.1-2 Hz and 120×120×40 km. The sources are double couples inside the 6-dimensional subspace of the detector (S2893-2894), the noise is white and the clutter is made of isolated Ricker pulses. **Correction:** run before the defence Tasks 1p (real noise) and 1v (3D Green); add non-DC sources and red microseismic noise; repeat Task 2A in the declared regime or redeclare the regime according to the evidence (local, short period).

**D3 (high). There is no deployable operating point.** In 2A-2C and in part (i) of 2D the operating point of each method is chosen on the evaluation windows (S2694-2696). The frozen rule 2A-e fails in its strict reading (realized rate of 2.3 times the budget, S2828-2839). The conformal threshold of 2D controls only empty windows. With it, MT-FW+SC emits 0.73-1.46 false events per window with weak clutter, between 3 and 8 times what QuakeMigrate emits (S3910-3914; verified in `analyze_2d_out.txt`, Table 7), and 78-92 % of those false events fall 5-20 km from a true event (S3920-3926). A good part of the recall gain of 2D-c is bought with those false events. **Correction:** a per-event calibration (threshold as a function of the residual energy after deflation, or per-event FDR control with conformal p-values conditioned on the presence of an event), calibrated on different networks and with windows that contain events; declare as primary synthetic endpoint the recall at realized false events per window with a frozen threshold, not the oracle.

**D4 (high). Weak statistical support at the replication unit.** 2A, 2B and 2D have five networks per cell and two seeds. The rules are decided on point estimates with "at least 2 of 3 cells". In 2A, at 0.1, four of twelve per-network intervals contain zero (S2874-2883). 2B with weak clutter passes at the minimum and holds only in 59-62 % of the replicates (S3087-3093). In 2D no rule survives Holm (M2621-2623; S3835-3837). **Correction:** bring all registered rules to 30 networks per cell (like 2C-c), decide on the lower bound of the hierarchical interval, set the sample size by power with the measured between-network $\sigma_d$, and report $\hat A_{12}$ per cell, which M14 promises and does not appear in Y.18-Y.21.

**D5 (high). Comparators at a structural disadvantage.** The frozen STA/LTA picks 9.6-16.8 % of the arrivals and GaMMA and PyOcto are throttled (S2855-2861); with optimistic picks, GaMMA is ahead at 0.1 in 4 of 6 cells with weak clutter (S2804-2808). QuakeMigrate was recalibrated only in weak clutter, in two cells, on 60 windows and without SC, with the LTA at the edge of the grid (S2308-2318, S2867-2871). Moreover, by construction, it cannot trade false alarms for recall on the FROC curve (S2357-2361), so the endpoint favours long-list methods. MALMI [110], [112], REAL and GENIE have not been run. **Correction:** recalibrate QuakeMigrate per regime and cell with the same budgeting as MT-FW; run MALMI and GENIE; include a picker within its domain (synthetic data with real noise or retrained PhaseNet); report an endpoint at a common false-event rate that does not depend on the length of the list.

**D6 (high). Integrity of the preregistration verifiable only for 2D.** Verified with `stat` in UTC. The hash chain of 2D is correct: the SHA-256 of the first 431 lines is `6e8510b3...`, that of the first 436 is `dc8cf708...` and that of the complete file is `b01a2eac...`, and `PREREGISTRO_hashes.txt` was born at 21:13:04.46. For 2B and 2C, by contrast, the registration was last rewritten on 27-09 at 03:54:55, so the file system cannot date a text within the file. The typed time of 2B, "14:30", coincides with the birth of `task2b_report.txt` (14:30:37), later than the first output of results (`analyze_2b_out.txt`, 14:26:15). The claim of R632-635 and M794-798 that the text of 2B predates 14:16:49 is an *Undemonstrated speculative assertion*. Moreover, `analyze_2b.py` was modified at 14:26:33, 18 s after its first output, and the output was rewritten at 14:26:43; the edit is not declared in `task2b_report.txt`. The registration file itself was born on 25-09 at 20:14:09, before the typed time of 2A (20:15:00) and even though it contains entries from 23-09, so the entries prior to 2A cannot be dated. **Correction:** reclassify 2B and 2C as "registered without a verifiable timestamp"; correct M790-798, which states that "for Tasks 2A and 2C the registered times precede every file of the task in the file-system record" (for 2C this is not verifiable); from now on, seal externally (signed commit in a public repository, OSF or an RFC 3161 timestamp).

**D7 (medium-high). Risk of transfer to the RSNC without mitigation within the detector.** MT-FW projects the three-component vectors onto the radiation subspace of a homogeneous medium. In the RSNC, single-station P polarization passes the gate at 1 of 5 stations and at CHI it is at chance level (S843-896). The orientation has not been audited. The velocity model does not support waveforms at 2 Hz (M2760-2765). The per-station reliability scheme, (9d)-(9e), was designed for the per-pick likelihood, not for the operator $G$ of MT-FW, which treats all stations equally except for the cap. **Correction:** introduce into $G$ a per-station orientation and reliability nuisance parameter (free horizontal rotation or weighting by $\hat r_s$); run MT-FW on catalogued events of the RSNC (recovery of location and mechanism) as a gate prior to G1'; restrict G1' to stations with audited orientation.

**D8 (medium-high). Scope unfeasible for a single doctoral student in 48 months.** Table 8 (M2893-2902) keeps the solver, SGT, beams, learned residual (up to 1 800 conditional GPU-h), optimal transport, UQ/SBI, single station, network design, two real domains and T1-T8. M2358-2361 itself acknowledges that most of the budget funds lines that do not decide the thesis. **Correction:** reduce the thesis to G1/G1' and T2, plus the observability certificate (B3, B1), and move families F, C-II and C-V, the learned residual and optimal transport to future work with one page at most.

**D9 (medium). Gap between theory and method.** The theory of Appendix AA is correct but covers the idealized convex problem ($S=0$, $\kappa=\infty$, $c=0$). For the deployed iteration it offers only AA.4(e)-(h): neither a rate nor monotone decrease with the cap active. Exact recovery is Conjecture AA.1. In the Proposition of 6.1 only the leading term is proven and the remainder is not bounded. The theory of C-II (certified exclusion) and C-V (admissible cone) is not developed. **Correction:** formulate and prove, for the model with bounded shifts $|\delta_{ib}|\le w$, existence, the optimality conditions, the certificate over shifts and a rate for SFW with an inexact oracle (extension of AA.4(c)); bound the remainder of (8b) at least in the homogeneous medium.

**D10 (medium). Coherence, length and citations.** There are contradictions and obsolete text (list in §5), and the evidence is restated six or more times. Citation [16] is misrepresented: M401-402 attributes to Puente Huerta et al. (2025) that GaMMA is "the most degraded", but the source points out that PhaseLink drops almost to zero and that GaMMA "outperforms PhaseLink but struggles". **Correction:** a main document of 15 000 to 20 000 words, with the evidence stated only once (Section 1.6) and referenced afterwards; correct [16]; execute and report the PRISMA protocol of Appendix E with counts, which M3326-3327 promises and does not appear.

---

## 4. Real strengths

1. **Demonstrable scientific integrity.** There are 35 registered errors (R, E1-E35), four positives withdrawn, registered rules that failed and are reported as such (2A-d with strong clutter, strict 2A-e, 2D-a) and a rewriting of the thesis based on its own negative result (Y.21a). This exceeds the usual standard of a doctoral proposal.
2. **Cryptographic traceability in 2D, independently verified.** The three registration hashes correspond exactly to prefixes of 431 and 436 lines and to the complete file; the hashes of `task2d_methods.py`, `task2d_run.py` and `task2d_validate.py` match `task2d_code_hashes.txt`; the only discrepancy (`analyze_2d.py`) is the declared edit.
3. **Exceptional numerical-engineering controls.** There is bit-for-bit regeneration of 1200 windows, an exact per-network control, detection and pinning of the OpenBLAS kernel (S3016-3035) and record-by-record equivalence of the certified oracle against the dense search (S4093-4141).
4. **Honest threshold-calibration study (Y.17).** It shows that the EEC of a $\chi^2_6$ field overestimates $P(M>u)$ by a factor of 20 to 60 and that strong clutter breaks every Gaussian calibration, and it adopts conformal calibration with a proven guarantee (AA.5(c)).
5. **Observability certificate with quantitative validation.** It gives $\rho=+0.9998$ over 120 geometries and extends the formula to gradient media through the take-off slowness (8d), with $\rho\approx0.999$ versus a 40-90 % degradation of the straight-ray analogue. It is a genuine mathematical contribution, publishable in GJI or Inverse Problems, independent of G1.
6. **Benchmark design that is hard to replicate in the literature.** Mismatched generator, station corrections estimated from a catalogue with error, unmodified QuakeMigrate 1.2.2, controls that break the geometry, ablation of the operator within the same solver, hierarchical bootstrap over 30 networks and Holm.

---

## 5. Inconsistencies with line number

| # | Location | Inconsistency |
|:-:|:--|:--|
| 1 | M2, M10 versus M1958-1960 and S3851-3853 | The title announces "a Convex Certificate", but the certificate covers only the dictionary that builds the heuristic oracle and the deployed detector is not certified over shifts. |
| 2 | M629-631, M737-740 and M1110-1121 versus M1918-1920 and S4418-4422 | G1 defines the detector as the certificate of the group BLASSO; the implemented statistic (peak search and cap) is neither that certificate nor a GLRT. |
| 3 | M63-64 versus M287-288 | The Abstract says that GL-FW-PK "stayed 0.11 to 0.30 below" without "in 11 of 12"; in one case it stays above (+0.013, S3792). |
| 4 | M401-402 versus the source [16] | The citation attributes to GaMMA the greatest degradation; the source attributes it to PhaseLink. |
| 5 | M557-559 and M3154-3156 versus M1526-1529, M1469-1477 and S1105 | Reciprocity and unique continuation appear as generators of the feasible set; 6.3 and V.2 say that reciprocity cannot be tested with station data, and 6.2 says that unique continuation is only motivation. |
| 6 | M1297 versus M742-744 and M752 | "This is the core of the proposal" (Section 6) contradicts "Everything else is subordinate [to G1]". |
| 7 | M2912-2914 versus M2926-2948 | "dropping the constraint theory, which is the contribution" contradicts the primary methodological and empirical contribution of 20.0. |
| 8 | M790-791 versus `stat` | "For Tasks 2A and 2C the registered times precede every file of the task in the file-system record": for 2C this is not verifiable, because the registration was last modified on 27-09 at 03:54:55 and its time is typed, like that of 2B. |
| 9 | M796-798 versus `stat` | "the preregistered status of Task 2B rests on file-system times only": the file system does not date the text of 2B; the typed time (14:30) coincides with `task2b_report.txt` (14:30:37), later than the results (14:26:15). |
| 10 | M3302 (Appendix D, G1') versus M2799-2809, M2709-2710 and M2791-2792 | G1' is defined "per window … in both seeds" on real data, where the rate is per day and there are no seeds; the list of comparators omits BOMP-PK, which M15 declares mandatory and M16 includes. |
| 11 | M2363 versus M2373-2374 | "Tables 3 to 5 constitute a planning model, not a measurement" contradicts "derives from the measured cost model of Tables 3 and 4". |
| 12 | M1711-1712 versus S125-126 | "degrading three times more slowly": the relative degradation is 66 % versus 29 %, a ratio of ≈2.3, not 3. |
| 13 | M1586-1587 | Broken sentence: "which is not assume and does not use". |
| 14 | M3100-3104 versus S3415-3419 | "exceeded … every comparator … in every cell": with Holm, in the 4-station cells with strong clutter the advantage over BOMP+SC (2024) and over QM-3C+SC, BOMP+SC and MTSTACK+SC (7331) is not significant. |
| 15 | S2115-2119 and S4495-4496 versus S4121-4125 | It is said that the certified bound in the pruned search is "not yet run" or that it is "registered work", but Z already measured it (it avoids 4 % of the evaluations). |
| 16 | S2644 versus `stat` of the registration | Registration time of 2A "20:15 UTC"; the registration file was born at 20:14:09. The registration times prior to 2D were typed (E34). |
| 17 | S3821-3822 versus `analyze_2d_out.txt` l.146 | "0.1 per cent of the replicates" versus "0% of replicates H" (rounding; minor). |
| 18 | R3-5 and table of R §2 | Paths `results/task2a/`, `results/task2b/`, etc. that do not exist; the files are in the root of `/home/claude/work`. |
| 19 | R108 versus R112-642 | "In the six cases": the register contains 35 errors. |
| 20 | R618-619 | "The English proposal was not edited in this update": obsolete, because M already incorporates the change. |
| 21 | M1419-1421 | The coercivity result of I ($\rho=+0.82$) is presented as one that "exceeds the registered B3 threshold", when B3 refers to the local margin (Appendix S); two endpoints are conflated. |

---

## 6. Verifications performed

**Numbers against result files (13 checks, all concordant except for the rounding of item 17 of §5).**

1. Table Y.31, weak clutter, 6 stations (7331): +0.49 [+0.25, +0.65] at 0.1 = `analyze_2a_out.txt`, Table 3.
2. Table Y.32, one-sided 8 stations (7331) at 0.1: 0.77/0.22, "NO" = Table 4 (thin).
3. 2A-d with strong clutter, "fails" (1/3 on 7331) = Rules (strong).
4. 2A-b, difference of +0.0513 at the limit = final line of `analyze_2a_out.txt`.
5. 2B with weak clutter: +0.0494 in the one-sided 7331 and rule held in 0.59/0.62 of the replicates = `analyze_2b_out.txt`.
6. Table Y.40, strong, 4 stations (7331): +0.099 [+0.029, +0.181], pHolm 0.044 = `analyze_2c_out.txt` l.177.
7. Table Y.44: -0.155 [-0.277, +0.032] and rule held in 11 % of the replicates = `analyze_2d_out.txt` l.83 and l.110.
8. Table Y.47: 10/45 = 0.222 at 6 stations (7331) at α = 0.1 = l.177.
9. Y.17, conformal: 0.007/0.037/0.095/0.239; 8 700-19 000 effective hypotheses versus 330 000-640 000 = `rft_report.txt`.
10. Y.16, K = 20: rms of 0.088 s and 4 stations with 0.51/0.64 = `station_corrections_report.txt`.
11. Y.15: MT-FW - QM-3C = +0.28 [+0.10, +0.43] = `qm_report.txt`.
12. Y.14: 0.420/0.143/0.041 arrivals picked = `gen_realistic_report.txt`.
13. Z: 215 982 of 225 444 evaluations and 13.1 s versus 55.9 s = `certified_fw_report.txt`.

Internal arithmetic verified: budget (1 668 ≈ 1 670 GPU-h; breakdown 140 + 605 + 708 + 215); Table 3 ($h$, $N_{\mathrm{dof}}$, steps, PFLOP, wall-clock time); constant 7.85 of (25), which with $n=100$ gives $\sigma_d=0.178$ for $\Delta=0.05$; quantiles of AA.5(a), $\chi^2_6$ with $1-\alpha/20925$ equal to 36.32 and 39.89 (the text says 36.3 and 39.9).

**File-system times (UTC, `stat`).** Registration born on 25-09 at 20:14:09 and modified on 27-09 at 03:54:55. 2A: the first file, `task2a_confirm.py`, is from 23:33:22 and the freeze of 2A-e (01:12:25) precedes the first evaluation (01:12:53); concordant. 2B: code at 14:16:49, first result at 14:26:15, `analyze_2b.py` edited at 14:26:33 and report at 14:30:37; the text of the registration cannot be dated. 2C: `task2c_run.py` at 15:36:47, after the typed time of 15:05; the text cannot be dated. 2D: registration sealed at 21:13:04, code from 02:24:59 to 02:29:23, hashes at 02:41:59, parts at 02:42:01 and result at 03:54:55; concordant and verified.

**References (5 checked).** [9] Shi et al., *Nat. Commun.* 17, art. 8140, 2026: correct. [109] Winder et al., *Seismica* 5(1), 2026, doi 10.26443/seismica.v5i1.1854: correct. [129] Vera Rodriguez et al., *GJI* 188(3):1188-1202, 2012: correct. [16] Puente Huerta et al., *Seismica* 4(2), 2025: bibliography correct, content misrepresented (§5, item 4). [19] VORA, arXiv:2607.10450: exists.

---

## 7. What would raise the score

**+5 points (up to ≈56), with rewriting and repositioning only, without new experiments:**

1. Retitle and reformulate G1 and the methodological contribution in terms of the timing-tolerant family with deflation, and declare T2 as an open methodological hypothesis with its rule (D1).
2. Resolve inconsistencies 1-21 of §5, correct [16] and cite and distinguish CAP, Phạm (2024) and PS-JMDS (2021).
3. Reclassify 2B and 2C as registrations without a verifiable seal and declare the later edit of `analyze_2b.py` (D6).
4. Reduce M to a maximum of 20 000 words and move families F, C-II, C-V, the learned residual and optimal transport to "future work" (D8, D10).
5. Redo the G1' row of Appendix D (per day, without seeds, with BOMP-PK) and close the sample design of G1' with the between-network $\sigma_d$.

Estimated effect: coherence from 4.0 to 6.5, feasibility from 5.0 to 6.5, novelty from 4.0 to 4.5 and reproducibility from 7.0 to 7.5, that is, +0.55 → ≈56/100.

**+15 points (up to ≈66), with new results:**

1. Implement the group BLASSO with off-the-grid per-station and per-phase shifts (SFW with debiasing refit) and test with a preregistered and externally sealed rule, on 30 networks per cell, whether it outperforms BOMP-PK and MT-FW on the windows of 2A. If it does not, report the null and present the contribution as benchmark plus theory (D1).
2. Extend AA to that model: existence, certificate over shifts, duality gap and rate with an inexact oracle (D9).
3. Run Tasks 1p and 1v and add non-DC sources and a second independent generator (D2).
4. Obtain a deployable operating point: per-event calibration on separate networks, with the realized false events per window over all windows within the budget (D3).
5. Run a real pilot of G1' on a subset of the RSNC with audited orientation, a threshold frozen on noise windows and false events per day, with MT-FW and BOMP-PK against QuakeMigrate and MALMI, and include per-station reliability and orientation in $G$ (D5, D7).

Estimated effect: novelty from 4.0 to 6.0, evidence from 4.5 to 7.0, design from 6.0 to 7.5, theory from 5.5 to 6.5, plus feasibility and coherence; ≈+1.5 → ≈66/100.

---

*Labels applied.* **Undemonstrated speculative assertion:** M173-176 (the BLASSO with shifts outperforms the heuristic), M794-798 and R632-635 (the text of 2B predates 14:16:49), M1709 ("Passing all of them simultaneously requires the clutter to be, in effect, an earthquake"), M2367-2368 ("the ordering between architectures … is robust to these factors"). **WEAK BASELINE:** QuakeMigrate calibrated only in weak clutter and two cells; GaMMA and PyOcto with a throttled picker; MALMI and GENIE absent. **UNDEMONSTRATED CONTRIBUTION:** MT-FW versus [129] with the same tolerance. **EXPERIMENTAL OVERFITTING:** a single generator, two seeds. **OVERCLAIMING:** title.
