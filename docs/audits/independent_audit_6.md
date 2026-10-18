# Independent audit 6 of the undergraduate thesis "Timing-Tolerant, Bounded-Influence Moment-Tensor Detection of Overlapping Earthquakes"

Examiner: independent committee member (did not write any of the work). Date of audit: 17-18 October 2026.
Object: repository clone at /home/claude/work/repo (github.com/juanjgc0609/timing-tolerant-mt-detection), last commit
095f630 (2026-10-17 21:10:47 -0500). Earlier audits (docs/audits/, docs/audit_status.md) were not read.

Main document examined: docs/thesis/Thesis.md (9,330 words, ASCII, no figures). Supporting material examined: the
evidence table and its generator, Appendices Y.18 to Y.33 and AA.9 of the evidence volume, the error register, the
sealed preregistration log and its translation, the result and verification files of every task cited in the thesis,
the R1 code and freeze list, README.md, data/README.md, docs/executive_summary.md and the git history.

Computation done for this audit (all light): recomputation of the hash chain against byte prefixes of the sealed log;
check of the append-only property across all seven commits of the public history; line alignment of the translation;
run of code/thesis_evidence_table.py; re-run of code/r1/r1_analyze.py on the archived R1 outputs; re-run of
code/theory_tail_check.py; check of the 297 file hashes of code/r1/FREEZE_R1.sha256 against the repository; reading
of the deployed MT-FW code (code/explore_contrib/c9_lib.py, c8_lib.py, c10_lib.py) and of the R0/R1 labelling code.

---

## 1. Hash chain and translation

### 1.1 Chain

File: preregistration/original_sealed_es/PREREGISTRO_hashes.txt (68 lines) against PREREGISTRO_1j_1k.txt (120,847
bytes, 1,220 lines, ends with a newline).

| Check | Result |
|:--|:--|
| Number of SHA-256 entries | 34 |
| Every entry equals the SHA-256 of a newline-terminated byte prefix of the log | yes, 34 of 34 |
| Prefix lengths non-decreasing | yes (39,694 to 120,847 bytes; strictly increasing) |
| Last entry equals the digest of the whole file | yes (dd785c02...a3ff7) |
| Registration of R1 (3b238aa4...402e416c) | prefix of 113,828 bytes, line 1152, timestamp line 63 = Sat Oct 17 18:49:11 UTC |
| Amendment 1 of R1 (d5d420a6...67125c1) | prefix of 117,909 bytes, timestamp 19:10:24 UTC |
| Result of R1 (dd785c02...) | whole file, timestamp 02:09:37 UTC on 18 October |

The chain is valid. Two observations on the format, neither of which breaks it:

1. The file is not self-describing. Reading the git history shows that each timestamp line precedes the hash it
   dates (for example, the log in commit 4d5a931 has hash 53b71c22..., which is line 62, preceded by the timestamp
   07:26:05 UTC of line 61). With that reading, the first hash (line 1) has no timestamp and line 32 (Mon Oct 5
   02:36:52 UTC) is an orphan timestamp with no hash after it. preregistration/README.md says only "UTC time and
   SHA-256 ... after each append". A reader who pairs hash and time the other way gets every time shifted by one
   entry. This should be documented.
2. It is a list of digests of successive prefixes, not a cryptographic chain in which each entry commits to the
   previous one. Because each digest covers the whole prefix this is equivalent for tamper evidence of the log, but
   the hash file itself can be rewritten wholesale. Its protection rests entirely on where the digests were
   disclosed outside the container (see 5.6).

### 1.2 Append-only property across the public history

For each of the seven commits (c5f45bd of 11 October to 095f630 of 17 October) the committed log and the committed
hash file are byte prefixes of the current ones (14 of 14 checks). So nothing in the log was rewritten after 11
October 2026, the date of the first public commit. Everything before that date (Tasks 1j to 2L and the R0
registration) is protected only by the self-timestamped chain and by the chat record. Note also that git author
dates are set by the committer and are not trusted timestamps; only the server-side push time on GitHub is external.

### 1.3 English translation

PREREGISTRATION_LOG_EN.txt has the same 1,220 lines as the sealed file, ASCII only, and blank lines at identical
positions (0 mismatches). Lines 1 to about 776 correspond one to one. From about line 777 the paragraphs are re-wrapped
inside each paragraph: in 243 lines the multiset of numeric tokens differs between the two files because a few words
moved to the neighbouring line (examples: lines 777, 790, 951, 1121, 1174). Paragraph boundaries and the starting line
of every entry agree (R1 registration line 1065, amendment line 1154, result line 1193 in both). The statement of the
README that "a line cited in one is the same line in the other" is true to within one line inside a paragraph, not
exactly. I found no substantive translation error in the sampled lines; the R1 entries in results/taskR1/log_entries/
are contained verbatim (whitespace-normalised) in the sealed log, with the placeholder __FECHA__ filled by the sealed
times.

Verdict on Section 1: the chain verifies; the format and the one-line drift of the translation should be documented.

---
## 2. Number checks

### 2.1 Table 1

`python3 code/thesis_evidence_table.py` run from the repository root exits with status 0 and its output is
byte-identical to lines 352 to 369 of Thesis.md (diff empty). Each of the 16 rows was also checked by eye against its
source line; all values are on the cited lines.

Two weaknesses of the generator itself. (a) It checks that each value is a substring of the source line, so "MET"
would also be found in a line that says "NOT MET", and "0/3" in "10/30". For E1 and E8 I checked that the lines end in
"-> MET" and not "NOT MET"; they do. (b) It certifies quotation, not computation: it proves that the thesis copies the
analysis files correctly, not that the analysis files are right. Recomputation is a separate question (Section 6,
reproducibility).

### 2.2 Re-runs

| Item | Command | Result |
|:--|:--|:--|
| R1 analysis | `python3 code/r1/r1_analyze.py eval results/taskR1/run/hyb_eval.json results/taskR1/run/hyb_cal.json OUT` (env of /home/claude/work/r1/env.sh) | 3 s; output byte-identical to results/taskR1/run/analysis_eval.txt, including the bootstrap |
| Proposition 4 check | `python3 code/theory_tail_check.py OUT` | 38 s; output byte-identical to results/theory/tail_check_out.txt |
| R1 freeze | SHA-256 of the 297 file entries of code/r1/FREEZE_R1.sha256 against code/ and the data folders | 297 of 297 match; hash of the list = 1669903f...28a2a3 as quoted in amendment 1 |
| Constants of Prop. 4 | scipy: Pr(chi2_3 > 36) = 7.49e-8; 14 x 9 x p = 9.4e-6; 7 x 640 x p = 3.4e-4 (one window in 2,981); bound at p = 1e-4 = 28.49 | all as stated |
| Prop. 1 numbers | 2c = 12; 12 sqrt(log(1 + (12500/6)^2)) = 46.9; sqrt(2 x 77.35) = 12.4, sqrt(2 x 74.91) = 12.2 | as stated (thresholds from results/task2M/analyze_2m_out.txt lines 162 and 331) |

### 2.3 Spot checks of numbers in the text (beyond Table 1)

| # | Claim in Thesis.md | Source checked | Status |
|:--|:--|:--|:--|
| 1 | 2A: MT-FW+SC ahead by at least 0.05 at both budgets in all twelve cells | results/task2a/analyze_2a_out.txt, rules 2A-a and 2A-b, 3/3 cells both seeds both regimes | correct |
| 2 | 2C-a: +0.077 to +0.519 | results/task2c/analyze_2c_out.txt lines 44, 83 | correct |
| 3 | 2D: MT-FW without its own peak search lost +0.064 to +0.412 | analyze_2d_out.txt lines 85, 148 | correct as point values; not stated that with Holm 2D-b is not met (lines 210, 271) |
| 4 | 2E: +0.000 to +0.039, every interval contains zero | analyze_2e_out.txt lines 67, 160 | correct |
| 5 | 2F: +0.043 to +0.074, surviving Holm | analyze_2f_out.txt lines 26-31 | correct |
| 6 | 2F with origin time: +0.024 to +0.045 | results/task2f/verif_2f/task2f_verification.txt line 70 | correct (verification, not the registered analysis) |
| 7 | 2I: +0.058 and +0.053; Huber about nine times as costly | analyze_2i_primary_out.txt line 82-83; task2i/README.txt line 116 (0.927 vs 0.103 s) | correct |
| 8 | 2L: +0.057, +0.062 (standard), +0.038, +0.030 (dense); worse terms in about 22 per cent of networks | analyze_2l_out.txt; verif_2l line 122 | correct |
| 9 | 2M secondary: +0.057 and +0.120 | analyze_2m_out.txt lines 344, 346 | correct |
| 10 | 2M breakdown: +0.033 at night, +0.300 between 7.5 and 15 h | task2M/README.txt line 19; posthoc/blocks_out.txt line 2 | correct (post hoc, labelled as such) |
| 11 | R0: at most 0.022 against up to 0.667 | results/taskR0/analysis_eval.txt LOC+T 30/h | correct for M2, the registered detector; M1 (also MT-FW, cap 6) reached 0.044 on N8 |
| 12 | R0 closure: +0.000 [-0.023, +0.021] (LOC+T), -0.050 [-0.101, -0.013] (LOC) | run_close/analysis_eval.txt lines 13, 16 | correct; it is M1 - M3, not M2 - M3 |
| 13 | R0 timing residuals 0.35 to 0.41 s against about 0.14 s synthetic | Appendix Y.30 line 25 of the section | consistent with the evidence volume (not recomputed) |
| 14 | C18: AUC capped 0.39 to 0.49, uncapped 0.69 to 0.74 | Appendix Y.32 (0.49/0.39/0.43/0.43; 0.71/0.69/0.72/0.74) | correct |
| 15 | C18: layered model 0.45 against 0.37 to 0.39 within 0.25 s | Y.32 | correct for N8 and U8 |
| 16 | C18j: 4 and 5 of 16 development events at 30/h | c18j_detect_f2_U8_analysis.txt (NOCAP 0.312 = 5/16) | correct |
| 17 | E52: 0.78 instead of 0.94 of an arrival's energy | results/code_audit/floor_bias_bh.txt | correct (0.783, 0.938) |
| 18 | R1 simulation: McNemar false-positive rate 0.000 to 0.018 at zero difference; power 0.13-0.17 and 0.49-0.69 (40 targets), 0.35-0.53 and 0.92-0.99 (100 targets) | results/taskR1/sim_r1_rules_out.txt | correct |
| 19 | R1 simulation: "with 9 to 14 targets a percentile bootstrap gives a false-positive rate of 5 to 13 per cent" for the -0.20 margin | same file, column NI(0.20) at delta -0.20 | DISCREPANCY: 0.046 to 0.128 with 9 targets, 0.034 to 0.040 with 14 targets; the range over 9 to 14 targets is 3.4 to 12.8 per cent. The registration's wording ("with 9 events") was correct; the thesis generalised it |
| 20 | R1: HYB 23, BASE 22, BESTPICK 26, P1 1 of 54; McNemar p 0.50 and 0.97; JMA 9/10, 9/10, 10/10, 1/10; CAL threshold gives 41.8 and 12.9 false per hour | re-run of r1_analyze.py | correct |
| 21 | R1: "The block-bootstrap intervals were +0.018 [+0.000, +0.058]" | analysis_eval.txt | minor: +0.018 is the bootstrap median; the point difference is +0.019 |
| 22 | R1: development BASE 16 of 30, lags 0.265 s and 0.244 s | amendment 1; dev/analysis_cal_dev.txt | correct as development values; the frozen run gave 17 of 30 and a P lag of 0.272 s (verification D2); the thesis does not say so |
| 23 | R1: with EVAL from 880 s, 22 against 21 of 52 | verif_r1 section 5 | correct (verification) |
| 24 | R1: only 29 of 54 targets detected by any method | verif_r1 section 7 | correct |
| 25 | R1 registration 18:49 UTC; registration pushed about one minute after the run started | hash file line 63; commit 0374b2f at 19:11:45 UTC; run start 19:10:33 UTC | correct |
| 26 | Fifty-six errors E1 to E56 | docs/artifact_manifest_and_error_log.md | correct (56 distinct headings) |
| 27 | "Eleven registered tests ... each with a sealed registration, frozen code, a single run and an independent verification" (Sections 1.3, 7.3) | log entries; results/task2b, task2c, task2d (no verif_* folders); code/ (code hash lists exist only for 2D to 2L); executive_summary.md section 7 ("2B, 2C, 2D: self-checks, no independent verification") | DISCREPANCY: false for 2B, 2C and 2D (no independent verification) and for 2A to 2C (no frozen code hash list, registrations before the chain). Ten earlier registered tasks of the same project (1j, 1k, 1o, 1r, 1s, 1t, 1u, 1w, 1x, 1y; 1x "not met in any condition") are not counted |
| 28 | Section 6.2: "Each task used two seeds never used before" | log entries of 2B, 2C, 2D, 2E (lines 325-332, 350-360, 409-416, 449-456) | DISCREPANCY: 2B, 2C, 2D and 2E re-use the 1,200 windows of 2A (seeds 2024 and 7331), regenerated bit for bit, after the 2A (and 2C) results were known. Only 2A, 2F, 2I and 2L used fresh seeds. The 2A seeds had also been used by Tasks 1j and 1k with other generators |
| 29 | Section 5.1: "the hash is sent to a third party (the chat record)" | preregistration/README.md ("sent to the author in the conversation"); registration text of R1 ("delivered to the user through the chat") | DISCREPANCY in wording: the recipient is the author; the chat platform holds the record. It is not a third-party timestamp |
| 30 | Section 6.2: "no deployable threshold was frozen" (2A) | rule 2A-e and log2a_freeze.txt (thresholds frozen on other networks); Y.18a | inconsistent: frozen thresholds were tested in 2A-e (met per cell, not met in the strict reading, realised rate up to 2.3 times the budget) |

Summary: every quantitative result that the thesis quotes from a result file is copied correctly (Table 1 and 26
other values). The discrepancies are in descriptive statements about the protocol (rows 27 to 30) and one
generalisation of a simulation result (row 19). The protocol discrepancies matter more than their size suggests,
because the protocol is the first claimed contribution.

---
## 3. Mathematics (Section 4 of the thesis and Appendix AA.9)

### 3.1 Statements and proofs

The deployed statistic was checked against the code: `c9_lib.bases` takes the left singular vectors of the radiation
operator (P and S blocks with geometric spreading and radiation factors), dropping columns below 1e-8 of the largest
singular value, so Q has orthonormal or zero columns as Lemma AA.9.1 assumes; the score is
`(y*y).sum(-1)/2` with `y = Q^T v`; the cap is applied to the norm of each three-component block after the
argmax-energy sample is chosen (`c8_lib._leaf`), which is what the propositions assume.

| Result | Statement correct? | Proof complete under stated assumptions? | Assumptions stated? | Comment |
|:--|:--|:--|:--|:--|
| Lemma AA.9.1 (contraction) | yes | yes | yes | elementary (orthogonal projector plus reverse triangle inequality) |
| Prop. AA.9.2 (influence of one block) | yes | yes; part (c) correctly shows linear growth for the identity | yes | the bound r(a)+r(b) is loose; the sharper ||phi(v_k) - phi(v'_k)|| is available and would change the numerical reading (3.2) |
| Prop. AA.9.3 (saturation) | yes | yes (direct computation) | yes | it holds only when every block of the hypothesis is at or above the cap; the text then applies it to "most stations", which the proposition does not cover |
| Prop. AA.9.4 (log map keeps ordering) | yes | yes (one line) | yes; AA.9 itself calls the equal-norm hypothesis restrictive | with equal block norms there is no inter-station amplitude pattern, so the result concerns the ordering of events by overall size, not the "amplitude pattern across stations" that the text of Prop. 2 says is lost. The two narratives are conflated in Section 4 |
| Prop. AA.9.5(a) | yes | yes (union bound over n_b L candidate samples, chi2_3 tail) | yes | |
| Prop. AA.9.5(b) | yes | yes. I checked each step: the map from white noise to the stacked candidates has operator norm at most kappa = max|K(f)|/||k||_2 (restriction of a convolution, kernel normalised to unit norm); the radial maps are 1-Lipschitz (for l_c, l_c'(t) <= (1+t^2/c^2)^(-1/2) follows from log(1+u) >= u/(1+u), and l_c(t) <= t); the read score is bounded by the maximum over the L^n_b choices; each f_s is kappa-Lipschitz, hence sub-Gaussian with proxy kappa^2 (Tsirelson-Ibragimov-Sudakov); E f_s <= sqrt(E f_s^2) <= sqrt(tr Q^T Q) <= sqrt(6) needs the transformed blocks centred, isotropic, of second moment at most 3 and uncorrelated across blocks, which the independence assumption gives; the maximal inequality adds kappa sqrt(2 n_b log L) | yes, including the separation of P and S of one station by more than the kernel support | correct and standard. The constant 28.5 against an empirical 7.18 shows it is a guarantee, not a calibration, as the text says |
| Prop. AA.9.5(c) | yes | yes (Lemma AA.9.1 over m blocks, each contributing at most 2c) | partly: the robust per-window normalisation of the data is itself affected by a transient, and this is not covered | |

The mathematics is correct. It is also elementary: Lemma 1, Propositions 1 to 3 are one- or two-line consequences of
the definitions, and Proposition 4 is a textbook application of Gaussian concentration and a maximal inequality. That
is appropriate for an undergraduate thesis, and the statements are honestly scoped ("The theory does not cover the
greedy deflation, the grid, the station terms or non-Gaussian, coloured noise").

The numerical check (code/theory_tail_check.py) reproduces. It tests the bound under the bound's own assumptions:
independent white-noise blocks, a random orthonormal Q (not the radiation subspace of any node), one hypothesis, no
scan over hypotheses, no deflation. It is therefore a check of the arithmetic and of the looseness of the bound, not
evidence about the null law of the deployed detector. The thesis says this in substance.

### 3.2 Is the connection to the experiments justified?

Four points where the text says more than the mathematics.

1. "The cap costs nothing in Gaussian noise (Proposition 4a)" (Thesis Section 4). Proposition 4(a) concerns the null:
   with no event, cap and identity coincide except with probability 3.4e-4 per window. It says nothing about power.
   With events present the cap discards amplitude even in Gaussian noise, and the registered data show a cost: in
   Task 2F, standard windows, strong clutter, MT-FW (cap) against MT-FW-NOCAP at 0.25 false events per window gives
   0.379 against 0.405 (seed 4242) and 0.351 against 0.383 (seed 8888), that is -0.026 and -0.032
   (results/task2f/analyze_2f_out.txt lines 19-20); and Task 2I found that the cap falls behind Huber when the
   amplitude spread is large. Appendix AA.9 has the correct wording ("costs nothing in false alarms"); the thesis
   drops the qualifier.
2. The numbers under Proposition 1. A transient can move the capped sigma "by at most 12" while "a real event at the
   detection threshold of Task 2M has sigma of about 12" (12.2 to 12.4). The bound and the threshold are the same
   size, so the bound does not by itself show that one transient cannot produce a capped detection: a noise
   hypothesis with sigma near its mean (3.5 in the simulation) plus 12 exceeds the threshold. Appendix AA.9 concludes
   that one transient "cannot decide [the ranking] of a capped one by itself", which does not follow. With the
   sharper bound (a clean noise block of norm about 2 to 3 replaced by a capped block of norm 6, change at most about
   8.5) the conclusion is borderline, and with m = 2 contaminated blocks (P and S of the same station, the usual case
   for a long transient) the margin of Proposition 4(c) is 17, above the threshold. The mechanism is plausible and the
   experiment supports it; the inequality does not prove it.
3. "Proved in the stated model (Propositions 1, 2, 4)" for the statement "The cap trades robustness to transients for
   saturation in continuous seismicity" (Table 2). What is proved is a bounded influence and a saturation identity
   that holds when every block saturates. That continuous seismicity puts every block of every competing hypothesis
   above the cap is an empirical premise; C18 measured its consequence (AUC) but not the premise (fraction of
   saturated blocks). "Proved" should be "explained, under a premise measured indirectly".
4. "as Proposition 4 predicts" for the concentration of the 2M gain in hours with transients (Section 6.3). The
   breakdown was written after the result and AA.9 after the breakdown; this is a post hoc explanation, not a
   prediction.

A fifth point is not a statement of Section 4 but bears on what it explains. In the deployed MT-FW the moment tensor of
each accepted atom is fitted to the capped read vectors and its prediction is subtracted from the uncapped data
(`c10_lib.mt_loop`: `vsel = _vectors(R, tpk, capv)`, `lstsq(G, vsel)`, `_deflate(R, G @ m, tp)`; the same in
`c9_lib.mtfw_c9` and `c8_lib.mtfw_general`). For an arrival of hundreds of noise deviations, as in real aftershocks,
the subtraction removes a prediction of cap size and leaves almost the whole arrival in the residual. Appendix AA.5(g)
acknowledges that "a capped event is not fully explained away", but the thesis neither states this in Section 3.3
("its least-squares moment tensor is estimated from the read vectors") nor lists it among the causes of the R0
failure and of the ghost detections in Section 7.2, where the ghosts are attributed to the point-source prediction of
the S coda. For R0 (capped M2) under-subtraction is a direct candidate mechanism; for R1 (uncapped P1) it does not
apply.

Verdict on Section 3: statements correct, proofs complete under stated assumptions, assumptions stated, depth modest.
The link to the experiments is qualitative and post hoc, and three sentences (the "costs nothing" claim, the
transient-versus-threshold reading, "proved" in Table 2) overstate it.

---
## 4. Claims against evidence

### 4.1 What the thesis gets right

The negative results are reported as negative, prominently and without hedging: the abstract states the Ridgecrest
recall (at most 0.022 against 0.667), the Hi-net recall (1 of 54) and the null hybrid result (23 against 22, p =
0.50); Section 1.3 states explicitly what is not claimed (that pick-free detection beats pick-based association on
real data, or that MT-FW beats BOMP-PK). The refutation of G1-M (2E) is reported as a refutation. Every major
qualification of 2A (STA/LTA picker, operating points chosen on evaluation, GaMMA ahead with optimistic picks) is in
the text. The hypotheses are restated after the real-data failures (Table 2), and the "conclusion" that the velocity
model does not matter, later reversed (E50), is told as a lesson. This is much better than the norm.

### 4.2 Overclaims and inaccuracies

| # | Statement (location) | Evidence | Assessment |
|:--|:--|:--|:--|
| C1 | "Eleven registered tests ... each with a sealed registration, frozen code, a single run and an independent verification" (1.3, 7.3; README "checked by an independent verifier") | 2B, 2C, 2D: no independent verification (no verif_* folder; executive_summary.md section 7 says "self-checks in its report; no independent verification"); 2A, 2B, 2C: no code hash list; 2A to 2C registered before the chain (declared in 5.1). Ten earlier registered tasks (1j to 1y) omitted from the count | false as written; the first claimed contribution is overstated |
| C2 | "Each task used two seeds never used before" (6.2) | 2B to 2E re-use the 2A windows | false for four of eight synthetic tasks; these four are same-data follow-ups designed after seeing earlier results on those windows |
| C3 | "the hash is sent to a third party (the chat record)" (5.1) | the hash is sent to the author through the chat | inaccurate; it is a record held by the chat platform, not a third-party attestation |
| C4 | "the gain belongs to the per-station bound and not to the particular pursuit, since block orthogonal matching pursuit with the same bound recovers almost all of it" (1.3, 6.3) | cited evidence is C18f, exploratory, one calibration seed, MT minus BOMP-CAP6 +0.026 [-0.009, +0.081] (an interval that does not exclude a difference of 0.08). Stronger evidence exists and is not cited: the registered 2F run contained BOMP-PK-CAP+SC as a descriptive arm, and MT-FW+SC minus BOMP-PK-CAP+SC is +0.004 to +0.015 in dense windows (analyze_2f_out.txt lines 60-91) | conclusion probably right, but argued from the weaker source. Conversely, the headline contrasts E4, E5, E8, E9 (capped or Huber MT minus uncapped BOMP-PK) confound the bound with the pursuit; in 2M no registered arm isolates the cap, so "Bounding the influence of each station added ... +0.228" (abstract) attributes the whole registered contrast to the bound on exploratory grounds |
| C5 | "The cap costs nothing in Gaussian noise" (4) | 2F standard windows, strong clutter: cap -0.026 and -0.032 against no cap at 0.25 | overstated (3.2, point 1) |
| C6 | "Proved in the stated model" (Table 2, row "new") | see 3.2, point 3 | overstated |
| C7 | G1 "Retained only as a conditional statement: it can hold where the travel-time error is below the reading tolerance, which neither real dataset satisfied" (Table 2) | the timing residual of the pick-free readings on the Hi-net swarm after amendment 1 is not reported anywhere I could find (Y.33, verif_r1); the R1 failure is attributed to S-dominated energy and coda ghosts, not to timing | the R1 half of "neither" is unsupported; and the conditional form is close to unfalsifiable as stated, because the only data that satisfy the condition are synthetic data generated to satisfy it |
| C8 | Section 7.2, ghosts: "the far-field point-source prediction of a 2 Hz wavelet does not predict a real S coda" | for the capped detector the deflation also subtracts a capped-amplitude prediction (3.2, fifth point) | incomplete diagnosis for R0 |
| C9 | "the pick-free detector added one event in fifty-four" (6.6) | verif_r1 4.4: the single HYB-only target is a P1 atom with 3 supporting picks | correct; good that it was traced |
| C10 | 2D "showed that the peak search, not the pursuit, carries most of the advantage" (3.2) | 2D-b not met under Holm (one of six cells per regime rejects); the exact-time arm was handicapped by E52 | the thesis discloses E52 but not the Holm failure; "showed" is stronger than the registered analysis |

### 4.3 Stale or contradictory documents

| Document | Statement | Current state |
|:--|:--|:--|
| docs/executive_summary.md | title "Epicentre determination and seismic source attribution as an inverse problem with constraints generated from elastodynamics"; section 2 "Pose epicentre determination as constrained optimisation where the feasible set is generated from the elastodynamic equation"; section 6 still argues G1 from synthetic data; section 9 ends "The forty-two errors are in ..." after saying fifty-six | the top status paragraph is current; the body is the doctoral-proposal framing. README calls it the "decision document" |
| data/README.md | Hi-net swarm "held out for a future registered real test ... Its waveforms have not been read by any analysis of this project" | R1 has read and spent it |
| data/README.md vs README.md | MatchLocate2 data "GPL-3.0; NIED Hi-net data" vs "under the NIED terms of use" | inconsistent licensing statement |
| preregistration/README.md | "Every confirmatory task (2A to 2L, R0, 2M)" | omits R1 |
| Thesis 6.2 vs Y.18a | "no deployable threshold was frozen" vs the frozen-threshold rule 2A-e | inconsistent |
| Thesis 5.1 | the freeze covers "every module it can import, the data files and the package versions" | the R1 freeze omits the PhaseNet source, r0_verif/pydeps, hyb/pydeps_pyrocko and the PyOcto model file (verification D5, Y.33); the thesis does not repeat this |

Verdict on Section 4: the scientific results are reported honestly, including the failures; the protocol is described
more uniformly than it was carried out, and two explanatory claims (C4, C5) lean on weaker or post hoc evidence.

---
## 5. Protocol of Task R1

### 5.1 Timeline reconstructed from the hash chain, file times and git

| UTC (17-18 Oct) | Event | Source |
|:--|:--|:--|
| 18:49:06 | registration text written (hyb/r1_reg_es.txt) | file mtime in /home/claude/work |
| 18:49:11 | registration hashed (3b238aa4...) | hash file line 63 |
| 18:52 to 19:09 | development on CAL: P1 v1 chunks (18:56), PhaseNet and associators on CAL (19:01), first hybrid and analysis (19:02-19:03), terms check (19:04), second P1 and hybrid (19:07), r1_lib.py last modified 18:59, r1_run.py and r1_hybrid.py 19:05, r1_queue.sh 19:09 | mtimes of /home/claude/work/r1/ and r1/dev/ |
| 19:10:24 | amendment 1 hashed (d5d420a6...) | hash file line 65 |
| 19:10:33.44 | FREEZE_R1.sha256 rewritten with identical content | verif_r1 1.2 |
| 19:10:33.47 to 19:20:56 | single run | verif_r1 1.2; run/queue.log |
| 19:11:45 | commit 0374b2f (registration, amendment, freeze) | git (committer date) |
| 19:25 to 02:08 | independent verification | r1_verif_indep file times |
| 02:09:37 | result entry hashed (dd785c02...), already containing the verification findings | hash file line 67 |

No file in the working folder that refers to the Hi-net waveforms has a modification time before the registration,
other than the registration text itself. This is consistent with the claim that no waveform was read before 18:49
UTC, with the caveat that modification times record only the last write.

### 5.2 What the documents disclose

The thesis and Appendix Y.33 disclose: the amendment and its reasons (P terms of +1.66, -1.70, +1.60 s caused by the S
arrival entering the P window at the near stations); that the amended P1 terms depend on the CAL picks; the
development read of [840, 880) s of EVAL (E56) and the false sentence of the amendment; the sensitivity analysis
without those 40 s (22 against 21 of 52); the choice of every operating point on EVAL; that the CAL threshold of HYB
gave 41.8 false detections per hour on EVAL; the ceiling of 29 of 54 detectable targets and the two M >= 2 targets
found by no method; the push of the registration about a minute after the run started; the low power and the
prediction that R1-a would very probably not be met. The independent verification reproduced every number with its
own code; my own re-run of the analysis reproduces the output byte for byte. On these points the protocol is
exemplary for an undergraduate project.

### 5.3 What the documents do not disclose, or disclose only outside the thesis

1. **The location criterion is not tested.** The 140 truth epicentres lie within 3.1 km of the frame centre (verif_r1
   6.1, section 7), and a hit needs an epicentre within 5 km. Any detection placed within about 2 km of the swarm
   centre is spatially compatible with every truth event, so LOC+T is in practice a test of origin time only, and
   "detect and locate" is not tested on this dataset. Not stated in the thesis.
2. **Lenient false-detection rule.** In `r0_analyze.label` and `r1_analyze.operating` a detection is false only if it
   is within tolerance of no truth event; a second detection of an already-hit event, or one near an event below
   magnitude 0.5, is neutral. With 84 truth events in the 1,120 s of EVAL, a detection placed on the swarm falls
   within 0.5 s of some truth event by coincidence about 7.5 per cent of the time (84 x 1 s / 1,120 s). The rule is
   the same for all methods, but it favours methods that emit many detections (HYB 277, P1 269 against BASE 29). The
   thesis states the neutral rule for small events but not the duplicate rule.
3. **Freeze coverage.** The PhaseNet source, the extra Python packages on PYTHONPATH (r0_verif/pydeps,
   hyb/pydeps_pyrocko) and r1/pyocto_jma2001.model were not hashed (verif_r1 D5; disclosed in the result entry and in
   Y.33, not in the thesis, whose Section 5.1 says the freeze covers "every module it can import").
4. **Run values differ from the amendment values.** The thesis quotes the development values (BASE 16 of 30, lag
   0.265 s); the frozen run gave 17 of 30 and 0.272 s because PhaseNet ran on the full record (verif_r1 D2). Harmless,
   but the thesis should quote the run.
5. **The result was sealed after the verification.** The result entry was appended seven hours after the run, with
   the verifier's findings already in it. The analysis is deterministic, so no number can have moved, but the chain
   holds no sealed record of the result as first seen.
6. **Power and the value of the test.** The development effect was one or two events in sixteen (0.06 to 0.12). From
   the registered simulation, with about 55 targets the McNemar test has power of roughly 0.15 to 0.25 at +0.10. The
   last unread real dataset was therefore spent on a test that would most likely fail even if the hybrid worked as
   in development. The registration says this, and the thesis does not misread the outcome ("no gain was
   demonstrated"), but "The prediction was right" presents an expected non-rejection as a confirmation; it carries
   almost no information about the hybrid.
7. **Weak baseline.** All ten GaMMA-1D configurations tied on CAL at 17 of 30, and the tie rule chose one that
   recalled 22 of 54 on EVAL, below the homogeneous PyOcto with raw picks (25) and BESTPICK (26). A positive R1-a
   would have shown a gain over a sub-optimal associator only. The thesis reports the numbers but not this reading.
8. **Who verified.** The "independent agent" of the verification is a separate AI agent (the repository commits are
   co-authored by Claude Opus 5.5, and the registration says the hash is "delivered to the user through the chat").
   The thesis does not say so (see defect D1 in Section 7).

### 5.4 Deviation E56

E56 is correctly logged, scoped and tested by sensitivity; I agree that it cannot change the registered decision
(the 40 s contain two EVAL targets, no constant was set from them, and the restricted analysis gives the same
discordance 1/0). The error is not the read itself but the sentence of the amendment that denied it while the next
clause stated the 880 s limit; the amendment containing that sentence was hashed nine seconds before the run.

Verdict on Section 5: the single run, the freeze of the code that ran and the verification are real and reproducible.
The undisclosed points 1, 2 and 6 bear on how much R1 can say; none changes its registered outcome.

---
## 6. Scores

Calibration: an undergraduate thesis (trabajo de grado) in systems engineering and mathematics in Colombia. The
reference is what a strong final-year student can be expected to deliver in one or two semesters, with a director;
not a publication and not a doctorate. For each criterion I state what 5, 7 and 9 mean at that level, then score.
The doctoral comparison is given in one line per criterion.

### 6.1 Novelty and positioning (20 %)

- 5: a competent application of known methods to a known problem, with a literature review that places it correctly.
- 7: a clearly positioned, non-trivial combination or benchmark that adds something identifiable to the literature,
  with prior work credited accurately and the contribution narrowed to what is new.
- 9: a contribution that could be submitted to a specialised journal as is, clearly new with respect to the field.

The detector is, by the thesis's own registered evidence, block orthogonal matching pursuit (Vera Rodriguez et al.
2012) with a peak search; per-station clipping is common practice in back-projection (the thesis says so); pick-free
station terms from energy peaks are a joint hypocentre-station-term inversion of a classical kind. What is new is the
way it was tested: a preregistered benchmark with equal-tolerance comparators and two held-out real datasets, and a
documented negative transfer to real data with measured causes. That has value, and the positioning is honest (the
thesis states it does not beat BOMP-PK or PhaseNet association). The literature review is short (about 30 references)
and omits the deep-learning detectors and associators a reader would expect next to PhaseNet (for example
EQTransformer, GENIE, HARPA, which docs/executive_summary.md does name) and template-matching practice beyond two
citations. Score: **6.5**.

Doctoral comparison: 4.5 (incremental over BOMP plus clipping; the central hypothesis is refuted on real data).

### 6.2 Theoretical rigor (15 %)

- 5: correct use and statement of standard results, with some derivations.
- 7: original propositions about the method, correctly stated with complete proofs and explicit assumptions.
- 9: non-trivial results (for example a calibrated null law of the deployed detector, or recovery guarantees under
  timing error) tightly connected to the experiments.

Section 4 and AA.9 meet the description of a 7 in form: correct statements, complete proofs, explicit assumptions
(Section 3.1). The results are elementary, and the connection to the experiments is overstated in three places
(Section 3.2). The evidence volume also contains a larger theory of the refuted group BLASSO (AA.1 to AA.8), which I
did not audit in detail and which does not bear on the defended results. Score: **6.5**.

Doctoral comparison: 5.0 (correct but elementary; nothing on the null law of the scan with deflation, which is what
an operating threshold needs).

### 6.3 Experimental design (20 %)

- 5: sensible experiments with at least one published baseline and some uncertainty quantification.
- 7: controlled comparisons with ablations, held-out evaluation data and statistics appropriate to the dependence
  structure.
- 9: all of the above plus preregistration, adequately powered confirmatory tests and more than one independent data
  source per claim.

Strengths: preregistration with decision rules and declared predictions; ablations that isolate the radiation operator
(2B), the peak search and the pursuit (2D); equal-tolerance comparators; Holm correction; hierarchical and block
bootstrap with the operating point re-chosen per replicate; two held-out real datasets; a simulation that replaced an
arbitrary non-inferiority margin before R1. Weaknesses: one synthetic generator family; operating points chosen on
the evaluation data for every method in every task; four of eight synthetic tasks re-used the windows of 2A after its
results were known (C2); the main registered contrasts confound the bound with the pursuit (C4); R1 had power of
about 0.15 to 0.25 for the effect seen in development and spent the last unread dataset; in R1 the location tolerance
is wider than the swarm and duplicates are never false (5.3). Score: **7.5**.

Doctoral comparison: 6.0 (the protocol is above the field's norm; the power and data-independence are not).

### 6.4 Empirical evidence (15 %)

- 5: results that support modest claims, with some uncertainty.
- 7: results that clearly support the stated claims, with intervals, on more than one dataset where the claim needs it.
- 9: strong, replicated evidence across independent datasets, with effect sizes that matter in practice.

The numbers are right (Section 2) and reproducible. The positive claim G1-R rests on one generator family and one
primary real-noise day of one network, with a smaller effect on two near-Gaussian days; its attribution to the bound
rather than the pursuit rests, in the thesis, on an exploratory seed (C4). G1-C rests on one synthetic task, mainly
one geometry. The negative real-data results are clear but come from 45 and 54 targets on two short records. The
evidence supports the narrowed claims of Section 7.3, not much more. Score: **6.5**.

Doctoral comparison: 5.0 (the only positive real-data evidence is noise from one network; the research question is
answered negatively on two small datasets).

### 6.5 Internal coherence and writing (10 %)

- 5: understandable, with inconsistencies a reader must resolve.
- 7: clear, consistent, well organised, claims matched to evidence.
- 9: polished, concise, with figures that carry the argument, and consistent across all documents.

The thesis is well organised and unusually candid; its prose is dense but precise. Against it: not a single figure
(no network map, no waveform, no FROC curve, no illustration of a ghost or of saturation) in a seismology thesis; the
protocol misstatements C1 to C3; the stale executive summary and data README (4.3); no formal front matter
(institution, programme, director) and no statement of how the work was produced (defect D1). Score: **6.0**.

Doctoral comparison: 5.0.

### 6.6 Feasibility (10 %)

- 5: the planned scope is only partly delivered, or the future work is vague.
- 7: the scope is appropriate to the level and was delivered; future work is concrete and realistic.
- 9: an ambitious scope fully delivered, with a realistic, resourced continuation.

The delivered thesis is of a scope an undergraduate can defend: one method family, a benchmark, two real tests, a
short theory. The volume behind it (189 scripts, a 700 KB evidence volume, 21 registered or confirmatory tasks in
about four weeks of log time from 23 September to 18 October 2026) is far beyond what one student produces in that
time unaided, and the commits are co-authored by an AI system; the committee cannot tell from the documents which
parts the candidate did. Future work items (onset estimates from a neural picker read at predicted times, empirical
deflation, a registered test of the logarithmic map, a second generator, external timestamping) are concrete and
feasible, but the thesis does not name the real data on which a next registered test would run, and both held-out
datasets are now spent. Score: **6.5**.

Doctoral comparison: 5.5 (the doctoral plan of Doctoral_Proposal.md was largely abandoned; what remains is feasible
but small).

### 6.7 Reproducibility (10 %)

- 5: code available with partial instructions; results reproducible only by the author.
- 7: code, data instructions with checksums, pinned environment; the main results regenerable by a third party with
  effort.
- 9: one-command or containerised reproduction, raw outputs archived, licence, no machine-specific paths.

Table 1 regenerates (as quotation); the R1 analysis and the theory check reproduce byte for byte; the R1 freeze list
matches 297 of 297 files; data sources have commits and checksums; the environment is pinned. Against: absolute paths
(/home/claude/work, /home/claude/realdata) throughout; no container or single entry point; raw outputs of 2F, 2I and
2L not archived (only analysis outputs; regeneration needs hours of CPU and the exact OpenBLAS kernel); the original
raw outputs of R0 were lost (E49) and replaced by the verifier's re-run; results/ holds English translations of what
some scripts print, so a re-run cannot be compared byte for byte for those tasks; PhaseNet code and auxiliary
packages not frozen; no licence. Score: **7.5**.

Doctoral comparison: 7.0.

### 6.8 Weighted total

| Criterion | Weight | Score (undergraduate) | Weighted | Score (doctoral, for comparison) | Weighted |
|:--|--:|--:|--:|--:|--:|
| Novelty and positioning | 0.20 | 6.5 | 1.300 | 4.5 | 0.900 |
| Theoretical rigor | 0.15 | 6.5 | 0.975 | 5.0 | 0.750 |
| Experimental design | 0.20 | 7.5 | 1.500 | 6.0 | 1.200 |
| Empirical evidence | 0.15 | 6.5 | 0.975 | 5.0 | 0.750 |
| Internal coherence and writing | 0.10 | 6.0 | 0.600 | 5.0 | 0.500 |
| Feasibility | 0.10 | 6.5 | 0.650 | 5.5 | 0.550 |
| Reproducibility | 0.10 | 7.5 | 0.750 | 7.0 | 0.700 |
| **Total (x 10)** | 1.00 | | **67.5 / 100** | | **53.5 / 100** |

Fixing defects D2 to D7 of Section 7 (text changes only) would move writing and design up by about half a point each
and the total to about 70; resolving D1 is a precondition for the grade rather than a matter of points.

---
## 7. Critical defects, by severity, with the smallest fix

| # | Severity | Defect | Smallest action that fixes it |
|:--|:--|:--|:--|
| D1 | critical (for grading, not for the science) | The thesis does not say how the work was produced. Twenty of the twenty-one commits are co-authored by an AI system (Claude Opus 5.5); the registration texts speak of "the author's decisions" and of a hash "delivered to the user through the chat"; the "independent agent" that verified R0, R1 and other tasks is a separate AI agent. A committee cannot tell which design decisions, code, analyses and text are the candidate's, and the volume of work (21 registered tasks in about four weeks) makes the question unavoidable. Most universities now require such a declaration | Add a one-page statement of contribution and tools: what the candidate decided and did (for example the four R1 decisions, the choice of hypotheses, the reading of results), what AI agents did (code, runs, drafting, verification), and describe the verifiers as separate AI sessions. Examine the candidate's command of the work in the oral defence |
| D2 | high | The first claimed contribution is overstated: "eleven registered tests ... each with a sealed registration, frozen code, a single run and an independent verification" (1.3, 7.3, README). 2B, 2C, 2D had no independent verification; 2A to 2C had no code hash list and predate the chain; ten earlier registered tasks of the same log (1j to 1y, at least one negative) are not counted | Replace the sentence by a table of all 21 registered tasks with columns: registration in chain (yes/no), code hash list (yes/no), single run, independent verification (yes/no), outcome. Count the Task 1 series or state why it is excluded |
| D3 | high | "Each task used two seeds never used before" (6.2) is false: 2B, 2C, 2D and 2E re-used the 1,200 windows of 2A after its results were known | Correct the sentence; label 2B to 2E as follow-up analyses on the 2A windows, designed after seeing 2A (and 2C for 2D, 2E) |
| D4 | medium-high | The theory is connected to the experiments more strongly than it supports: "The cap costs nothing in Gaussian noise" (contradicted by 2F standard windows, strong clutter, -0.026 and -0.032 against no cap); the bound of 12 is read as protection although it equals the threshold sigma of about 12.4; "Proved in the stated model" for the trade-off; "as Proposition 4 predicts" for a post hoc breakdown | Four sentence edits: "leaves the false-alarm tail unchanged" for "costs nothing"; report the 2F cost; state that the bound and the threshold are of the same size and that the protection shown is relative to the identity; "explained" for "proved" and "predicted" |
| D5 | medium | The deflation of the capped detector fits the moment tensor to capped vectors and subtracts a cap-sized prediction from uncapped data (c10_lib.mt_loop, c9_lib, c8_lib); large real arrivals are therefore barely removed. This is conceded in AA.5(g) but absent from Section 3.3 and from the diagnosis of R0 and of the ghosts (7.2) | One sentence in 3.3 describing the deflation as implemented, and one in 7.2 naming under-subtraction as a candidate cause for R0; optionally an exploratory check on the Ridgecrest CAL segment with the fit on uncapped vectors |
| D6 | medium | R1 limitations not stated in the thesis: the 5 km tolerance exceeds the 3.1 km extent of the truth list, so location is not tested; duplicate detections of a hit event are never false; the freeze did not cover the PhaseNet code, auxiliary packages and the PyOcto model; the quoted 16/30 and 0.265 s are development values; the expected power for the development effect was about 0.15 to 0.25; BASE was below the homogeneous PyOcto | One paragraph in 7.4 listing these, and the run values in 6.6 |
| D7 | medium | The attribution of the bounded-influence gain to the bound rather than the pursuit is argued from exploratory C18f, while the registered 2F descriptive arm BOMP-PK-CAP+SC (MT-FW minus BOMP-PK-CAP +0.004 to +0.015) is not cited; the abstract attributes all of the registered 2M contrast (+0.228, capped MT-FW against uncapped BOMP-PK) to "bounding the influence" | Cite the 2F arm; rephrase the abstract as "capped MT-FW exceeded uncapped BOMP-PK by +0.228; on synthetic data and in one exploration a capped BOMP-PK recovers most of the gap" |
| D8 | medium-low | Restated G1 (Table 2) says neither real dataset had timing errors below the tolerance; no timing residual of the pick-free readings on the Hi-net swarm is reported; and the conditional form is hard to falsify | Report the residual scale of the P1 readings on R1 CAL or limit "neither" to Ridgecrest; state the conditional hypothesis with a numerical threshold that a future registered test could refute |
| D9 | low-medium | Stale and inconsistent satellite documents: executive summary (title and body in the doctoral framing; "forty-two errors"); data/README (Hi-net "not read"); prereg README (omits R1); licence statements of the Hi-net data disagree; 2D "showed" without the Holm failure; "no deployable threshold was frozen" against rule 2A-e; "5 to 13 per cent with 9 to 14 targets" (actually 3.4 to 12.8); "+0.018 [...]" is a bootstrap median | Update or archive the executive summary as a historical document; correct the other lines |
| D10 | low-medium | No figures anywhere in the thesis or the repository; no front matter (institution, programme, director, date of defence) | Add four to six figures: the networks and events of R0 and R1; one FROC curve per real test; a waveform panel showing an S-coda ghost; the capped and uncapped score at true hypocentres (C18) |
| D11 | low | Reproduction mechanics: absolute paths; no container or single entry point; raw outputs of 2F, 2I and 2L not archived; translated analysis outputs; the Table 1 script checks substrings ("MET" inside "NOT MET"); no licence | A path configuration module; archive the part files of 2F, 2I and 2L or state their regeneration time and command; match whole tokens in the table script; choose a licence |
| D12 | low | The hash file is not self-describing (timestamp precedes its hash; first hash untimed; orphan timestamp at line 32); the translation drifts by up to a line inside paragraphs after line 777 | Two sentences in preregistration/README.md |

---

## 8. Simulated committee decision

**Decision: approved with mandatory revisions (aprobado con modificaciones), not proposed for distinction.**

Reasoning of the simulated committee (director, internal examiner, external examiner):

The science is sound and honestly reported. The candidate asked a clear question, tested it under a protocol that is
stricter than most published work in seismological method development, accepted that the answer on real data is no,
and measured why. The defended results (a registered benchmark, a modest bounded-influence effect confirmed on
synthetic data and on one real-noise day, and a quantified account of the failure on two real records) are supported
by reproducible numbers; every number I checked against its source is copied correctly, and the two re-runs I made
reproduce byte for byte. The mathematics is correct, if elementary.

Distinction is not proposed because the novelty is modest (the detector is, by the thesis's own test, block OMP with a
peak search), the positive real-data evidence is one day of one network, and the description of the protocol, which
is the first claimed contribution, is more uniform than the record (D2, D3).

Mandatory before deposit: D1 (statement of contribution and tools), D2, D3 and D4. Recommended: D5 to D10.

Questions for the oral defence, chosen to establish the candidate's own command of the work (D1):

1. Derive the bound of Proposition 4(b) on the board, and explain where the term kappa sqrt(2 n_b log L) comes from
   and why it makes the bound four times too loose.
2. Using the numbers of Task 2M, does Proposition 1 guarantee that one transient cannot produce a capped detection?
   What would you need to show it?
3. In MT-FW as implemented, what is subtracted from the data after a capped atom is accepted, and what does that
   imply for a magnitude 3 aftershock at Ridgecrest?
4. Why was rule R1-a almost certain not to be met even if the hybrid worked as in development, and what test would
   you register instead?
5. Which of the four R1 decisions were yours, and what alternatives did you consider?
6. Why are 2B to 2E not independent of 2A, and does that change the confidence you place in the conclusion of 2D?

If the candidate answers these well, the committee would grade the work in the upper part of the passing range for
an applied-mathematics undergraduate thesis; on this rubric, 67.5/100 as submitted and about 70 after the text
revisions.
