# Independent audit 7 of the undergraduate thesis "Timing-Tolerant, Bounded-Influence Moment-Tensor Detection of Overlapping Earthquakes"

Examiner: independent committee member (no part in the work). Object: repository state of 18 October 2026, commit
b4ae9b7 (github.com/juanjgc0609/timing-tolerant-mt-detection). Main document `docs/thesis/Thesis.md`; supporting
material as listed in the README. Earlier audits (`docs/audits/independent_audit_1.md` to `_6.md`,
`docs/audit_status.md`) were deliberately not read, so this judgement is independent of them. No file other than this
one was modified. All computations were light and were run either from the repository root or from a copy of
`code/` and `results/` in a scratch directory (to make sure nothing silently depends on the original working directory).

---

## 1. Hash chain and translation

**Method.** For every newline position in `preregistration/original_sealed_es/PREREGISTRO_1j_1k.txt` the SHA-256 of the
byte prefix ending there was computed, and each of the hashes in `PREREGISTRO_hashes.txt` was looked up in that set.
The English translation was compared with the Spanish file line by line (blank-line positions) and paragraph by
paragraph (multiset of numeric tokens). The git history was checked to see whether every committed version of the log
and of the hash file is a prefix of the final one.

**Results.**

| Check | Result |
|:--|:--|
| Number of hashes in the chain | 35 (first one without a time line; one time line, 5 Oct 02:36:52, without a hash, as the README declares) |
| Every hash is the digest of a newline-terminated byte prefix of the log | yes, 35 of 35 |
| Prefix lengths non-decreasing | yes, strictly increasing, from 39,694 bytes (line 431) to 121,757 bytes (line 1229) |
| Last hash equals the digest of the whole file | yes (5342c46b...) ; file ends with a newline |
| First hash falls inside the registration of Task 2D (log lines 409 to 437) | yes (prefix ends at line 431), consistent with Table 3 "partly" for 2D |
| R1 registration (log line 1065) chained at 18:49:11 UTC, amendment at 19:10:24 UTC, result at 02:09:37 UTC (18 Oct), post hoc entry at 02:36:40 UTC | yes (prefixes end at lines 1152, 1191, 1220, 1229) |
| English translation: same number of lines (1229) | yes |
| Blank lines at identical positions | yes (0 mismatches) |
| Numeric tokens per paragraph identical | 63 of 64 paragraphs; the exception is the word "sha256" in "las sumas sha256" (no number differs) |
| Numbers within a paragraph can drift by one line | yes, about 250 lines; this is declared in `preregistration/README.md` |
| All committed versions (8 commits) of the log and of the hash file are prefixes of the final files | yes |

**Interpretation.** The chain is internally consistent and the log is append-only across the git history. Two limits
should be stated more plainly than the thesis does. First, the chain is a list of digests of prefixes, not a
cryptographic chain anchored outside the author's control: whoever holds the file can rewrite it and recompute every
digest, so the chain by itself proves only that the stored hashes and the stored log agree. The claim of Section 5.1
that "the hash chain proves that entries were not edited after they were chained" holds only to the extent that each
hash was recorded somewhere external at the time. Second, the only external records in the repository are the git
commits, and they postdate the runs: the first commit (c5f45bd, 11 October 2026) comes after Tasks 2A to R0; the 2M
registration first appears in commit ca4185e (17 Oct 00:02 UTC), after the 2M result had been logged (16 Oct 19:59 UTC);
the R1 registration and amendment first appear in commit 0374b2f (19:11:45 UTC), 72 s after the single run started
(19:10:33 UTC). Commit timestamps are moreover set by the author's machine. The thesis says that the hashes were sent
to the author through the chat before each run, but that chat record is not in the repository, so a third party cannot
check the order of registration and run for any task. The thesis declares the self-timestamping and the R1 push delay;
it does not say that no task has a third-party-verifiable pre-run anchor. See defect D4.

## 2. Number checks

### 2.1 Table 1

`python3 code/thesis_evidence_table.py` run from the repository root exits with status 0. Its output is byte-identical
to Table 1 of the thesis (diff of the extracted Markdown table) and to `docs/thesis/evidence_table.md`. I read the
source line of each of the 16 rows to check that the claim text describes what the line contains, not only that the
values occur in it (the script checks token occurrence only). All 16 rows are correctly described, with two caveats:

- E4, E5, E6, E8, E9 and E10 are measured under the hit rule LOC (epicentre within 5 km, origin within 3 s), not under
  LOC+T; the table says "origin time required" for E7, E11 and E13 to E15, and says nothing for the bounded-influence
  rows. A reader of Table 1 will assume one hit rule throughout (defect D1).
- E1 lists point estimates. The registered rule 2A-a is a point rule, but several cell intervals contain zero (for
  example the four-station cell of seed 2024 at 0.1 false events per window: +0.14 [-0.01, +0.34]).

### 2.2 Other numbers of the thesis

Thirty-six further numbers or groups of numbers were traced to their result files. "OK" means the value in the
thesis equals the value in the file (to the printed precision).

| # | Thesis statement (section) | Source | Status |
|:--|:--|:--|:--|
| 1 | 2A: best comparator in every cell is QuakeMigrate (QM-3C+SC); 2A-e met per cell, not in the strict reading (6.2) | `results/task2a/analyze_2a_out.txt` tables 2 and 4, rules lines 61-66 | OK |
| 2 | 2B: loses in 2 of 3 cells (weak clutter) and 3 of 3 (strong) (6.2) | `results/task2b/analyze_2b_out.txt:53,104` | OK |
| 3 | 2C-a: +0.077 to +0.519 (6.2) | `results/task2c/analyze_2c_out.txt:44,83` | OK |
| 4 | 2D-b: +0.064 to +0.412; after Holm only the one-sided 8-station cell of one seed per regime (6.2) | `results/task2d/analyze_2d_out.txt:85,148,210,271` | OK |
| 5 | 2E: +0.000 to +0.039, every interval contains zero (6.2) | `results/task2e/analyze_2e_out.txt:67,160` | OK |
| 6 | 2F: +0.043 to +0.074, Holm rejects (6.2) | `results/task2f/analyze_2f_out.txt:26-31` | OK |
| 7 | 2F standard windows, strong clutter: 0.379 and 0.351 capped vs 0.405 and 0.383 uncapped (4) | `analyze_2f_out.txt`, std strong rows | OK |
| 8 | 2F descriptive arm: BOMP-PK-CAP 0.243 and 0.245 vs MT-FW 0.251 and 0.252, "within 0.007 and 0.008" (6.3) | `analyze_2f_out.txt`, dense both rows | Values OK; the two gaps are 0.008 (4242) and 0.007 (8888), order swapped |
| 9 | 2F with origin time: +0.024 to +0.045 (6.2) | `results/task2f/verif_2f/task2f_verification.txt:70` | OK |
| 10 | 2I: +0.058 and +0.053; non-inferior to the cap within 0.03 (6.2) | `results/task2i/analyze_2i_primary_out.txt:82-89` | OK |
| 11 | 2I: Huber costs about nine times as much (6.2) | Appendix Y.26 (evidence volume line 4879) | OK |
| 12 | 2L: +0.057, +0.062 standard; +0.038, +0.030 dense (6.2) | `results/task2l/analyze_2l_out.txt:113-121` | OK |
| 13 | 2L: S2 terms worse than S0 in about 22 per cent of networks (6.2) | `results/task2l/verif_2l/task2l_verification.txt:122` | OK |
| 14 | 2M secondary days: +0.057 and +0.120 (6.3) | `results/task2M/analyze_2m_out.txt:344,346` | OK |
| 15 | 2M breakdown: +0.033 at night, +0.300 between 7.5 and 15 h; Huber +0.10 to +0.15 in quiet noise (6.3) | `results/task2M/posthoc/blocks_out.txt` | OK |
| 16 | 2M capped threshold 77.4, so sigma about 12.4 (4) | `analyze_2m_out.txt` (TRANSFER line of seed 2718); sqrt(2 x 77.4) = 12.44 | OK |
| 17 | Transient of 12,500 noise deviations, "the largest found" (4) | `results/task2M/log_entries/ver_2M_en.txt` | Value OK; the verifier reports "a transient of ~12500", not that it was the largest |
| 18 | R0: at most 0.022 vs up to 0.667 (6.4) | `results/taskR0/analysis_eval.txt:5` | OK |
| 19 | R0 closure: +0.000 [-0.023, +0.021] and -0.050 [-0.101, -0.013] (6.4) | `results/taskR0/run_close/analysis_eval.txt:13,16` | OK |
| 20 | C18: AUC capped 0.39 to 0.49, uncapped 0.69 to 0.74 (6.5) | `results/explorations/c18/c18b_{N4,N6,N8,U8}.json` (REG 0.394-0.494, NOCAP 0.687-0.740) | OK |
| 21 | C18i: layered 0.45 against homogeneous 0.37 to 0.39 within 0.25 s (6.5) | `results/explorations/c18/c18i_timing_s2_{N8,U8}.json`, key `models` | OK as stated, but see defect D8: with the event term (key `_evterm`) the values are 0.427 vs 0.427 (N8) and 0.417 vs 0.387 (U8) |
| 22 | C18j: 5/16, 5/16, 1/16 (E12) | `results/code_audit/c18j_loeo_bh.txt:8` | OK |
| 23 | R1: HYB 23, BASE 22, BESTPICK 26, P1 1 of 54; p = 0.50 and 0.97 (6.6) | `results/taskR1/run/analysis_eval.txt` | OK |
| 24 | R1 block bootstrap +0.018 [+0.000, +0.058] and -0.052 [-0.163, +0.000] (6.6) | same file | OK |
| 25 | R1 JMA-only: 9 of 10 (HYB, BASE), 10 of 10 (BESTPICK), 1 of 10 (P1) (6.6) | same file | OK |
| 26 | R1 deployable threshold: 41.8 and 12.9 false per hour (6.6) | same file | OK |
| 27 | R1 homogeneous PyOcto with raw picks: 25 of 54 (6.6) | same file (R0-pyocto, LOC+T 30/h) | OK |
| 28 | R1 development: BASE 16, HYB 16 (17 with registered terms), P1 2 of 30 (6.6) | `results/taskR1/dev/analysis_cal_dev.txt`, `_v1.txt` | OK |
| 29 | R1 run lags 0.272 / 0.244 s, CAL 17 of 30; amendment values 0.265 / 0.244 s (6.6) | `results/taskR1/verif_r1/r1_verification.txt` sections 4.5, 6.6 | OK |
| 30 | R1 restricted to origins from 880 s: 22 vs 21 of 52; 29 of 54 targets found by any method (6.6) | `r1_verification.txt` sections 5, 7 | OK |
| 31 | R1 post hoc timing: P 0.22 s, 55 % within 0.25 s, 36 % beyond 0.5 s; S 0.17 s, 76 %, 20 % (7.2) | `results/taskR1/posthoc/timing_cal_out.txt` | OK |
| 32 | Rule simulation: bootstrap NI false positives 3.4 to 12.8 %; McNemar 0.000 to 0.018 at zero; power 0.13-0.17 and 0.49-0.69 (40 targets), 0.35-0.53 and 0.92-0.99 (100 targets) (5.3) | `results/taskR1/sim_r1_rules_out.txt` | OK |
| 33 | Theory: kappa = 2.146; 9.4e-6 per hypothesis; 3.4e-4 per window; quantiles 7.18 / 7.18 / 6.78 vs bound 28.50 (4) | `results/theory/tail_check_out.txt`; recomputed: Pr(chi2_3 > 36) = 7.49e-8, 126 x = 9.44e-6, 4480 x = 3.35e-4, bound 28.49 | OK |
| 34 | Proposition 1 numbers: 12 (cap), about 47 (log), 12,500 (identity) (4) | recomputed: 2 x 6 x sqrt(log(1 + (12500/6)^2)) = 46.9 | OK as an upper bound; see defect D11 |
| 35 | E52: exact-time readers keep 0.78 instead of 0.94 of the energy (6.2) | `docs/artifact_manifest_and_error_log.md` E52 | OK |
| 36 | Fifty-six errors E1 to E56 (1.3, 9) | error register, headings E1 to E56 all present | OK |

No number of the thesis was found to differ from its source. The problems found are of selection and labelling
(rows 8, 17, 21 and the hit rule of Table 1), not of transcription.

### 2.3 Table 3 against the log and the result folders

| Task | Table 3 claim | What the log and folders show | Verdict |
|:--|:--|:--|:--|
| 2A | not chained, no hash list, single run, verified | registration at log line 251 (before the chain); `results/task2a/verif_2a/` and `task2a_report.txt`; no code hash file | OK |
| 2B | not chained, no hash list, single run after a stopped first attempt, not verified | `log2b_attempt1.txt`, `task2b_control_failed_attempt1.txt`; no verification folder | OK |
| 2C | not chained, no hash list, single, not verified | no hash file, no verification folder | OK |
| 2D | partly chained, hash list, single, not verified | first hash inside the registration (2.1); `code/task2d_code_hashes.txt` | OK, but the current `code/analyze_2d.py` does not match its frozen hash (110aa893 vs 8dd5013b); the change is declared in the artifact manifest and in `results/task2d/README.txt` (one diagnostic line) and is not mentioned in Table 3 |
| 2E | chained, hash sent, hash list, single, verified | log line 449; `code/task2e_code_hashes.txt` matches the repository (5 of 5); `results/task2e/verif_2e/` | OK |
| 2F | same | hash list matches (26 of 26); `verif_2f/`; log line 552 | OK |
| 2I | same, resumed after two interruptions | hash list matches (27 of 27); log lines 590-604 | OK |
| 2L | same | hash list matches (37 of 37); `verif_2l/`; log line 664 | OK |
| R0 | same, resumed after an interruption | `code/r0/R0_FREEZE.txt` matches for every file present (27; 7 external); `results/taskR0/README.txt` lines 10-11 | OK |
| 2M | same, "one lane resumed after a declared failure (E48)" | `T2M_FREEZE.sha256` matches for every file present (26; 19 external); besides the REAL lane of seed 3141 (E48), the GEN arms were stopped and resumed after a container restart (`run_*_resume.log`, README lines 11-13) | Understated: three resumptions, the extra ones on descriptive arms only |
| R1 | same, single, verified | `code/r1/FREEZE_R1.sha256` has SHA-256 1669903f...28a2a3 as quoted; 291 of its 297 file hashes match the repository copies and the data files, the other 6 are PhaseNet checkpoint files outside the repository; run 19:10:33 to 19:20:56 UTC without interruption (`run/queue.log`, verification 1.2) | OK |

"Hash sent before the run" cannot be checked from the repository for any task (Section 1).

### 2.4 Reruns

| Item | Command | Result |
|:--|:--|:--|
| Table 1 | `python3 code/thesis_evidence_table.py` (repository root) | exit 0, identical |
| R1 analysis | `python3 -I code/r1/r1_analyze.py eval results/taskR1/run/hyb_eval.json results/taskR1/run/hyb_cal.json OUT.txt` from a scratch copy of `code/` and `results/`, with the public Hi-net catalogues in place | 2.6 s, output byte-identical to `results/taskR1/run/analysis_eval.txt` (also under Python 3.13 and NumPy 2.5, not only the pinned versions) |
| Numerical check of Proposition 4 | `code/theory_tail_check.py` after editing one hard-coded path | 39 s, output byte-identical to `results/theory/tail_check_out.txt` |
| Figures 1 to 3 | `code/thesis_figures.py results/taskR1/run OUTDIR` after editing one hard-coded path | Figures 1 and 2 build with the system Python; Figure 3 fails without `pyproj` (imported through `qm_wrapper.py`), which is not in `requirements-frozen.txt`; in the project virtual environment (which has pyproj 3.7.2) all three PNG files are byte-identical to `docs/thesis/figures/` |

## 3. Mathematics (Section 4 and Appendix AA.9)

Notation as in the thesis: blocks v_k in R^3, k = 1..n_b, n_b = 2 n_s; Q_h with orthonormal (or zero) columns;
phi radial, phi(x) = r(|x|) x/|x|; sigma_phi = |Q_h^T phi(v)|.

**Lemma AA.9.1 (contraction).** Correct and complete: Q_h Q_h^T is an orthogonal projector, so |Q_h^T a| <= |a|, and
the reverse triangle inequality gives the Lipschitz statement.

**Proposition AA.9.2 (influence of one block).** Correct. The chain |sigma(v) - sigma(v')| <= |phi(v_k) - phi(v'_k)|
<= r(|v_k|) + r(|v'_k|) is valid; (a) and (b) follow from r_c <= c and the monotonicity of l_c; the proof of (c)
(sigma(v) - sigma(v') >= gamma |e| - 2 |Q^T v'|) is correct. Two remarks on the numbers attached to it. (i) The quoted
"at most 12 with the cap, about 47 with the logarithmic map" use the bound with both blocks at the large amplitude M.
For the situation the text describes, one clean block replaced by a transient, the same proposition gives
r(M) + r(|v_k|), that is c + |v_k| (about 8 for a clean block of norm about 2) for the cap and about 23.5 + 2 for the
logarithmic map. The sentence "one transient can move a capped score by about the threshold itself" (12 against 12.4)
is therefore pessimistic by about a factor of 1.5; the qualitative conclusion stands. (ii) The bound is a one-line
consequence of the triangle inequality; it is correctly presented as elementary.

**Proposition AA.9.3 (saturation).** Correct and trivial: |r_c(v)|^2 = sum min(|v_k|, c)^2 <= n_b c^2, and if every
block is saturated r_c(v) = c u. The "consequently" clause is a restatement. The thesis adds an interpretation that
the proposition does not prove: that the amplitude pattern across stations "is what distinguishes the radiation of one
real source from a coincidence". The proposition shows that amplitudes are discarded when all blocks saturate; it does
not show that the discarded information is what matters, which is the role of the C18 measurement (AUC 0.39 to 0.49
capped against 0.69 to 0.74 uncapped). The thesis makes this link explicitly, so the strength is acceptable. One
missing observation: for any radial map other than the identity, phi(G m) is in general not in the column space of G
unless all block norms are equal, so even a noise-free single event is scored below its full energy by the cap and the
logarithmic map. This is the precise reason for the cost in isolated windows that Section 4 reports from Task 2F
(0.379 against 0.405), and it would be worth one sentence.

**Proposition AA.9.4 (logarithmic map).** Correct under its hypothesis (equal block norms). The hypothesis is very
restrictive (it excludes exactly the cross-station amplitude pattern of Proposition AA.9.3), and the appendix says so.
The thesis title of the proposition, "keeps the amplitude ordering", should read "keeps the ordering in event size
for events with equal block norms"; nothing is proved about the cross-station pattern under the log map, which the map
also distorts.

**Proposition AA.9.5 (null tail and contamination).** I checked each step.

- (a) Correct: if every candidate has norm at most c the cap is the identity on every block; a union bound over the
  n_b L candidates, each |x|^2 ~ chi2_3, gives the probability. The consequence drawn in the thesis, "the cap leaves the
  Gaussian false-alarm tail unchanged", should be stated as what (a) gives: the two scores differ on an event of
  probability at most 9.4e-6 per hypothesis (3.4e-4 per 40 s window), so any tail probability changes by at most that
  amount. Since the operating budgets of the tasks are 0.1 to 0.25 false events per window, this is enough for the
  practical claim, but "unchanged" is not literally true for tail probabilities below 1e-5 (and the event "some block
  exceeds c" is strongly positively correlated with a large score, so in the far tail the two scores do differ).
- (b) Correct under the stated assumptions. For a fixed choice s of one candidate per block, f_s is the composition of
  a linear map of the white noise with operator norm at most kappa (convolution with k/|k|_2, then sample selection),
  a block-wise radial map that is 1-Lipschitz (the constant of a radial map is max(sup r', sup r(t)/t); for l_c,
  l_c'(t) <= (1 + t^2/c^2)^(-1/2) follows from log(1+u) >= u/(1+u), which I verified), and the contraction. The
  maximum F of the L^(n_b) functions f_s is kappa-Lipschitz and dominates the read score. The Gaussian concentration
  inequality, the bound E f_s <= sqrt(E f_s^2) <= sqrt(tr(Q^T C Q)) <= sqrt(6) with C block-diagonal and C_kk <= I_3,
  and the maximal inequality E max_s (f_s - E f_s) <= kappa sqrt(2 n_b log L) are all standard and correctly combined.
  Two assumptions are used but not stated: that the three components of a block are independent (the proof calls the
  block "isotropic"), and that the normalisation is exact (the deployed detector divides by a robust estimate computed
  from the data, which is not a Lipschitz map of the noise). Neither is needed in full: E[phi phi^T] <= E[x x^T]
  holds for any centred Gaussian block because r(t)/t <= 1, so only "unit variance per component" is required for the
  sqrt(6) term. The numerical check (random 6-dimensional Q in R^42, independent blocks 40 samples apart, 2 x 10^5
  replicates) confirms the bound with a factor of about four of slack.
- (c) Correct: m changed blocks move phi(v) by at most 2c sqrt(m) in norm, and the coupling argument gives the
  threshold statement for one hypothesis.

**Connection to the experiments.** The thesis states the right strength in most places: Section 4 says the
propositions were written after the experiments, describe a trade-off, and do not cover deflation, the grid, station
terms or non-Gaussian coloured noise; Table 2 says "consistent with", not "proved". Two places overstate it.
Section 8 calls the bounded-influence result "explained by theory", and Section 7.3 "explained by a bound that does
not depend on the amplitude of the transient". What the propositions explain is why an unbounded score can be moved
arbitrarily by one block and a capped one cannot; that is the definition of a cap restated as an inequality. They do
not predict the size of the gains (+0.04 to +0.07, +0.23), do not model the operating-point mechanism that actually
produced the 2M result (the uncapped BOMP-PK threshold rose to 438-2228 score units on the transient day and its recall
fell to 0.003 in the transient hours), and do not cover HUBER-N, which gave the largest 2M gain. "Consistent with"
would be the correct verb throughout.

**Overall judgement on the mathematics.** For an undergraduate thesis in applied mathematics the section is correct,
honestly scoped and numerically checked. It is not deep: Lemma 1 and Propositions 1 to 3 are one-line consequences of
the definitions, and Proposition 4 is a careful but standard application of Gaussian concentration and a union bound.
No proof is wrong. The defects are of statement precision (D10, D11, D12) and of the verb used to connect them to the
data.
I did not audit AA.1 to AA.8 (theory of the refuted group BLASSO detector), which the thesis does not rely on.

## 4. Claims against evidence

### 4.1 What is reported well

The negative results are reported as negative and with the same prominence as the positive ones: the refutation of
G1-M (2E), the failure of 2D-a (MT-FW equals BOMP-PK with the same peak search), the R0 result (-0.383), the 1 of 54 of
the pick-free detector in R1 and the null R1-a are all in the abstract. Section 1.3 states what the thesis does not
claim. Table 2 restates every hypothesis after the real-data failures. The error register is real and specific (E52,
E53, E56 change the reading of earlier explorations and are cited in the text). The development leakage of R1 (E56)
is in the thesis, the evidence volume and the log. This is unusually candid for a thesis at any level.

### 4.2 Overclaims and imprecise claims

1. **Hit rule of the bounded-influence results (D1).** The headline numbers for G1-R (+0.043 to +0.074 in 2F, +0.058
   and +0.053 in 2I, +0.228 and +0.241 in 2M, and the C18f +0.026) are under LOC. The thesis gives the origin-time value
   for 2F (+0.024 to +0.045) but not for 2I (+0.031 and +0.019, from the 2I verification entry in the log) nor for 2M
   (+0.085 [+0.053, +0.120] and +0.104 [+0.080, +0.134], `analyze_2m_out.txt` 2M-c lines). Meanwhile every real-data
   result (R0, R1) is under LOC+T. The abstract compares "+0.228" on 2M with "0.022 against 0.667" on R0 without saying
   that the two are measured under different hit rules. The 2M registration does make LOC primary, so this is not a
   deviation, but it is a material omission in the reporting.

2. **"Confirmed on real noise" (D2).** The 2M primary result is one day of noise (2010-12-31, all 72 parts), with
   synthetic events injected by the same generator family. The registered contrast is a capped statistic against an
   uncapped one under a common operating rule; in the hours with transients the uncapped BOMP-PK recalled 0.003 because
   its operating threshold was pushed above the transients. The 2M verifier wrote that this "is the mechanism of the
   hypothesis and not a general win over BOMP-type detectors", and reported dependence-aware intervals of [+0.132,
   +0.314] (hours clustered) and [+0.088, +0.314] (moving block) next to the registered [+0.148, +0.298]. The thesis
   keeps the narrower interval in the abstract and Table 1, and calls the noise "real noise of independent days" in
   Section 1.3 although the confirmatory evidence is one day (the two near-Gaussian days give +0.057 [+0.004, +0.121],
   with origin time +0.037 [-0.003, +0.076]). The thesis also states that clipping each station before stacking is
   common practice in back-projection, without a citation, and no registered task compares the cap with such standard
   preprocessing (despiking, gap and glitch removal, per-station clipping before the projection). Exploration C18f
   (capped BOMP-PK recovers the gap) is the closest, and it shows exactly that the gain is a generic robustification.
   The honest form of the contribution is "a measured size, inside a moment-tensor subspace score, of a known
   robustification", which Section 7.3 nearly says; Section 8 ("what survives as a contribution is the bound on the
   influence of each station, confirmed on real noise and explained by theory") overstates it.

3. **Restated G1 (D3).** Table 2 "retains" G1 "as a conditional, refutable statement: it can hold where at least 90 per
   cent of the energy-peak readings of calibration events, with held-out station terms, lie within the reading
   tolerance". The threshold of 90 per cent was introduced in a post hoc log entry (18 Oct 02:36 UTC) after both real
   datasets were known to fail it, and it is not shown anywhere that the synthetic benchmark, where G1 was "supported",
   meets it when measured the same way (the thesis quotes a timing error of "about 0.14 s" for the generator, a
   different quantity). The condition is also confounded: on synthetic data the comparators were fed by an STA/LTA
   picker (and with optimistic picks GaMMA was ahead in four of six weak-clutter cells at 0.1 false events per window),
   while on real data they were fed by PhaseNet. The R0 and R1 failures are therefore explained at least as well by the
   change of comparator as by the timing. As written, the restatement protects the hypothesis rather than exposing it.
   It should either be measured on the synthetic windows with the same procedure and stated with the comparator
   condition, or labelled as a conjecture.

4. **Conclusion on synthetic data (D7).** Section 8 says the detector "beats pick-based association on a synthetic
   benchmark where the timing error is smaller than the tolerance". In 2A the best comparator in all twelve cells was
   QuakeMigrate, a migration method, not a pick-based associator; the associators were starved by an STA/LTA picker;
   operating points were chosen on the evaluation windows; and 2D showed that the advantage belongs to timing-tolerant
   block pursuit in general (BOMP-PK equals MT-FW). Section 6.2 and Table 2 contain all of these qualifiers; the
   conclusion, which is what most readers retain, drops them.

5. **Protocol uniformity (D5, D6).** The README's status paragraph says "Eleven tests were registered before their data
   were analysed, run once with frozen code and checked by an independent verifier". Table 3 of the thesis shows that
   2A to 2C had no code hash list, 2B to 2D had no verification, and 2A to 2C predate the chain. Section 5.1 of the
   thesis opens with "Every confirmatory test followed the same steps" and lists freezing and verification among them,
   before Table 3 corrects it. The abstract's "Every confirmatory claim was registered before its data were analysed" is
   also too strong for 2B to 2E, which were registered after 2A (and 2C) had analysed the same windows with other
   methods; Section 6.2 says so, the abstract does not.

6. **Medium at Ridgecrest (D8).** Section 6.5 and the lesson of Section 7.5 ("the first conclusion that the velocity
   model does not matter ... was reversed when the window was widened") rest on the fractions 0.45 against 0.37 to 0.39.
   These are the readings without the event term. The evidence volume itself argues (Y.32, C18a) that the relevant
   measure for a detector is the one with the event term removed, because the detector estimates the origin time. With
   the event term, in the same files, the layered and homogeneous media give 0.427 and 0.427 on N8 and 0.417 and 0.387 on
   U8: the "reversal" is mostly an origin-time shift, and the medium still matters little. The stronger reading was
   selected for the thesis.

7. **Smaller items.** Section 3.2 states that 2D showed the peak search "carries most of the advantage over exact-time
   sparse inversion" without the E52 qualifier that Section 6.2 adds (the exact-time readers were handicapped).
   Section 6.6 and the R1 registration call the velocity model a "28-point approximation" while the 13 points listed in
   the log and in `code/explore_contrib/c18_tt.py` are the ones used (D9). Section 4 calls the 12,500-deviation
   transient "the largest found in the real noise of Task 2M"; the verification says only that one such transient
   passed the quality control. Section 6.2 says the associators and QuakeMigrate were "all fed by an STA/LTA picker";
   QuakeMigrate uses STA/LTA onset functions, not picks.

### 4.3 Consistency across documents

| Statement | Thesis | README | Evidence volume (Y) | Executive summary | Verdict |
|:--|:--|:--|:--|:--|:--|
| R1 numbers (23/22/26/1, p 0.50/0.97) | 6.6 | table | Y.33 | 17 Oct status | consistent |
| R0 numbers and closure | 6.4 | table | Y.30 | 10 and 16-17 Oct status | consistent |
| 2M numbers | 6.3 | status paragraph | Y.31 | 16-17 Oct status | consistent |
| Protocol steps per task | Table 3 | "run once with frozen code and checked by an independent verifier" for all eleven | per task | not stated | README contradicts Table 3 |
| Model of R1 | "28-point" | not stated | "28-point" | not stated | contradicts the 13 points of the log and code |
| Status of G1 | conditional on 90 % | not restated | not restated | "R1 gives no evidence for G1" | thesis-only post hoc restatement |
| Scope of the evidence volume | n/a | n/a | header still says "the eight preregistered synthetic tasks ... and three [explorations] that followed the eighth" | n/a | stale header; Y.30 to Y.33 exist |
| Doctoral proposal | "kept as the record of the original plan" | same | n/a | n/a | banner present; its subtitle "all evidence to date is synthetic" is stale but clearly labelled historical |

The executive summary is consistently framed as a dated decision record and its latest paragraph agrees with the
thesis. Appendix Y.33 agrees with the thesis in every number I compared.

### 4.4 Positioning in the literature

Section 2 covers pick-based association, subspace detection, sparse moment-tensor inversion and robust statistics,
but omits the closest waveform-migration detectors and the hybrids that already exist: the source-scanning algorithm
(Kao and Shan 2004), LOKI (Grigoli et al. 2013), BackTrackBB (Poiata et al. 2016), and MALMI (Shi et al. 2022), which
migrates neural-picker probabilities. Future-work item (1) of Section 8 ("reading the pick probability of a neural
picker at the predicted time instead of the matched-filter energy") is essentially what MALMI and probability-based
migration already do; it should be positioned against them. EQTransformer is not mentioned among the pickers. The
reference list inherits the numbering of the proposal's bibliography and mixes numbered and author-year citations,
which is acceptable for a draft and not for a deposited thesis.

## 5. Protocol of Task R1

**Sequence, as recorded.**

| Time (UTC, 17-18 Oct 2026) | Event | Evidence |
|:--|:--|:--|
| 17 Oct, 00:02 | Hi-net checksums committed (data downloaded, not analysed) | commit 7fc3789, `data/checksums/hinet_matchlocate2_example.sha256` |
| 17 Oct, before 18:49 | Rule simulation (`code/hyb/sim_r1_rules.py`) replaces the -0.20 non-inferiority margin by an exact McNemar superiority test | `results/taskR1/sim_r1_rules_out.txt`; numbers checked in 2.2 row 32 |
| 18:49:11 | Registration chained (hash 3b238aa4...) | log line 1065; chain prefix ends at line 1152 |
| 18:49 to 19:10 | Development on CAL (PhaseNet, associator grids, P1, hybrid) | `results/taskR1/dev/` (two analyses: registered terms, amended terms) |
| 19:09:24 | Freeze list written; latest frozen file 19:09:14 | verification 1.2 |
| 19:10:24 | Amendment 1 chained (hash d5d420a6...), quoting the freeze hash 1669903f... | log line 1154 |
| 19:10:33 | Freeze file rewritten with identical content 0.03 s before the run; single run starts | verification 1.2 and D5 |
| 19:11:45 | Registration and amendment first committed (commit 0374b2f) | git |
| 19:20:56 | Run ends without interruption | `results/taskR1/run/queue.log` |
| 18 Oct, 02:09:37 | Result and independent verification chained | log line 1193 |
| 02:36:40 | Post hoc timing analysis chained, declared as such | log line 1222 |

**What is sound.** The decision rule was validated by simulation before registration, and the validation changed it
for the better (a margin with no justification and an anticonservative bootstrap were replaced by a test with known
size). The amendment was made before the freeze, is described, was chosen by an alignment criterion, and in
development it made the hybrid slightly worse (17 to 16 of 30), which argues against recall-driven tuning. The freeze
list covers every module the run could import and the data; all 291 hashes that refer to files in the repository or the
public data match. The run was single and uninterrupted. The verification re-derived every number with independent
code, including the 277 hybrid detections and their pick support, and it found E56, which was then logged and reported.

**Deviations and weaknesses, and whether they are disclosed.**

| Item | Disclosed in the thesis? | Comment |
|:--|:--|:--|
| Development read EVAL samples [840, 880) s; the amendment states the opposite (E56) | yes (6.6, 7.4, E56) | No constant set on them; the sensitivity (22 vs 21 of 52) is reported. Correctly handled. |
| Amendment quotes development values (lag 0.265 s, pick terms, CAL 16/30) that differ from the run (0.272 s, 17/30) | yes (6.6) | Cosmetic. |
| Freeze does not cover the PyOcto model file, the PhaseNet source code and two auxiliary package folders | yes (7.4) | The model file is not in the repository either; a third party has to rebuild it. |
| Freeze file rewritten 0.03 s before the run (same content) | no (only in the verification, D5) | Harmless; one sentence would close it. |
| Registration public only after the run started; the chat record with the hashes is not in the repository | partly (7.4 says "pushed about one minute after the run started, while their hashes had been sent before it") | The claim that the hashes were sent before the run is not checkable from the repository (Section 1). |
| PhaseNet run once on the full record, so CAL picks (and the CAL choice of BASE) were produced with EVAL samples in the normalisation context of the last CAL chunk | registered in the amendment (item 3); not discussed in the thesis | Minor; it is the registered procedure, but it is a second, smaller route by which EVAL data touched a CAL decision. |
| R1-b compares with a BESTPICK selected on EVAL, while BASE was selected on CAL | yes ("optimistic for the comparators") | Correctly framed. |
| Hybrid design selected on 16 Ridgecrest events | yes (registration limits; 7.4) | |
| Power: about 0.15 to 0.25 for the gain seen in development; the registration says a negative R1-a "does not by itself refute the hybrid" | yes (5.3, 7.4) | A design weakness, not a protocol deviation: the last unseen dataset was spent on a test that could not distinguish "no effect" from "the effect expected". The thesis reports "no gain demonstrated", which is the right phrase. |
| Truth list from template matching with the JMA events; 29 of 54 targets attainable; two M >= 2 targets possibly artefacts | yes (6.6, 7.4) | |
| Post hoc timing entry, "requested by the sixth audit", used to restate G1 with a 90 % threshold | yes as post hoc; the threshold itself is not justified | See D3. |

**Verdict on R1.** Nothing material is undisclosed. The two points that are not in the thesis (the freeze file rewrite
and the PhaseNet normalisation context) cannot have changed a registered number. The real limitation of R1 is its
design (power, truth list, 33 minutes, seven stations), and the thesis states it.

## 6. Scores

### 6.1 Calibration and scores

The anchors below are for an undergraduate thesis (trabajo de grado) in systems engineering and mathematics in
Colombia. A 5 is a passing thesis with clear weaknesses in that criterion, a 7 is a good thesis, a 9 is exceptional
work that would stand out in any undergraduate programme.

**Novelty and positioning (20 %).**
5: a known method applied to a new case, with an adequate literature review. 7: a clearly positioned adaptation or
combination with a measured contribution beyond replication and a literature review that covers the direct
competitors. 9: an original method or finding of publishable interest, positioned against the full state of the art.
*Assessment.* The pick-free detector is a modest variant of block OMP over moment-tensor blocks, and the registered
tasks show that its own elements either add nothing over BOMP with the same tolerance (2D) or fail on real data (R0,
R1). The surviving positive result (a per-station bound) is, by the thesis' own account, common practice in
back-projection, and its measured gain is against an uncapped comparator that collapses with transients. What is
genuinely uncommon is the preregistered, falsification-oriented benchmark and the measured account of why a pick-free
detector does not transfer to real sequences. The literature review misses the closest migration detectors and
the existing probability-migration hybrids (4.4). **Score 6.0.**

**Theoretical rigor (15 %).**
5: standard theory used correctly, derivations mostly correct with gaps. 7: correct and complete proofs of non-trivial
statements under stated assumptions, with numerical checks. 9: original, non-trivial theory, fully rigorous.
*Assessment.* All proofs of AA.9 are correct; Proposition AA.9.5 is a careful application of Gaussian concentration
with a numerical check; assumptions are mostly stated and the scope is honest. The results are elementary and partly
tautological (Lemma 1, Propositions 1 to 3), two assumptions are implicit, and the connection to the data is
occasionally phrased as "explained" rather than "consistent with". **Score 7.0.**

**Experimental design (20 %).**
5: a sensible experiment with a baseline and some uncontrolled confounds. 7: controlled comparisons with appropriate
baselines, ablations, held-out data and uncertainty quantification. 9: preregistered, adequately powered, replicated,
with strong baselines and more than one real dataset.
*Assessment.* Preregistration, frozen code, single runs, two fresh seeds per synthetic task, ablations (2B, 2D, the
2F mechanism arm), hierarchical and block bootstraps, Holm correction, and a decision rule validated by simulation are
well above the undergraduate norm. Against that: the synthetic generator is matched to the detector (double-couple
sources inside its subspace, homogeneous medium, Gaussian noise); the synthetic comparators were fed by a weak picker;
operating points were chosen on evaluation data for every method; 2B to 2E reuse the 2A windows; the real tests are
one hour, 33 minutes and one day of noise, and R1 was knowingly underpowered and spent the last unseen data.
**Score 8.0.**

**Empirical evidence (15 %).**
5: results reported that partly support the claims. 7: results that support the stated (suitably limited) claims, with
uncertainty, negative results included. 9: strong, replicated evidence on real data for the main claims.
*Assessment.* The main hypothesis is refuted on real data and the thesis says so; that is not penalised. The evidence
for what survives is real but thin: G1-R rests on synthetic data and one day of real noise with synthetic events, under
the lenient hit rule, against a comparator that fails for a known reason; G1-C on one synthetic test driven by one
network geometry; the hybrid on 54 targets with one discordant event. Every number checked is right (Section 2).
**Score 6.5.**

**Internal coherence and writing (10 %).**
5: readable, with contradictions between sections or documents. 7: clear and consistent, claims match evidence with
few slips. 9: polished and fully consistent.
*Assessment.* The thesis is clear, concise for its scope, and honest; Tables 1 to 3 do real work. Slips: the README
status paragraph contradicts Table 3; Section 5.1 and the abstract overstate protocol uniformity; the conclusion drops
qualifiers that Section 6 states; the hit rule is unlabelled in Table 1; the "28-point" model; tables numbered out of
order (Table 3 precedes Table 2); mixed citation styles and a bibliography numbered after another document.
**Score 7.0.**

**Feasibility (10 %).**
5: scope too large or only partly executed. 7: scope executed with the available resources; the candidate can be
expected to command it, with some risk. 9: well scoped, fully executed, candidate clearly in command.
*Assessment.* Everything planned for the thesis was executed, with modest computing resources, and the scope was cut
honestly when the real-data results came in. The size of the record (eleven registered tasks, an evidence volume of
about 700 kB, fifty-six logged errors) is far beyond an undergraduate thesis and was reached with an AI system that,
by the candidate's own statement, wrote most of the code and drafted the documents. That is disclosed and not
penalised; it does, however, make the candidate's command of the material the main risk, and that can only be checked
at the defence. **Score 7.0.**

**Reproducibility (10 %).**
5: code available, rerunning needs substantial undocumented effort. 7: code, data provenance and a pinned
environment; the main results rerun with path edits. 9: one-command reproduction, all dependencies pinned, verified by
a third party.
*Assessment.* Table 1, the R1 analysis, the theory check and the three figures were rerun and matched byte for byte
(Section 2.4); freeze lists match the repository; data sources carry commits and checksums. Gaps: hard-coded absolute
paths; `pyproj` and QuakeMigrate (1.2.2 per the docstring of `qm_wrapper.py`) missing from `requirements-frozen.txt`;
the PyOcto model file of R1 not in the repository; GaMMA and PhaseNet to be installed by hand; no README in
`results/taskR1/`; no license. **Score 8.0.**

### 6.2 Weighted total

| Criterion | Weight | Score (undergraduate) | Weighted |
|:--|--:|--:|--:|
| Novelty and positioning | 0.20 | 6.0 | 1.200 |
| Theoretical rigor | 0.15 | 7.0 | 1.050 |
| Experimental design | 0.20 | 8.0 | 1.600 |
| Empirical evidence | 0.15 | 6.5 | 0.975 |
| Internal coherence and writing | 0.10 | 7.0 | 0.700 |
| Feasibility | 0.10 | 7.0 | 0.700 |
| Reproducibility | 0.10 | 8.0 | 0.800 |
| **Total** | 1.00 | | **7.025, that is 70.3 / 100** |

### 6.3 The same work judged as a doctoral proposal (for transparency)

- Novelty and positioning: 4.0 (the core method is refuted on real data; the surviving element is a known robustification; key literature missing).
- Theoretical rigor: 5.0 (correct but elementary; no theory of the deflation, the grid or real noise).
- Experimental design: 6.5 (exemplary protocol, but tiny real tests, a matched synthetic generator and an underpowered final test).
- Empirical evidence: 4.5 (no positive result on real seismicity; G1-R on one day of noise).
- Internal coherence and writing: 6.0 (the same slips weigh more at that level).
- Feasibility: 5.0 (the real-data programme that a doctorate would need is not yet shown to be viable).
- Reproducibility: 7.0 (the same gaps).
- Weighted total as a doctoral proposal: 5.325, that is 53.3 / 100.

## 7. Defects, ordered by severity, with the smallest fix

Severity: **major** = changes how a reader should weigh a main claim; **moderate** = an overstatement or inconsistency
a committee would require to be corrected; **minor** = precision, presentation or reproducibility convenience. No
defect found changes a registered number or a registered decision.

| ID | Severity | Defect | Where | Smallest fix |
|:--|:--|:--|:--|:--|
| D1 | major | The bounded-influence results (E4 to E6, E8 to E10, abstract "+0.228") are under the epicentre-only rule LOC, unlabelled, while the real-data results are under LOC+T; the origin-time values for 2I (+0.031, +0.019) and 2M (+0.085 [+0.053, +0.120], +0.104) are not reported | Abstract, 1.3, Table 1, 6.2, 6.3 | Add the hit rule to the Claim column of every Table 1 row and one sentence in 6.2 and 6.3 giving the LOC+T values of 2I and 2M; say "(epicentre rule)" after +0.228 in the abstract |
| D2 | major | G1-R "confirmed on real noise ... and explained by theory" overstates: the primary is one day of noise with injected synthetic events; the comparator's recall falls to 0.003 in transient hours; dependence-aware intervals ([+0.088, +0.314]) are not reported; no comparison with standard per-station clipping or despiking, which the thesis calls common practice without a citation | Abstract, 1.3 item 2, 6.3, 7.3, 8 | In 6.3 state "one day", give the moving-block interval, and quote the verifier's "not a general win over BOMP-type detectors"; in 8 replace "confirmed on real noise and explained by theory" by "measured on one day of real noise and consistent with Propositions 1 and 4"; cite a back-projection clipping reference |
| D3 | major | Restated G1 (90 % of readings within 0.25 s) is post hoc, its threshold is arbitrary, it is not measured on the synthetic benchmark with the same procedure, and it ignores the change of comparator (STA/LTA-fed on synthetic, PhaseNet-fed on real data) | Table 2, log entry of 18 Oct 02:36 | Either run the same held-out timing measurement on the calibration seed and add the comparator condition, or relabel the row "Refuted on real data; a conjecture for future work: ..." and drop "Retained" |
| D4 | moderate | Pre-run registration is not verifiable by a third party for any task: the first commit postdates 2A to R0, the 2M registration was first committed after its result, R1's 72 s after its run started; the chat record with the hashes is not in the repository; Section 5.1 says the chain "proves that entries were not edited after they were chained" | 5.1, 7.4, Table 3 | Add an export of the chat messages that carry the hashes (with their timestamps) to `preregistration/`, or rewrite the sentence as "the chain shows that the stored log and hashes agree; their timing rests on the chat record, which is not public"; add a column "earliest external anchor" to Table 3 |
| D5 | moderate | README status paragraph says all eleven tests were run with frozen code and independently verified, contradicting Table 3 | README.md, "Status in one paragraph" | Replace the first sentence by "Eleven tests were registered before their run; which of them also had frozen code and an independent verification is listed in Table 3 of the thesis" |
| D6 | moderate | "Every confirmatory test followed the same steps" (5.1) and "Every confirmatory claim was registered before its data were analysed" (abstract) are contradicted by Table 3 and by the reuse of the 2A windows in 2B to 2E | Abstract, 5.1 | "The protocol below was followed in full from Task 2E on; Table 3 states which steps each earlier task had" and, in the abstract, "registered before its run" |
| D7 | moderate | Section 8 says the detector "beats pick-based association" on synthetic data without the qualifiers of 6.2 (STA/LTA-fed associators, QuakeMigrate the best comparator, operating points on evaluation, BOMP-PK equal) | 8 | Add one clause: "against associators fed by an STA/LTA picker and QuakeMigrate, with operating points chosen on the evaluation data, and no better than block pursuit with the same tolerance" |
| D8 | moderate | The C18i medium comparison reports the reading without the event term (0.45 vs 0.37 to 0.39); with the event term, which Y.32 itself calls the relevant one, the media are equal on N8 (0.427 vs 0.427) and differ by 0.03 on U8; the "reversal" lesson of 7.5 is overstated | 6.5, 7.5 | Report both readings in 6.5 and soften 7.5 to "the search window hid part of the effect of the medium; with origin time absorbed, the medium changes the fraction by at most 0.03" |
| D9 | moderate | The R1 velocity model is called a "28-point approximation" but 13 points are listed and used; not checked against the source | 6.6, Y.33, registration | State "the first 13 points (to 65 km) of the 28-point approximation", or correct the count; add an erratum line to the log if the registration text was wrong |
| D10 | minor | Theory wording: "leaves the Gaussian false-alarm tail unchanged" should be "changes any tail probability by at most 9.4e-6 per hypothesis (3.4e-4 per window)"; independence of the three components and exact normalisation are used but not stated; Proposition 3's title overstates (equal block norms) | 4, AA.9.4, AA.9.5 | Reword three sentences and add the two assumptions to the statement of Proposition 4 |
| D11 | minor | Proposition 1 numbers (12, 47) are the two-large-block bound; for one transient replacing a clean block the bound is about 8 (cap) and 25 (log), so "about the threshold itself" is pessimistic | 4, AA.9.2 numbers | Give r(M) + r(|v_k|) for a typical clean block next to the current numbers |
| D12 | minor | One missing observation: any non-identity radial map moves a noise-free event out of the moment-tensor subspace unless its block norms are equal, which is the mechanism of the cost in isolated windows | 4 | One sentence after Proposition 2 |
| D13 | minor | Table 3 understates 2M resumptions (GEN arms also stopped and resumed after a container restart) and omits the post-freeze change of `analyze_2d.py` (declared elsewhere) | Table 3 | Edit the two cells |
| D14 | minor | Reproducibility gaps: absolute paths; `pyproj` and QuakeMigrate 1.2.2 not in `requirements-frozen.txt` (Figure 3 fails without pyproj); PyOcto model of R1 not in the repository; no README in `results/taskR1/`; no license; the Hi-net section of `data/README.md` gives neither the target path nor the checksum command | requirements-frozen.txt, data/README.md, code/r1 | Add the two packages, add `code/r1/pyocto_jma2001.model` or its one-line build command, add a short `results/taskR1/README.txt`, choose a license, complete the Hi-net block |
| D15 | minor | Literature: no source-scanning (Kao and Shan 2004), LOKI (Grigoli et al. 2013), BackTrackBB (Poiata et al. 2016), MALMI (Shi et al. 2022) or EQTransformer; future-work item (1) is close to MALMI | 2, 8 | Add a short paragraph and four or five references; position future work (1) against MALMI |
| D16 | minor | Presentation: Table 3 precedes Table 2; mixed citation styles; bibliography numbered after the proposal; evidence-volume header stale; "within 0.007 and 0.008" order swapped; "the largest transient found" not supported by the verification; "all fed by an STA/LTA picker" for QuakeMigrate; Section 3.2 lacks the E52 qualifier | various | Copy-edit pass |

## 8. Simulated committee decision

**Composition (simulated).** Examiner A, applied mathematics and statistics; examiner B, seismology; examiner C (this
audit), methods and reproducibility.

**Deliberation.** Examiner A finds the statistics of the protocol (exact FROC, hierarchical and block bootstraps, Holm,
the simulation that replaced the non-inferiority margin) the strongest part of the work and the proofs of Section 4
correct but elementary; A would score about 72. Examiner B notes that the method the thesis set out to establish does
not work on real seismicity, that the surviving result is a known robustification measured against a comparator that
fails for a known reason, and that the literature review misses the closest migration methods; B credits the
diagnosis of the real-data failures (timing of energy peaks, S-coda ghosts, saturation) as genuinely informative and
would score about 64. Examiner C confirms that every number checked is right and that the main artefacts rerun byte for
byte, and scores 70.3 on the fixed rubric. All three agree that the negative results are reported with exemplary
candour and that the AI assistance is disclosed precisely, and that the disclosure itself is not a ground for a lower
grade. They also agree that, because the code and drafts were produced largely by an AI system, the oral defence must
establish that the candidate can (i) reproduce the proof of Proposition AA.9.5(b) at the board, (ii) explain why
BOMP-PK collapses on the 2M transient day and why that does not generalise, (iii) justify the four R1 design decisions
the log attributes to him, in particular the replacement of the non-inferiority margin, and (iv) explain the E56
deviation and why it does not change R1-a.

**Decision.** Approved (aprobado), with required corrections before final deposit: D1 to D9 (all are text changes or
small additions; none requires a new experiment, except that D3 may alternatively be resolved by one measurement on
the calibration seed). D10 to D16 are recommended. Indicative grade 70 / 100 (about 3.5 on the 0 to 5 scale), to be
confirmed by the defence; the committee would lower it if the candidate cannot answer points (i) to (iv), and would not
raise it above about 75 even with an excellent defence, because the novelty and the real-data evidence limit the work.
No distinction (meritoria or laureada) is recommended: the protocol is exceptional for this level, but the scientific
contribution that survives is modest and its strongest positive result rests on one day of noise under the lenient hit
rule.
