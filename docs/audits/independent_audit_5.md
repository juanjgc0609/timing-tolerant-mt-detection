# Independent Audit 5 of the doctoral proposal "Timing-Tolerant, Bounded-Influence Moment-Tensor Detection of Overlapping Earthquakes"

Auditor role: independent doctoral committee member, not involved in producing the work. Repository audited:
/home/claude/work/repo (clone of github.com/juanjgc0609/timing-tolerant-mt-detection, state of 16 October 2026).
Previous audits (docs/audits/, and the audit sections of docs/audit_status.md) were not read. No heavy computation
was run; only file reading, hashing and small arithmetic checks.

Rubric (fixed): Novelty and positioning 20 %; Theoretical rigor 15 %; Experimental design 20 %; Empirical evidence
15 %; Internal coherence and writing 10 %; Feasibility 10 %; Reproducibility 10 %. Each criterion 0-10; total =
weighted sum x 10.

## 1. Verification of the preregistration hash chain

Procedure: for every SHA-256 listed in preregistration/original_sealed_es/PREREGISTRO_hashes.txt, search for a byte
prefix of PREREGISTRO_1j_1k.txt that ends at an end of line and has that hash; check that the prefix lengths are
non-decreasing and that the last hash equals the hash of the whole file.

Result: PASS. All 28 hashes in the chain are SHA-256 digests of newline-terminated byte prefixes of the sealed log,
their lengths are strictly increasing (39694, 40221, 41377, 46364, 46795, 47540, 50804, 51545, 52280, 55405, 56184,
56925, 57567, 61718, 62586, 63478, 69972, 71371, 71737, 73217, 74020, 89695, 90177, 91359, 92309, 96091, 98160,
99405 bytes), and the last hash (c5cbe8e4...91b3) equals sha256 of the whole 99405-byte file. The prefix sizes
47540, 56925, 57567 and 63478 quoted in preregistration/README.md are indeed states in the chain.

Limitations found during the check (these matter for how much weight the preregistration can carry):

1. The first hash (prefix 39694 bytes) ends at line 431 of the log, i.e. inside the registration of Task 2D. The
   registrations AND results of Tasks 1j to 2C (log lines 1-430, including Task 2A, 2B and 2C, which carry the
   headline "+0.09 to +0.49" claim) were written before any hash existed. For 2A-2C the chain certifies only that
   the text has not changed since 26 Sep 2026 21:13 UTC, not that the registration preceded the run. The
   preregistration/README.md sentence "Every confirmatory task (2A to 2L, R0, 2M) was written into an append-only log
   before its run" is therefore attested by the chain only from 2D onward (and by the conversation-delivered hashes
   only from 2E onward, as the README itself says).
2. The chain is self-maintained in the same working environment: timestamps are written by the same process that
   writes the log, and there is no external anchor (public commit, OSF/AsPredicted registration, RFC 3161 or
   OpenTimestamps). The hashes delivered "to the author in the conversation" from 2E on are the only external
   witness, and they are not reproduced in the repository in a checkable form.
3. PREREGISTRO_hashes.txt contains one timestamp line without a hash ("Mon Oct 5 02:36:52 UTC 2026" immediately
   followed by "Sat Oct 10 00:00:50 UTC 2026"), so the line structure of the chain file is not perfectly regular.
   Cosmetic, but a committee reviewer checking by hand will notice it.
4. The 2M registration (prefix 89695 bytes, 16 Oct 16:32 UTC) and its result (96091 bytes, 19:59 UTC) and
   verification (98160 bytes, 20:22 UTC) all fall within roughly four hours on the same day as this audit.

## 2. What the proposal claims, as I read it

The proposal (docs/proposal/Doctoral_Proposal.md, 4455 lines, about 56,000 words) and its evidence volume
(Supplementary_Evidence_Appendices_G_to_AA.md, 7920 lines, about 109,000 words) defend one primary hypothesis, G1:
a timing-tolerant moment-tensor detector with deflation, applied to three-component waveforms before picking, detects
more events than published associators and migration locators at an equal realized false-event rate (Doctoral_Proposal.md
l. 863-872, master table l. 1135). Three methodological hypotheses are attached: G1-M (off-grid group BLASSO with
station time shifts beats block OMP with equal tolerance; refuted, l. 1136), G1-R (bounded influence per station, by a
cap or a noise-scaled Huber fit, raises recall in dense sequences; confirmed twice on synthetics and once with real
noise, l. 1137) and G1-C (pick-free station terms raise origin-time-matched recall; confirmed once on synthetics,
l. 1138). The first registered real-data test (R0, one hour of Ridgecrest 2019) was strongly negative (l. 511-521;
Y.30). The thesis statement (l. 177-181) and the closing of Section 1.6 (l. 532-535) now say that what is defended is
"the benchmark, its positive and negative results and the measured conditions under which pick-free physics helps".

The work is unusually honest: negative results are registered and reported, every positive is followed by an
adversarial verification that narrows it, and 49 implementation and design errors are logged
(docs/artifact_manifest_and_error_log.md, E1-E49). That honesty is the strongest asset of the project and I credit it
under design and reproducibility. It does not, however, substitute for a contribution that survives contact with real
data, and the documents have not yet been restructured around what the evidence now supports.

## 3. Verification of numbers stated in the text against result files

Each row: the claim, where it is stated, the result file and line, and the verdict. Line numbers refer to the files as
they are in the repository on 16 October 2026.

| # | Claim (where) | Result file : line | Verdict |
|:--|:--|:--|:--|
| 1 | R0 pooled difference -0.383 [-0.481, -0.310] (README l. 18, 36; Proposal abstract l. 79; Y.30 l. 6851) | results/taskR0/analysis_eval.txt : 13 | matches |
| 2 | Thesis detector (M2) recalled "at most 0.022" at 30 false/h, LOC+T (README l. 16; abstract l. 78) | analysis_eval.txt : 5 (M2: 0.022, 0, 0, 0) | matches (note: M1 reaches 0.044 on N8 on the same line) |
| 3 | PhaseNet + PyOcto "up to 0.667" (README l. 17; abstract l. 79) | analysis_eval.txt : 5 (PyOcto N6 0.667) | matches |
| 4 | Per-subnet N8 -0.578, N6 -0.667, N4 0.000, U8 -0.289 (results/taskR0/README.txt) | analysis_eval.txt : 13 | matches |
| 5 | M6 - M2 = +0.000 [-0.016, +0.019]; LOC rule pooled -0.450 [-0.563, -0.385] (Y.30 l. 6852-6853) | analysis_eval.txt : 43 and 16 | matches |
| 6 | Timing errors 0.35 to 0.41 s after corrections (README l. 19; Y.30 l. 6842) | results/taskR0/R0_DEV_LOG.txt : 8, 31 (0.358 and 0.412 s; "0.35-0.41 s") | matches |
| 7 | Mislocation 8 to 18 km in median, 12.8 km for M2 on N8, 1.5 km for PyOcto (abstract l. 80; Y.30 l. 6866-6867) | results/taskR0/R0_DEV_LOG_POSTFREEZE.txt : 20-21 (17.9, 12.8, 9.6, 18.2, 8.2; PyOcto 1.5) | matches |
| 8 | Semivariogram 0.25 to 0.42 s^2 with source separation (Y.30 l. 6871) | R0_DEV_LOG_POSTFREEZE.txt : 4 (0.254 ... 0.421) | matches |
| 9 | 2M-a cap - BOMP-PK +0.228 [+0.148, +0.298]; 2M-b HUBER-N - BOMP-PK +0.241 [+0.179, +0.293] (README l. 20, 37; exec. summary l. 15-16; Y.31 l. 7062) | results/task2M/analyze_2m_out.txt : 12-13, 343, 345 | matches |
| 10 | Recalls 0.325 / 0.339 / 0.097 for cap / HUBER-N / BOMP-PK (Y.31 l. 7063-7064, Table Y.91) | analyze_2m_out.txt : 7, 11, 10 | matches |
| 11 | HUBER-N - cap +0.014 [-0.038, +0.056]; origin-time gains +0.085, +0.104; standard windows +0.423, +0.346 (Y.31 l. 7064-7066) | analyze_2m_out.txt : 14, 49-50, 86-87 | matches |
| 12 | Secondary days +0.057 [+0.004, +0.121] and +0.120 [+0.070, +0.168]; HUBER-N - cap +0.062 [+0.012, +0.107] (exec. summary l. 17; Y.31 l. 7067-7068) | analyze_2m_out.txt : 181-183, 344, 346, 357 | matches |
| 13 | Post hoc: cap gain +0.033 at night, +0.300 from 7.5 to 15 h; HUBER-N +0.123 at night (exec. summary l. 19; Y.31 l. 7082-7084) | results/task2M/posthoc/blocks_out.txt : 1-2 | matches |
| 14 | Event-free windows: 0.000 / 0.028 / 0.042 false events per window for cap / HUBER-N / BOMP-PK (Y.31 l. 7089-7090) | analyze_2m_out.txt : 162, 166, 165 | matches |
| 15 | C18a fraction within 0.25 s: N8 0.59, 0.60, 0.58, 0.59, 0.59; U8 0.61-0.63; N6 0.53-0.59; N4 0.53-0.55; with event term 0.62 / 0.73 / 0.59 / 0.54 (Y.32 l. 7133-7137; exec. summary l. 22) | results/explorations/c18/c18a_summary.txt : 1-4 | matches (Y.32 says "130 readings" on N8; the file reports n = 128 for the event-term column; immaterial) |
| 16 | C18b/c AUC with cap 0.49, 0.39, 0.43, 0.43; no cap 0.71, 0.69, 0.72, 0.74; cap 12 0.61-0.64; cap 24 0.69-0.74; LOG6 0.66-0.72; argmax within 5 km 6-31 % (Y.32 l. 7150-7158; exec. summary l. 24-25) | results/explorations/c18/c18c.log (last 24 lines) and c18b.log (last 16 lines) | matches |
| 17 | C18d cap 24 minus cap 6 on N8, LOC rule, +0.125 [+0.043, +0.300]; on U8 -0.062 (Y.32 l. 7163-7164) | c18/c18d_analysis.txt : 6, 10 | matches |
| 18 | C18e cap +0.148 [+0.056, +0.261], LOG6 +0.133 [+0.052, +0.224], HUBER-N +0.175 [+0.095, +0.265], cap 24 +0.066 [+0.010, +0.141]; LOG6 - cap -0.015 [-0.099, +0.043]; cap 24 - cap -0.082 [-0.183, -0.014]; LOC+T LOG6 - cap +0.015 [-0.022, +0.060] (Y.32 l. 7170-7173) | c18/c18e_analysis.txt : 9-14, 28 | matches |
| 19 | 2A: "+0.09 to +0.49 in all 12 cells" over the best of GaMMA, PyOcto, QuakeMigrate (README l. 31) | results/task2a/analyze_2a_out.txt : 61, 137 (cap 0.1: 0.09 ... 0.49) | matches as a number |
| 20 | Same range attributed to "eight preregistered and independently verified tasks" (README l. 11-13) and "in eight preregistered tasks, such detectors exceeded GaMMA, PyOcto and QuakeMigrate by 0.09 to 0.49 in every cell" (abstract l. 65-66) | analyze_2a_out.txt is the only output that compares with GaMMA and PyOcto; 2C/2D compare with QuakeMigrate only; 2E-2L compare with BOMP-PK (the proposal itself says so at l. 1135); verification folders exist only for 2A, 2E, 2F, 2I, 2L, 2M (results/*/verif_*), none for 2B, 2C, 2D | does NOT match: the range comes from one task, and three of the eight tasks have no independent verification folder |
| 21 | 2F: +0.072 and +0.074 (thin), +0.051 and +0.043 (strong), 96 networks per cell (Proposal l. 345-347) | results/task2f/analyze_2f_out.txt : 6, 26-30 | matches |
| 22 | 2I: +0.058 [0.039, 0.076], +0.053 [0.037, 0.070]; non-inferiority -0.008, -0.002 with bounds -0.023, -0.018; absolute dense recall 0.235 and 0.251 (Proposal l. 371-386) | results/task2i/analyze_2i_primary_out.txt : 82-88, 24, 29 | matches |
| 23 | 2L: +0.057 [+0.034, +0.081], +0.062 [+0.036, +0.083] standard; +0.038, +0.030 dense; BOMP +0.068, +0.071 (Proposal l. 481-484; exec. summary) | results/task2l/analyze_2l_out.txt : 113-119, 146-148 | matches |
| 24 | 2E: +0.000 to +0.039, all intervals include 0 (README l. 32; Proposal l. 329-331) | results/task2e/analyze_2e_out.txt : 67, 160 | matches |
| 25 | 2D-a: -0.155 to +0.091 (Proposal l. 313) and 2C-a +0.08 to +0.52 (l. 297-299) | results/task2d/analyze_2d_out.txt : 83, 146; results/task2c/analyze_2c_out.txt : 44, 83 (0.077 to 0.519) | matches |
| 26 | "forty-five errors were documented (E1 to E45)" (executive_summary.md l. 827) | docs/artifact_manifest_and_error_log.md : 973-1023 (E46-E49 exist) | does NOT match (stale) |
| 27 | Y.32 code list includes `c18e_analyze.py` (Supplement l. 7117) | not present in code/explore_contrib/ (c18e_run.py and c18d_analyze.py are) | does NOT match (missing file) |
| 28 | Y.30 cites output `r0/run/analysis_eval.txt` (Supplement l. 6822) | the file is at results/taskR0/analysis_eval.txt | path mismatch (minor) |

Summary: 24 of 28 checked items match their result files exactly; the numerical pipeline from outputs to text is
reliable. The four mismatches are attribution, staleness and file-inventory errors rather than wrong numbers, but item
20 sits in the first paragraph of the README and of the abstract, which is where a committee reads first.

## 4. Scores on the fixed rubric

### 4.1 Novelty and positioning (weight 20 %): 4.0 / 10

What is new is thin once the project's own ablations are applied, and the proposal says so itself.

- The deployed detector is block orthogonal matching pursuit over moment-tensor blocks (Vera Rodriguez, Sacchi and Gu
  2012, ref. [129]) with a +-0.25 s peak search. Task 2D-a failed in both regimes: with equal tolerance BOMP-PK matches
  MT-FW (-0.155 to +0.091, results/task2d/analyze_2d_out.txt l. 83, 146), and the proposal concedes that "the advantage
  of the detector over [129] and [122] is attributed to the peak search alone" (Doctoral_Proposal.md l. 311-315).
- The one genuinely new mathematical object, the off-grid group BLASSO with station time shifts (G1-M), was refuted
  (Task 2E, +0.000 to +0.039, analyze_2e_out.txt l. 67, 160). The theory of Appendix AA is the theory of that refuted
  object (Doctoral_Proposal.md l. 416-421, 1001-1002).
- The confirmed element, G1-R, is per-station clipping or a Huber fit, which the proposal itself calls "known practice in
  robust stacking" and "standard robust statistics" (l. 445-448). Task 2F shows that adding the same cap to BOMP-PK
  (BOMP-PK-CAP+SC 0.174/0.274 against MT-FW+SC 0.187/0.289 in thin 4242, results/task2f/analyze_2f_out.txt l. 12)
  recovers most of the gain, so the effect belongs to clipping, not to MT-FW.
- G1-C (pick-free station terms) is, by the proposal's own account, the pick-free version of classical joint
  hypocentre and station-term determination (l. 491-493); its gain is confined to one network geometry.
- The positioning against published associators rests on a synthetic comparison in which a calibrated STA/LTA picker
  "starves" GaMMA and PyOcto (l. 282-283; risk table l. 3391); on real data with PhaseNet the order reverses by 0.38
  (results/taskR0/analysis_eval.txt l. 13).

Positive: the literature review is broad (Section 3, 140+ references), the gap is articulated, and the reframing of the
contribution as a preregistered benchmark with its negatives is honest and defensible as a methods contribution to
seismology. For a doctorate in applied mathematics, however, the novel content that survives is a measured design
result about clipping inside a known pursuit, plus a negative result. That is a 4.

### 4.2 Theoretical rigor (weight 15 %): 5.0 / 10

- Appendix AA (Supplement l. 7191-7920) is careful: every statement lists its hypotheses, proofs are given, conjectures
  are labelled (Conjecture AA.1, l. 7861), and AA.8 (l. 7892-7919) states plainly what is not covered. Proposition
  AA.2 (dual certificate at zero equals the GLRT), AA.4 (Frank-Wolfe on the grid, with the implemented deflation
  covered only by parts e-h) and AA.5 (false-alarm control, with only the conformal part c surviving real noise) are
  correct in form and appropriately scoped.
- The results are largely adaptations of known off-the-grid sparse recovery theory (BLASSO, sliding Frank-Wolfe, group
  non-degeneracy) to 6-dimensional blocks; I did not find a result that a referee in sparse recovery would regard as new.
- The theory does not cover the statistic that is actually deployed and confirmed: neither the peak search, nor the
  cap, nor the Huber fit, nor LOG6. Lemma AA.3 (bounded influence of one station, l. 7447) is elementary. AA.8 admits
  that "nothing here bounds false events in windows that contain events, which is what the FROC endpoint counts"
  (l. 7918-7919). The theory track for the deployed statistic is promised (Section 19, T2, Doctoral_Proposal.md l. 3447)
  but not delivered.
- Sections 5, 6, 9 and 10 (elastodynamics, five constraint families, identifiability, differentiable wavefield) are a
  programme of mathematics loosely connected to the evidence; the proposal acknowledges that the prototype tested is not
  the model (l. 2976-2982).
- No theory explains the most important empirical finding of the project, the saturation of the capped score in
  continuous seismicity (C18b/c AUC 0.39-0.49) versus its robustness to transients (2M). That trade-off is exactly the
  kind of question a robust-statistics analysis (influence function, breakdown, power under a contaminated alternative)
  should answer.

### 4.3 Experimental design (weight 20 %): 6.5 / 10

Strengths, which are real and above the norm for a proposal:
- Append-only registration before runs, code frozen by SHA-256, fresh seeds per task, hierarchical bootstrap over
  networks and windows, Holm correction, explicit non-inferiority and equivalence margins, and decision rules written
  in advance (e.g. PREREGISTRATION_LOG_EN.txt l. 251-311 for 2A, l. 776-923 for 2M). The hash chain verifies (Section 1).
- Registered ablations that could and did refute the author's preferred story (2B, 2D-a, 2E), and adversarial
  verifications that narrowed positives (origin-time rule, scope to one cell in 2L).
- R0 was a genuine out-of-sample real-data test with a development/evaluation split fixed before the waveforms were
  read (log l. 676-700), and the C18 gate for the next real test is fixed in advance (Supplement l. 7183-7187).

Weaknesses that a committee would press:
- Operating points are chosen on the evaluation windows in every synthetic task (Doctoral_Proposal.md l. 1091-1097);
  the deployable rule 2A-e fails in its strict reading and the conformal thresholds of 2D control false alarms only in
  empty windows (1.2 to 2.2 false events per window with events). No result at a deployable operating point exists.
- Comparator fairness. In 2A the pick-based associators are fed an STA/LTA picker that discards most arrivals; the best
  comparator was QuakeMigrate in all 24 cases (log l. 316-317). In 2M the comparator is BOMP-PK with unbounded energies,
  which collapses to a recall of 0.097 overall and 0.003 in the daytime block (results/task2M/posthoc/blocks_out.txt
  l. 2); no robustified baseline (per-station noise normalization, AGC, clipped BOMP, which existed as BOMP-PK-CAP in
  2F) was registered. The verification itself notes that the comparison "is a capped statistic against unbounded
  energies under one common threshold ... and not a general win over every detector of the BOMP type" (Y.31
  l. 7101-7103).
- 2M primary day selection. The December day was made the sole primary replicate after its noise was characterized as
  transient-rich (Y.31 l. 7031-7039), with the mechanism (transients sink BOMP-PK) already known from C16. That is a
  legitimate test of a mechanism, but it is designed to confirm; the cells share hours, so the intervals are too narrow
  (Y.31 l. 7111-7112).
- One synthetic generator family for all synthetic tasks (l. 442-443); the promised second generator (Section 19) has
  not been run.
- R0 is one hour, 45 target events, one network, a homogeneous medium and a scalar instrument response; the N4
  subnetwork has zero recall for every method at the primary operating point (analysis_eval.txt l. 5), so it contributes
  zero to the pooled contrast and dilutes it.
- The C18 gate (5 of 16 CAL events with origin time at 30 false/h on N8 or U8) is defined on the same 15 minutes that
  were used for R0 development, C17 and C18, which invites overfitting to 16 events, and the threshold (0.31) is about
  half of what PyOcto achieved on N8 in EVAL (0.600). Passing it would not make G1 likely to hold.
- 2A-2C have no hash attestation of registration before the run (Section 1 above), and 2B, 2C, 2D have no independent
  verification folder.

### 4.4 Empirical evidence (weight 15 %): 4.0 / 10

- Primary hypothesis G1: supported only on synthetics, at operating points chosen on evaluation, against starved
  pick-based baselines; its first real-data test was negative by a large and precise margin, -0.383 [-0.481, -0.310],
  with the thesis detector at 0.022 recall against 0.600-0.667 for PhaseNet + PyOcto on N8 and N6
  (analysis_eval.txt l. 5, 13). Among all pick-free arms recall lay between 0.000 and 0.044 at the primary rule.
- G1-R: two synthetic confirmations of modest size (+0.04 to +0.07; +0.02 to +0.045 with origin times matched; absolute
  dense recall 0.22 to 0.29) and one confirmation with real noise and synthetic signals against an unbounded comparator.
  Its real-noise gain is a property of transient-rich hours (+0.033 [-0.044, +0.128] at night, blocks_out.txt l. 1).
- On real seismicity the confirmed element is harmful: with the cap, the score at the true hypocentre of a real M >= 2
  event is not larger than the maximum of a null window (AUC 0.39-0.49, c18c.log), while without the cap it is
  0.69-0.74; and BOMP-PK (no cap) placed atoms closer to the events than the thesis detector on N8 (median 8.2 km
  against 12.8 km, R0_DEV_LOG_POSTFREEZE.txt l. 21; LOC recall at 30/h 0.089 against 0.022, analysis_eval.txt l. 8).
  The project's positive (robustness to transients) and its real-data failure (saturation in seismicity) are two faces
  of the same clip, and at present only the failure has been observed on real earthquakes.
- G1-C: one synthetic confirmation, concentrated in one geometry, conditional on a catalogue of known error.
- C18 shows that the R0 failure is not fixed by a better medium or a better score; there is no candidate yet that
  detects real events.

The evidence is high quality in the sense of being honest and reproducible; as support for the thesis hypothesis it is
weak, and the strongest single datum (R0) is against it.

### 4.5 Internal coherence and writing (weight 10 %): 4.0 / 10

- Length: about 165,000 words across the proposal and its supplement, plus 63,000 more in the executive summary, method
  design and task plan (wc -w). A committee cannot read this. Section 1.6 alone (l. 257-535) is a 280-line narrative
  of task history with paragraphs of 30-40 lines.
- Stale or contradictory parts after R0:
  - The master table status of G1 (l. 1135) still reads "Numerical evidence, synthetic only" and does not mention R0.
  - Section 13.0c "Real-data arm" (l. 3030-3034) and Section 16 (l. 3312-3362) still describe the RSNC decimation test
    G1' as the decisive next step with catalogue station corrections and the full operator, without the C18 gate.
  - Section 19, months 1-6, experimental track E1 (l. 3446) still schedules "the first real-data decimation test of G1
    on the RSNC with MT-FW and BOMP-PK" and the gate "first real-data decimation result for G1 reported by month 6,
    whatever its sign" (l. 3455), whereas Section 1.6 (l. 534-535) and Y.32 (l. 7181-7187) say the next real test
    waits for a candidate to pass the gate.
  - executive_summary.md §8 "Framing recommendation" (l. 805-823) still speaks of "its pending test on real data" and
    "Tasks 2A to 2F"; §9 (l. 827) counts 45 errors while the register has 49; §7 (l. 767-805) labels 2B, 2C, 2D
    "EXECUTED AND VERIFIED" although no independent verification output exists for them.
  - README l. 11-13 and abstract l. 65-66 attribute the 2A range to eight tasks (item 20 of Section 3).
- Structure: in the supplement, Y.31 and Y.32 (l. 7015, 7114) are placed after Appendix Z (l. 6889), breaking the
  lettering; the title still advertises "a Real-Data Test Plan" on a "decimated national network" (l. 10-15) while the
  first real test was on a different network and was negative.
- Positive: terminology is consistent, every number I checked is traceable, limitations are written next to claims,
  and the evidence taxonomy (l. 1026-1037) is a good device.

### 4.6 Feasibility (weight 10 %): 5.0 / 10

- The student has shown very high throughput (eight synthetic preregistered tasks, two real-data tasks, 18
  explorations) and the experimental core runs on two CPU cores; the immediate tasks are cheap.
- The 48-month plan (Section 19, Table 8, l. 3444-3453) still carries an elastodynamic solver (SPECFEM3D or SeisSol),
  Gaussian beams, a strain Green's tensor database, an optimal-transport misfit, a symmetry-exact learned
  representation and uncertainty theory, none of which the evidence motivates any longer; Section 1.5 admits that no
  evidence lies in the declared operating regime (l. 248-255). The risk row "Scope too large for the schedule"
  (l. 3401) is answered only by "two-track plan".
- The blocking scientific problem is identified (reading of arrivals in a dense sequence, Y.32 l. 7176-7187) but there
  is no candidate that works, and the bar to beat on real data (PhaseNet + PyOcto at about 0.6 recall) is high. The
  realistic path to a positive G1 on real data is uncertain within the plan.
- Dependence on a supervisor whose profile alignment is argued in Section 21, and on RSNC data access and a validated
  Colombian velocity model (Section 16 gate), adds risk that the plan does not cost.

### 4.7 Reproducibility (weight 10 %): 7.0 / 10

- Strong: code for every task and exploration, frozen code hashes (code/task2d..2l_code_hashes.txt,
  code/r0/R0_FREEZE.txt, results/task2M/T2M_FREEZE.sha256), analysis outputs for every number checked, data provenance
  with commits and checksums (data/README.md), a verifiable hash chain, pinned library versions and the BLAS kernel
  (README l. 66-68), an error register with 49 entries, and independent verification scripts for six tasks.
- Weak: paths are hard-coded to /home/claude/work (README l. 70-72); no lockfile, container or environment file;
  185 scripts in a flat layout with Spanish comments; no license (README l. 83); raw outputs of R0 and C16 were lost
  (E49) and the R0 contrast M1 - M3 must be re-run; `c18e_analyze.py` cited in Y.32 is absent; bit-for-bit
  reproduction depends on the OpenBLAS kernel (E31); the registry has no third-party seal; the git history starts on
  11 October 2026 (after every task except 2M), and the 2M registration was committed only after its result (commit
  ca4185e, 2026-10-16 19:02 -0500), so the public repository could have sealed 2M and did not.

### 4.8 Total

| Criterion | Weight | Score | Weighted |
|:--|--:|--:|--:|
| Novelty and positioning | 0.20 | 4.0 | 0.80 |
| Theoretical rigor | 0.15 | 5.0 | 0.75 |
| Experimental design | 0.20 | 6.5 | 1.30 |
| Empirical evidence | 0.15 | 4.0 | 0.60 |
| Internal coherence and writing | 0.10 | 4.0 | 0.40 |
| Feasibility | 0.10 | 5.0 | 0.50 |
| Reproducibility | 0.10 | 7.0 | 0.70 |
| **Total** | | | **5.05 x 10 = 50.5 / 100** |

## 5. Critical defects and the smallest action that fixes each

Ordered by how strongly each would drive a committee to require major revision.

**D1. The primary hypothesis was tested on real data, failed by a wide margin, and has not been reformulated.**
G1 remains "the only hypothesis whose outcome can make the thesis fail" (Doctoral_Proposal.md l. 1099-1101) and "If G1
fails on real data, the thesis fails as a claim about detection" (l. 997-998). R0 is that failure on its first real
contact (-0.383 [-0.481, -0.310]); the proposal reports it in Section 1.6 and the abstract but keeps G1 as the primary
hypothesis with status "Numerical evidence, synthetic only" (l. 1135) and keeps the RSNC decimation as the month-6 test
(l. 3446, 3455). *Smallest fix:* rewrite G1 and its master-table row so that its status includes R0, and choose
explicitly between (a) demoting G1 to a conditional hypothesis that is tested only after the C18 gate is passed, with
the benchmark and G1-R/G1-C as the defended core, or (b) restating G1 for a regime where pick-based association is
expected to fail (sparse networks, low SNR, isolated events) and registering a real test in that regime. Update
Sections 13.0c, 16 and 19 accordingly. This is a text change of a few pages.

**D2. The surviving methodological novelty is small and is not covered by the theory.**
The detector equals known block OMP with a peak search (2D-a), the novel convex method is refuted (2E), and the
confirmed element is clipping/Huber weighting, transferable to BOMP (2F, BOMP-PK-CAP). Appendix AA is the theory of
the refuted object. *Smallest fix:* deliver the first piece of track T2 now: a short appendix that states the capped,
Huber and LOG6 statistics as M-estimators of the block amplitude, derives their influence functions and their
null-distribution tail (for the false-alarm threshold), and proves or at least quantifies the trade-off observed
between robustness to single-station transients (2M) and saturation under many simultaneous arrivals (C18b/c). That
would turn the empirical G1-R into a mathematical contribution and connect AA to what is deployed.

**D3. The real-noise confirmation of G1-R (Task 2M) uses a comparator that no practitioner would deploy.**
BOMP-PK with unbounded energies collapses to 0.097 recall (0.003 in the daytime block) on the primary day; the gain
of +0.228 is largely the gain of not being destroyed by transients. *Smallest fix:* on the frozen 2M parts and noise,
run, declared as post hoc and descriptive, BOMP-PK with per-station noise normalization (rolling MAD or STA/LTA
normalization) and BOMP-PK with the same cap; report the cap and HUBER-N against those. Register that comparison as
the primary contrast of the next G1-R test.

**D4. The synthetic evidence for G1 against published associators is not a fair or deployable comparison.**
GaMMA and PyOcto are starved by an STA/LTA picker in 2A; QuakeMigrate was the best comparator in all cells; operating
points are chosen on evaluation; with optimistic picks GaMMA leads at 0.1 in four of six weak-clutter cells; and on
real data PhaseNet-fed associators win. Yet "+0.09 to +0.49 over GaMMA, PyOcto and QuakeMigrate" is the first number
in the README and the abstract, attributed to eight tasks. *Smallest fix:* remove that sentence from the README and
abstract or restate it as "in one preregistered synthetic task (2A), with an STA/LTA picker that starves the
pick-based associators and operating points chosen on evaluation", and correct the "independently verified" wording
for 2B-2D.

**D5. The documents are too long and partly stale.**
About 165,000 words in proposal and supplement; stale statements in the master table, Sections 13.0c, 16 and 19, and
the executive summary §7-§9; Y.31 and Y.32 placed after Appendix Z. *Smallest fix:* cut the main proposal to at most
60 pages (problem, the hypotheses as they now stand, a one-page evidence table generated from the result files,
methods of what is deployed, the next tests, a 36-48 month plan), move the task-by-task history into the supplement,
and fix the specific stale passages listed in Section 4.5.

**D6. The preregistration is not externally sealed, and the opportunity to seal Task 2M was missed.**
The hash chain is internally consistent but self-timestamped; 2A-2C precede the first hash; the public repository
existed from 11 October but received the 2M registration only after its result. *Smallest fix:* before the next
registered task, push the registration as a signed commit or tag to the public repository (or an OSF / AsPredicted
registration) before any run, and stamp the current chain head (c5cbe8e4...91b3) with OpenTimestamps or an RFC 3161
service today; correct preregistration/README.md to say that hash evidence of prior registration exists from 2D on.

**D7. The plan is not feasible as written and is not aligned with the blocking problem.**
The 48-month plan keeps the PDE solver, beams, learned representation and transport branches although Section 1.5
admits no evidence lies in their regime, and does not budget the one problem that blocks real-data transfer (reading
of arrivals in dense sequences). *Smallest fix:* move families F1-F5 and the learned representation to "optional
extensions", and make the first 12 months: gate work on the reading of arrivals; a second real domain; the T2 theory
of the bounded-influence statistic; the second synthetic generator; and the registered real test conditional on the
gate.

Non-critical but expected to be raised: the missing `c18e_analyze.py`; hard-coded paths and no lockfile; no license;
the lost raw outputs of R0 and C16 (E49); the N4 subnetwork of R0 with zero recall for every method; AUCs in C18
reported without intervals on 16 events; null windows in C18b/c that certainly contain uncatalogued events below M 2.

## 6. Simulated decision

**Decision: major revision.**

A reasonable committee would not reject, because the candidate has demonstrated the core competences of independent
research to an unusual degree: formulating falsifiable hypotheses, registering them, running fair ablations against
their own method, reporting negatives, and testing on real data early. It would not accept, because the primary
hypothesis has failed its first real-data test and has not been reformulated, the methodological novelty that
survives is small and outside the theory, and the documents are far too long and partly inconsistent. Minor revision
is not plausible because D1, D2 and D7 require rethinking what the thesis claims, not editing.

My estimate of what a reasonable committee would give on this rubric is 50 to 58 out of 100, centred near 53; members
from seismology would weigh the honest real-data test and the benchmark more (upper end), members from applied
mathematics would weigh the absence of theory on the deployed statistic more (lower end).

## 7. Actions that would raise the score most toward 75/100, in order

Gains are in rubric points (out of 100) and assume the action is done well; they are estimates, not promises, and
several depend on outcomes.

| # | Action | Criteria affected | Estimated gain |
|:--|:--|:--|:--|
| 1 | Restructure the proposal around what the evidence supports: reformulate G1 after R0 (D1), correct the README/abstract attribution (D4), cut the main text to at most 60 pages with a single evidence table generated from the result files, and fix every stale passage listed in Section 4.5 and D5. | Coherence 4 -> 7 (+3.0); Novelty/positioning 4 -> 4.5 (+1.0); Feasibility 5 -> 5.5 (+0.5) | about +4.5 |
| 2 | Re-plan Section 19 (D7): drop or defer the PDE/learned/transport branches, centre months 1-12 on the arrival-reading gate, a second real domain, the T2 theory and the second generator, and make the real test conditional on the gate with explicit fall-backs. | Feasibility 5 -> 7 (+2.0); Coherence +0.5 (+0.5) | about +2.5 |
| 3 | Deliver the theory of the deployed statistic (D2): M-estimator form, influence functions and null tail of cap, Huber and LOG6 with peak search; a quantitative account of robustness to single-station transients versus saturation under many arrivals; a false-alarm statement usable for the threshold. | Theory 5 -> 7 (+3.0); Novelty 4 -> 5 (+2.0) | about +5.0 |
| 4 | Fair comparator for G1-R (D3): on frozen 2M data, BOMP-PK with per-station noise normalization and with the cap, reported as post hoc; then a registered G1-R test with that comparator as primary, on the second generator family and on a second real-noise network. | Design 6.5 -> 7.5 (+2.0); Evidence 4 -> 4.5 (+0.75) | about +2.75 |
| 5 | Pass the C18 gate with a new reading of arrivals and run the registered real test on the unseen Hi-net swarm against PhaseNet + PyOcto, with a deployable operating point frozen on development data and false events counted over all windows. A positive or a narrowed non-inferiority would change the evidence score most; an honest negative still helps design and coherence. | Evidence 4 -> 6 (+3.0) if positive, +1.0 if negative; Novelty +1 to +2 (+2 to +4) if positive; Design +0.5 (+1.0) | +2 to +8 (outcome-dependent) |
| 6 | External sealing and packaging (D6 and Section 4.7): OpenTimestamps the current chain, register future tasks publicly before running, add a lockfile or container, relative paths, a license, a one-command script per table, and restore `c18e_analyze.py`; archive raw outputs outside the container. | Reproducibility 7 -> 8.5 (+1.5); Design +0.5 (+1.0) | about +2.5 |
| 7 | A deployable operating point on synthetics: thresholds frozen on held-out networks for a false-event budget over all windows (events included), with the decoy calibration of C13 or a conformal rule on event windows, registered. | Design +0.5 (+1.0); Evidence +0.5 (+0.75) | about +1.75 |

Taking actions 1-4 and 6-7 (all under the student's control and none dependent on a lucky outcome) gives roughly
50.5 + 19 = 69 to 70. Reaching 75 requires, in addition, a real-data result that is at least non-inferior to PhaseNet +
PyOcto in some declared regime, or a theoretical contribution on the bounded-influence statistic strong enough to stand
on its own (action 3 done at publication level would be worth more than estimated above).

## 8. Note on scope of this audit

R0 closure (M1 and M3 re-run, declared post hoc) is in progress and not reported; I did not penalise its absence beyond
noting that the contrast with BOMP-PK on R0 is incomplete (Y.30 l. 6884-6885) and that the raw outputs were lost (E49).
No computation was re-run; numbers were checked against stored analysis outputs only. I did not audit Appendices G to
X, Y.1 to Y.17, or method_design.md and tasks_and_real_data_plan.md beyond the passages cited.
