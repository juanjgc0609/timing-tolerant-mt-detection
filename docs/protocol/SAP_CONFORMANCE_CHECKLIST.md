# Conformance checklist: does the analysis match its registration?

Used by the review thread on every analysis that claims ADR-0008 (G4, R2, R3) before its result enters a document.
Each line is answered PASS, FAIL or N/A with the evidence (file and line, command and output). A FAIL on a line marked
[STOP] means the registered decision cannot be reported as made; other FAILs are deviations to log (protocol,
Section 8). Report with `docs/protocol/templates/verification.md`, named as a reimplementation by a separate AI
session (CR-19).

## A. Anchoring and order
1. [STOP] The registration file's SHA-256 equals the one in the tag `reg/<task>/registration-<n>`; the Release
   `created_at` (GET /repos/{owner}/{repo}/releases/tags/{tag}) precedes the first archived output of the run.
2. The `.ots` proof verifies (`ots verify`); for R2 and R3 the OSF registration exists and its frozen copy matches.
3. [STOP] The freeze tag precedes the run; every file of the freeze list hashes to its listed value.
4. Synthetic seeds equal the derivation of ADR-0004 item 4 from the tag commit and Release time.

## B. Data
5. [STOP] No evaluation waveform was read before the freeze (logs, file access times, development read ranges);
   margins included (E56).
6. Causal split: every estimate from the test network uses records ending at least 1 h before evaluation (SAP 3.2).
7. Eligibility was produced by the frozen script; excluded station-hours, events and blocks match its output (SAP 3.3).
8. Regimes and cells were computed from the catalogue or injection design, not from method outputs (SAP Section 6).
9. Pretrained-weight provenance table present; overlaps declared (SAP 3.6).
10. The use ledger shows no earlier configuration on the evaluated segments (SAP 3.5).

## C. Methods
11. [STOP] Every method's configuration in the run equals its frozen configuration.
12. Tuning ledger: same number of configurations per method, same seed list, same development data and objective;
    ledger hash in the freeze (ADR-0008 item 9).
13. Comparator set equals the registered set, including the model's-picker arm if the picker was fine-tuned.

## D. Scoring
14. [STOP] Labels recomputed independently with LOC+T-V (5 km, 0.5 s, dominance 0.5, one to one, score order) agree
    with the archived labels.
15. [STOP] The duplicator control was run and did not gain recall in T at one-sided 0.025.
16. Operating points recomputed at floor(budget x hours) false detections, ties grouped, agree.
17. Prevalence weights equal the registered table.

## E. Statistics
18. [STOP] Blocks as registered (definition, count); each sign-flip p-value reproduced (exact or with the registered
    seed and number of sign vectors), the IUT p-value is the largest pairwise one.
19. [STOP] Fixed sequence respected: H2 reported as tested only if H1 was rejected; secondary Holm family only after
    both; no hypothesis added, dropped or reordered.
20. Intervals: stratified block bootstrap with the registered replicates and seed, operating points recomputed per
    replicate; indicative flag applied (fewer than 30 blocks, atom above 5 %).
21. Budget check reproduced with the quasi-Poisson score rule; stratified by parent event; e-value check (R3).
22. Futility look (R2 only) taken at the registered point, if at all, and reported as non-binding.
23. Every registered number in the result entry is traced to a file and line of the frozen analysis output.
24. Number re-derived by a person (R2: ten discordant targets chosen by the seed): which, and the result.
