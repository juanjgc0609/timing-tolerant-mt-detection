# Review of the phase 0 to 7 plan (18 October 2026)

A review by an AI session of the diagnosis and plan that produced phase 0 (commits `6a24a14`, `311f4b4`). It was
checked against the code, the error register, the thesis, the results of the other work of the same day and the two
papers cited on PhaseNO. Its decisions are recorded in ADR-0004 (revised), ADR-0006 and ADR-0007 and in the roadmap.

## 1. Claims checked

| Claim | Status | Evidence |
|:--|:--|:--|
| Exact-time readers truncate the sample index; mean energy loss 22 % | Correct, narrower in scope | E52: 0.78 captured (worst 0.42) against 0.94 at the nearest sample. Only the exact-time readers of 2C and the MT-FW-EX ablation of 2D; tasks with the +-4-sample search share a shift of at most one sample |
| 3 Hz deflation kernel not updated; import-order dependence | Correct | E51; C18h to C18j at 3 Hz only |
| R1 development read 40 s of evaluation | Correct | E56; 22 against 21 of 52 from 880 s |
| Proposal not implemented (elastodynamics, adjoints, PINN) | Correct | `docs/TRACEABILITY.md`, Section 2 |
| Reproducibility incomplete | Correct | 88 files in `code/` hard-code `/home/claude/work`; no package, tests or CI; `requirements-frozen.txt` lacks numba, threadpoolctl, seisbench, quakemigrate |
| 347 Python files, 1089 JSON, 63 shell scripts | Correct | 63 shell scripts in `code/`, 70 in the repository |
| Baseline 53.3 | Correct as a record | `docs/audits/independent_audit_7.md:470`; an AI audit, not a committee |
| R1: 36 % of P and 20 % of S readings off by more than 0.5 s; 25 of 54 targets found by no method | Correct, post hoc | Thesis Section 7.2 and Section 6. The targets are the Match and Locate output, so "found by Match and Locate" is true by construction |
| PhaseNet run with peak distance 0.5 s and threshold 0.3 | Correct | `code/r0/r0_picks.py`; untuned defaults |
| Timing error 0.35 to 0.41 s as the cause to address with 3-D times | Incomplete | Homogeneous medium with truncated readings (E50); the layered medium of C18i and C18j removes much of it; the C18j gate is met exactly at threshold and fails with every event of the window held out (E53) |
| PhaseNet trained on about 700,000 arrivals | Imprecise | 623,054 P and S picks (Sun et al. 2023) |
| PhaseNO 26,176 against 21,748 events at Ridgecrest; about 94 % of the reviewed catalogue | Correct | Sun et al. 2023: 10,673 of 11,389 SCSN events, against about 88 % for PhaseNet |
| Kong et al. 2026 confirm that PhaseNO beats single-station pickers | Overstated | arXiv 2510.15281: PhaseNet has the higher P-pick F1 in all four tests; PhaseNO wins association F1 in two of four and loses one; worse with few stations; overlap not analysed |
| The commit time makes "registered, then run" checkable by a third party | Wrong | Commit dates are written by the committer; see ADR-0004 (revised) |

## 2. Problems of the plan

1. **The external anchor was not external.** Fixed by ADR-0004 (revised).
2. **Phases 1 and 2 were in the wrong order.** The scoring code applies no origin-time test under LOC although the
   thesis states one (40 % of LOC hits in 2M have origin errors above 0.5 s), and duplicates were neutral. G1-R, the
   only confirmed result, rests on that rule. Phase 2 now runs first; ADR-0007 fixes the rule for future tests.
3. **The phase 1 gate is costly.** The raw records of 2F and 2I are not stored; the frozen environment does not install
   on Python 3.13; moving `code/` would break its paths. The roadmap now decides the lock first and keeps `code/` in
   place.
4. **Straw-man risk in the central argument.** The 0.5 s merging of close arrivals is the `find_peaks` distance chosen
   for PhaseNet, not a property of PhaseNet. The comparators-at-their-best rule requires tuning it.
5. **Closest prior work missing.** EQNet (Zhu et al. 2022) already shifts and stacks learned features by travel
   times; MALMI migrates picker probabilities. The novelty to be tested is the moment-tensor consistency, the joint
   decomposition of overlapping events and the bounded-influence statistic (the register's E30 pattern).
6. **Target, metric and institution.** The audit score moved from 53.3 to 70.3 by recalibration, so it is not a
   progress metric. The work is formally an undergraduate thesis, and tying it to phase 6 could delay the degree. See
   ADR-0006. No human director is named, so the person part of the verification rule cannot yet be met.
7. **"Majority of regimes" was undefined and circular.** It is now fixed in the atlas card before the atlas runs.
8. **The schedule did not fit one person on two cores.** The atlas has a core comparator set and the theory a single
   primary line.

## 3. Items the plan omitted

- E53 and the fragility of the C18j gate.
- The missing origin-time test of LOC (being logged in the error register with the rescoring).
- The description of the generator in Section 6.2 of the thesis.
- The shared Ricker family and far-field operator between generator and detector (the mechanism part is measured by
  exploration C19).
- The failure of the deployable threshold in R1 (41.8 false detections per hour against 30).
- The two undetected R1 targets of magnitude 2 or more.
- The truth problem: no real catalogue is reliable in the overlapping regime; semi-synthetic injection with an event
  library disjoint from any deflation library.
- Power: R1 had power 0.11 to 0.13 for a gain of 0.05 (P6-power). Phase 6 needs hundreds of targets.
- The infrastructure the protocol cites (`data/cards/`, `experiments/`) did not exist; it now exists.
- A Colombian sequence (SGC) as a confirmatory candidate; the availability of an SGC FDSN service is not verified.

## Sources

- Sun, H., Ross, Z. E., Zhu, W., Azizzadenesheli, K. (2023). Phase Neural Operator for multi-station picking of
  seismic arrivals. arXiv 2305.03269.
- Kong et al. Evaluating multi-station phase picking algorithm Phase Neural Operator (PhaseNO) on local seismic
  networks. arXiv 2510.15281; Geophysical Journal International, doi:10.1093/gji/ggag242.
- Zhu, W., Tai, K. S., Mousavi, S. M., Bailis, P., Beroza, G. C. (2022). An end-to-end earthquake detection method for
  joint phase picking and association using deep learning. JGR Solid Earth.
