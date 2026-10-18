# Artifact manifest and implementation error log

Complete inventory of what was produced, with the status of each script and **the errors found**.
The root folder contains the versions whose result is reported in the proposal; the subfolder
`superseded/` contains the defective intermediates, included on purpose so that you can audit
what failed and how it was detected.

Reading rule: a script in the root supports a number that is in the proposal. A script in
`superseded/` produced a number that must **not** be used.

**Location of the files.** The paths in this manifest are relative to the root of the delivery folder
`outputs/`, not to the `docs/` folder where this file is located. The documents are in
`outputs/docs/`, and the four independent audits in `outputs/docs/auditorias/`. The scripts of the
pre-studies G to V, the drafts `Research_Proposal_..._REV2.md` to `REV12.md` and the subfolder
`superseded/` are in the previous delivery, `uploads/outputs/`. The scripts from Task 1o onwards
are in the root of `outputs/`, and the results of Tasks 2A to 2F and 2I, with a copy of the preregistration and of
its hashes, in `outputs/results/task2a/` to `outputs/results/task2f/` and `outputs/results/task2i/` (the
independent verifications of 2E, 2F and 2I in `results/task2e/verif_2e/`, `results/task2f/verif_2f/` and
`results/task2i/verif_2i/`). The exploration of seed 999 that
preceded 2F and 2I and the one that followed 2I (C9 and C10), with their scripts, outputs and reports, is in
`outputs/explore_contrib/`. In the working
folder (`/home/claude/work`) those same result files are in the root, without per-task
subfolders, except `verif_2e/`, `verif_2f/`, `verif_2i/` and `explore_contrib/`. The results of Task 2L, with
its preregistration, its hashes and its independent verification (`verif_2l/`), are in `outputs/results/task2l/`,
and the exploration C11 that preceded it, in `outputs/explore_contrib/`. The explorations C12, C13 and C14, subsequent
to 2L, are also in `outputs/explore_contrib/`, with their scripts, reports and output summaries in `out/`.

---

## 1. Documents

| File | Content |
|:--|:--|
| `Doctoral_Proposal.md` | Main proposal, titled *Timing-Tolerant, Bounded-Influence Moment-Tensor Detection of Overlapping Earthquakes: A Preregistered Benchmark and a Real-Data Test Plan*: hypotheses G1, G1-M (refuted in Task 2E), G1-R (confirmed on synthetic data in Tasks 2F and 2I) and G1-C (confirmed on synthetic data in Task 2L), design, budget, plan and Appendices A to F |
| `Supplementary_Evidence_Appendices_G_to_AA.md` | Evidence volume: 16 pre-studies (G to V), statistical audit (W), published baselines (X), controls of the catalog endpoint, Tasks 2A to 2F, 2I and 2L, the exploration that preceded 2F, the exploration C8 that preceded 2I, the explorations C9 and C10, the exploration C11 that preceded 2L and the explorations C12 to C14 that followed it (Y, up to Y.29), first real measurement (T.1), certified pruning (Z) and theory (AA) |
| `audits/independent_audit_1.md`, `_2.md`, `_3.md`, `_4.md` | The four independent audits (56, 51, 51 and 57.8 out of 100), summarised in `audit_status.md` |
| `audit_status.md` | Q1 audit of the previous version and strengthening plan, with the status after the corrections |
| `tasks_and_real_data_plan.md` | Nine tasks for your machine, closed directions, methodological warning |
| `Research_Proposal_..._REV2.md` to `REV12.md` | Revision history. REV2 to REV9 belong to the previous formulation (Bayesian association with PDE); REV10 onwards, to the constraints framework. They are kept only for traceability |

---

## 2. Scripts whose result is in the proposal

| Script | Appendix | What it measures | Status |
|:--|:--|:--|:--|
| `c1_identifiability_prestudy.py` | G | Rank and identifiability, minimum station count | valid |
| `c0_coherence.py` | H.1 | Resolution between two true sources | valid |
| `c0_clutter.py` | H.2 | Discrimination margin under clutter | valid |
| `ci_observability.py` | I.1-2 | Geometry → coercivity → empirical error | valid |
| `ci_time.py` | I.3 | Temporal half of the control condition | valid |
| `ci_lemma2.py`, `ci_lemma3.py` | I.4 | Verification of the leading-order predictor | valid |
| `ci_hetero.py` | I.5 | Heterogeneous extension with take-off slowness | valid |
| `cattrib2.py` | J | Attribution by family, ring geometry | valid |
| `cattrib3.py` | J | Idem, geometry of scattered distances | valid |
| `c3_reciprocity.py` | K.1 | Medium-independent kinematic admissibility | valid |
| `overlap4.py` | K.2 | Overlap as an association difficulty | valid (negative) |
| `a3d.py` | L | Degradation under model error, continuous optimisation | valid |
| `a2_minimal.py` | M | Minimal sufficient set, A1 and A2 | valid |
| `a2b.py` | M | Stress with coda-like clutter | valid, with caveat |
| `a2c.py` | M | Adversarial stress with decoy | valid |
| `a4_b2.py` | N | Cost per evaluation (A4) | valid for A4 |
| `b2b.py` | N | Network design by observability | valid, saturated endpoint |
| `f1b.py` | O | Sufficiency ladder, gate F1 | valid |
| `c3_single.py` | P | One station, several events | valid |
| `b1_interfaces.py` | Q | Certificate with interfaces and head wave | valid |
| `d1_coupling.py` | R | Detection-association coupling | valid |
| `b3b.py` | S | Observability maps where discriminating is difficult | valid |
| `e1_critical.py` | T | Nearby coherent clutter | valid |
| `a_i_b.py` | U | Association without a velocity model, untuned | valid |
| `c1c.py` | V.1 | Inter-station rank per frequency bin | valid (negative) |
| `b1_moment.py` | V.3 | Physicality of the moment tensor | valid (negative) |
| `octree_hypothesis_space.py` | V.4 | Octree over the hypothesis space | valid (negative) |
| `a_i_tuned_attempt.py` | U | Tuned azimuth estimator | **DOES NOT WORK**, see §4 |
| `task1_fetch.py`, `task1_analyze.py`, `task1_decide.py` | T.1 | Polarisation coherence of real clutter, CM network | valid |
| `diagnose_polarization_config.py` | T.1 | Per-station control: band, window, rotation, reflection, stability, magnitude and depth | valid |
| `task1_verify.py` | T.1 | Polarisation estimator on synthetic data with known back-azimuth | valid |
| `task1c_network_survey.py` | Task 1c | Map of polarisation reliability across the network, event windows only | tested on synthetic data, not run on real data |
| `task1e_polarization_reliability.py` | Task 1e | G1 under the measured reliability and per-station gate | valid |
| `task1_clutter_robustness.py` | Task 1f | Clutter verdict by subpopulation | tested on synthetic data, not run on real data |
| `task1jk_controls.py` | Y.1, Y.3 | Control of ties (1j) and empty windows with absolute threshold (1k); GaMMA and PyOcto tuned in calibration | valid |
| `task1j2_noise.py` | Y.2 | Polarisation versus the best continuous kinematic score, with pick-time and velocity error | valid, exploratory |
| `PREREGISTRO_1j_1k.txt` | Y | Decision rules written before each run, with dated amendments | registration |
| `diag_pol.py` | Y.4 | Discriminability of the synthetic polarisation channel (AUC, azimuthal hit) | valid |
| `task1o_consistent.py` | Y.5 | Consistent generator (detection and time according to signal-to-noise), KIN-LR, PHYS-LR (9d), WAVE-E, WAVE-P | generator valid; the reading of its endpoint was corrected in Y.6 (E21) |
| `task1o_baselines.py` | Y.5 | GaMMA and PyOcto tuned on the same regenerated windows | replaced by `task1r_baselines.py` (E21, E22) |
| `analyze_1o.py` | Y.5 | Table Y.5 | valid |
| `task1r_fair.py`, `analyze_1r.py` | Y.6 | Fair repetition of 1o: FROC at realised false events, negative control with permuted cosines, exact condition, CI by paired bootstrap | valid |
| `task1r_baselines.py` | Y.6 | GaMMA and PyOcto with the domain of the in-house methods, minimum picks from 3, oracle on the evaluation, GaMMA with fixed seed | valid |
| `check_gamma_det.py` | Y.6 | Variability of GaMMA without seed and identity of the two code paths with seed | valid |
| `task1s_gamma_pol.py`, `analyze_1s.py` | Y.7 | Polarisation likelihood as a filter on the GaMMA events, with negative control | valid |
| `task1t_wave.py`, `analyze_1t.py` | Y.8 | Pick-free waveform migration: energy, projection onto the predicted polarisation and control with the direction of another station | valid |
| `task1u_hetero.py`, `analyze_1u.py` | Y.9 | Unmodelled medium: per-station delays and polarisation deviation that the homogeneous model does not predict | valid |
| `froc_exact.py` | Y.6 to Y.9 | Exact FROC curve from the greedy list of each window (E23) | valid, verified against the original emission |
| `analyze_1r200.py`, `make_tables_y.py` | Y.6 to Y.9 | Rules with 200 windows and appendix tables generated from the JSON files | valid |
| `task1w_mech.py`, `analyze_1w.py` | Y.10 | Double-couple sources; controls that take the geometry of another node; WAVE-PT; replication with new windows | valid |
| `chk_gamma_loc.py`, `chk_gamma_tr.py` | Y.7 | Picks per source, GaMMA failures by distance to the nearest event, effect of the origin-time filter | valid; measured with the clutter of E26 |
| `task1x_mtglrt.py`, `analyze_1x.py`, `cal_capv.py`, `chk_nms.py` | Y.11 | Coherent moment-tensor detector (dual certificate of the group BLASSO), Frank-Wolfe with deflation, calibration of the per-station cap, sensitivity to the suppression radius | valid; the GaMMA column uses the clutter of E26 |
| `task1y_clutter.py`, `analyze_1y.py` | Y.12 | Same comparison with consistent clutter (weak and strong), control of the complete method with the geometry of another node | valid |
| `bnb_octree.py`, `bnb_scaling.py` | Z | Certified pruning of the hypothesis space and its scaling | valid |
| `gen_realistic.py` | Tasks 1p, 2A | Generator mismatched with the detectors: linear-gradient medium with curved rays, Ricker source function with variable central frequency and phase, scattered coda; with g = 0 and VS0 = 3 km/s it reproduces the times and amplitudes of `generate_y` | valid; self-test included |
| `real_picker.py` | Tasks 1p, 2A | PhaseNet (original 2019 model, not retrained) and classical STA/LTA on the waveforms, calibrated only on seed 999, against `degrade_y` on the same arrivals | valid; PhaseNet outside its domain (E27) |
| `qm_wrapper.py` | Tasks 1v, 2A | QuakeMigrate 1.2.2 end to end through its public API (archive, travel-time table, detect, trigger); variants QM (P on the vertical, S on the horizontals) and QM-3C (three components), calibrated on seed 999 | valid, exploratory |
| `station_corrections.py` | Tasks 1z, 2A | Per-station and per-phase time corrections from catalogued events; residual and recall of MT-FW without correction, with K = 5, 10 and 20 events and with the exact correction | valid |
| `rft_threshold.py` | Task 1z | Threshold of the certificate maximum: EEC, Bonferroni, GEV, Gumbel and conformal against the realised rate in 7,200 noise-only and clutter windows | valid; withdraws the guarantee of Section 7.5 (E28) |
| `certified_fw.py` | Tasks 1z, 2A | Certified oracle of MT-FW: Lipschitz energy bound, box variant and per-subspace envelope; records identical to `fw_record` | valid |
| `chk_aa_bound.py` | AA | Proposition AA.3 (cell bound) on the implemented statistic, raw and deflated residual | valid |
| `chk_aa_ec.py` | AA | Euler characteristic densities of chi-squared fields, symbolic derivation against Worsley | valid |
| `task2a_confirm.py`, `analyze_2a.py` | Task 2A, Y.18 | Preregistered confirmatory benchmark: mismatched generator, STA/LTA picks, station corrections, GaMMA, PyOcto and QuakeMigrate against MT-FW+SC; rules 2A-a to 2A-e | executed and verified; results in `results/task2a/`, independent verification in `results/task2a/task2a_report.txt` and its scripts in `results/task2a/verif_2a/` |
| `task2b_ablation.py`, `analyze_2b.py` | Task 2B, Y.19 | Ablation of the radiation pattern: MT-FW-FREE+SC, with the radiation operator replaced by the identity inside the same solver, on the windows, networks and corrections of 2A regenerated bit for bit; exact control against the records of 2A; rule 2B | executed; holds in both regimes (weak at the minimum, 2 of 3 cells per seed); results, controls and incident E31 in `results/task2b/` and `results/task2b/task2b_report.txt`; `analyze_2b.py` was edited 18 s after its first output to add one exploratory line per regime, without changing rules, tables or verdict (E36) |
| `task2b_kernel_check.py` | Task 2B, Y.19 | Recomputes with the default OpenBLAS kernel the 18 networks of 2A written with the SkylakeX kernel and compares records and recall | valid, exploratory; the kernel does not change any recall (E31) |
| `task2c_methods.py` | Task 2C, Y.20 | BOMP+SC (block orthogonal matching pursuit in the form of Vera Rodriguez et al., 2012), MTSTACK+SC (mechanism-corrected stack in the form of Anikiev et al., 2014) and GL-FW+SC (group BLASSO on the grid by fully corrective Frank-Wolfe with duality gap), in the matched-filter domain and at the exact predicted times, with the same dictionary of moment-tensor blocks | valid; inherited from a previous session, reviewed and corrected before running (E32) |
| `task2c_run.py` | Task 2C, Y.20 | Runner: validation, smoke test and runs of part (i) (windows of 2A regenerated bit for bit, with exact control per network) and of part (ii) (25 new networks per cell, seed and regime with MT-FW+SC, QM-3C+SC, BOMP+SC and MTSTACK+SC); sets, checks and records the OpenBLAS kernel of each network | executed; inherited and corrected before running (E32); results in `results/task2c/` |
| `task2c_validate.py` | Task 2C, Y.20 | Validations V1 to V4 on seed 999: dictionary identities against explicit atoms, recovery of a noise-free source, KKT and gap of GL-FW, MTSTACK equal to the projection energy of WAVE-MT | valid; written because the inherited runner imported it and it did not exist (E32) |
| `analyze_2c.py` | Task 2C, Y.20 | Rules 2C-a, 2C-b and 2C-c; bootstrap by windows, by network and hierarchical; Holm over 36 comparisons per regime; tables | executed; 2C-a and 2C-c hold in both regimes; 2C-b: GL-FW+SC below MT-FW+SC in 23 of 24 cases (E33) |
| `task2c_explore.py` | Task 2C, Y.20 | Diagnosis of the collapse of the exact-time methods: ideal generator, sources on the nodes, oracle corrections | valid, exploratory and written after seeing part (i); does not support any verdict |
| `task2d_methods.py` | Task 2D, Y.21 | Reading at the energy peak within ±0.25 s ("PK", the `node_vectors` rule of MT-FW) and the five methods: BOMP-PK+SC, MTSTACK-PK+SC and GL-FW-PK+SC (those of 2C read at the peak; GL-FW-PK with a growing dictionary and gap with respect to the realised dictionary), MT-FW-EX+SC (without peak search) and MT-FW-NOCAP+SC (without cap) | valid; hash recorded at 02:41:59 UTC before any result (`task2d_code_hashes.txt`) |
| `task2d_run.py` | Task 2D, Y.21 | Runner: validation, smoke test, part (i) (windows of 2A regenerated bit for bit, with exact control per network and fixed OpenBLAS kernel) and part (ii) (200 empty windows of seed 999 per cell and regime, maxima of MT-FW+SC and QM-3C+SC, reproduction control of 2A) | executed; hash prior to the results; results in `results/task2d/` |
| `task2d_validate.py` | Task 2D, Y.21 | Validations V1 to V6 on seed 999: tabulated picks against `node_vectors_c`, exact reduction to 2C without peak, explicit atoms with shifted times, independent paths, KKT and gap of GL-FW-PK, injected time error (V6, descriptive) | valid; V1 to V5 correct |
| `analyze_2d.py` | Task 2D, Y.21 | Rules 2D-a, 2D-b and 2D-c with their readings fixed in the docstring; hierarchical bootstrap (networks and windows) and by windows; Holm over 36 comparisons per regime; tables and diagnostics | executed; 2D-a does not hold in any regime (E35), 2D-b holds, 2D-c holds per cell and not in the strict one with weak clutter; one diagnostic line changed after the part files existed, without touching rules or tables (hash before 110aa893…73a7, after 8dd5013b…d2b4) |
| `task2d_explore.py`, `task2d_explore_sum.py` | Task 2D, Y.21 | E1: where the false detections of MT-FW+SC fall with the conformal thresholds in windows with events; E2: V6 repeated with 90 windows per cell (injected time error of ±0.2 s) | valid, exploratory and written after seeing the results; does not support any verdict |
| `task2e_methods.py` | Task 2E, Y.22 | Continuous atom model (position, origin and one shift per station and phase), exact band-limited correlation, continuous Gram, oracle (statistic of MT-FW-NOCAP), refinement and sliding by L-BFGS-B with analytical gradients, convex step with gap; TS-SFW+SC and BOMP-PK-SLIDE+SC (`tssfw_record`), `mtfw_sc_locs` | valid; hash recorded at 15:38:31 UTC before the evaluation (`task2e_code_hashes.txt`); the independent verification found no implementation error |
| `task2e_run.py` | Task 2E, Y.22 | Runner: validation, calibration of f on seed 999, time estimate and evaluation on 30 networks per cell (windows of 2A and 2C regenerated bit for bit, OpenBLAS kernel fixed and checked); refuses to evaluate if the registration hash is not the latest and does not contain the note with the frozen f | executed; results in `results/task2e/` |
| `task2e_validate.py` | Task 2E, Y.22 | Validations V1 to V5 on seed 999: continuous identities, gradients against finite differences, convex step against an explicit dictionary and FISTA, recovery of a noise-free source (V4), reduction to BOMP-PK | V1, V2, V3, V5 and the basin of V4 correct; the criterion of V4 from the grid does not hold (12 of 90), declared as a property of the method |
| `analyze_2e.py` | Task 2E, Y.22 | Rules 2E-a (with Holm) and 2E-b, exact FROC, hierarchical bootstrap, recall tables, descriptive differences, location error and convergence diagnostics | executed; 2E-a and 2E-b do not hold in any regime (G1-M refuted); tested before the evaluation only on the calibration parts, without changes afterwards |
| `task2e_explore.py` | Task 2E, Y.22 | Displacement of the atoms by the sliding, changes of hit and actual times | exploratory, written after the results; does not support any verdict |
| `explore_contrib/xc.py`, `explore_contrib/analyze.py`, `explore_contrib/t_check.py` | Y.23 | Common harness of the exploration (seed 999 only, networks `task2a_confirm.network(999, celda, j)` with new j, streams with TAG ≥ 3000 not used before), recall with hierarchical bootstrap and paired differences, check that it reproduces `setup_net` | valid, exploratory |
| `explore_contrib/c1_run.py`, `c1_analyze.py` | Y.23 (C1) | Grids of 3, 6 and 9 km with BOMP-PK+SC, TS-SFW+SC and MT-FW+SC | exploratory; null |
| `explore_contrib/c2_run.py`, `c2_analyze.py`, `c2_ablate.py`, `c2_ablate_analyze.py` | Y.23 (C2), Task 2F | Dense windows (3 to 5 events at < 9 km, `draw_dense`) and ablation of the cap, of the P/S cross term and of the deflation (MT-FW-NOCAP, BOMP-PK-CAP, BOMP-PK-SEQ); `c2_run.draw_dense` and `c2_ablate.PKDictCap` were frozen for 2F | exploratory; its candidate was confirmed in 2F; the mechanism reading "the cap resolves dense sequences" was too broad (E37) |
| `explore_contrib/c3_run.py`, `c3_analyze.py` | Y.23 (C3) | Statics of 0.3 to 0.8 s with and without catalog corrections, velocity 7.5 % slower, wider peak search, shared shifts and online station terms | exploratory; null with corrections |
| `explore_contrib/c4_run.py`, `c4_analyze.py`, `c5_run.py`, `c5_analyze.py` | — | False alarms in windows with events (C4) and location accuracy and uncertainty (C5) | written and never run; `c5_run.py` uses TAG 3006, which C6 already used: change it before running it |
| `explore_contrib/leak_test.py`, `leak_analyze.py`, `leak2.py`, `leak2_analyze.py` | Y.23 (masking) | Masking of arrivals inside BOMP-PK and MT-FW: pilot with 12 variants and replication with the variant fixed on networks 30-35 | exploratory; doubtful, weak clutter only |
| `explore_contrib/c6_lib.py`, `c6_validate.py`, `c6_run.py`, `c6_ceil.py`, `c6_analyze.py`, `c6_percell.py`, `c6_ttresid.py`, `c6b_run.py`, `c6b_analyze.py` | Y.23 (C6, C6b), Task 2F | 1D models (homogeneous, oracle gradient, calibrated by JHD, two layers) inside BOMP-PK and MT-FW; ceiling CEIL-GD; attribution run that removes one generator mismatch at a time; `c6_lib.ModelOps` and `c6_lib.mtfw_dense` were frozen for 2F | exploratory; null (the margin persists with the exact medium); `c6_lib` validated against the registered detectors (≤ 3e-13) |
| `explore_contrib/c7_lib.py`, `c7_validate.py`, `c7_resolution.py`, `c7_run.py`, `c7_analyze.py`, `c7_power.py` | Y.23 (C7) | Subspaces invariant to the source pulse (PH, PHF), pulses outside the family (Gabor, Berlage), resolution diagnostic, power | exploratory; negative; by-product: without the P/S cross term BOMP-PK is reproduced in 2400 of 2400 windows |
| `explore_contrib/power.py` | Y.23, Task 2F | Power of the confirmatory designs from the pilot SD of the paired difference | valid; reproduces the pilot SDs of the registration (0.0383 and 0.0353) |
| `task2f_run.py` | Task 2F, Y.24 | Runner: controls on seed 999 and run of 2F on seeds 4242 and 8888 (dense arm j = 0..95, TAG 4201; standard arm j = 96..111, TAG 4202), with MT-FW+SC, BOMP-PK+SC, BOMP-PK-CAP+SC and MT-FW-NOCAP+SC; Prescott kernel checked per process | executed only once after the freeze (hashes of 26 modules at 03:49:08 UTC, `task2f_code_hashes.txt`); results in `results/task2f/` |
| `task2f_control.py` | Task 2F, Y.24 | Controls on seed 999: bit-for-bit reproduction of the dense pilot and of standard records, identity MT-FW-NOCAP = BOMP-PK-SEQ, fast FROC against the reference one, time estimate | valid; all controls correct (`task2f_control.txt`) |
| `analyze_2f.py` | Task 2F, Y.24 | Rules 2F-a (with Holm), 2F-b and 2F-c (TOST), exact FROC with one threshold per method over the pooled windows, hierarchical bootstrap, per-cell tables and times | executed once, frozen before data from 4242/8888 existed; 2F-a, 2F-b and 2F-c hold; the hit rule without origin time is lax in dense windows (E38) |
| `explore_contrib/c8_lib.py` | Y.25 (C8) | Copy of `gen_realistic.generate_real` that returns the amplitude factors and allows them to be fixed; IID worlds (factors drawn in each window, the generator unchanged) and PERS worlds (factors fixed per network, shared by catalog and windows); catalog amplitude terms per station and phase; MT-FW engine with interchangeable statistic (caps 3, 6, 8, 12 and no cap, ORACLE, SITE-CAT, HUBER-N, HUBER-S, GAIN, GAIN-FREE, NORM-TR) | valid, exploratory; seed 999 only, networks 60-79 and 80-89, TAG 3801 to 3810 |
| `explore_contrib/c8_sitedev.py` | Y.25 (C8) | Choice of the estimator of catalog amplitude terms against the true factors on networks 60-63 (correlation 0.42 with the per-event ratio, 0.74 with the radiation-averaged one) | exploratory; chosen without detection results, on networks that are also evaluated (slightly favours SITE-CAT in PERS) |
| `explore_contrib/c8_validate.py` | Y.25 (C8) | Validations V1 to V8: copy of the generator equal to the original bit for bit, catalog and corrections equal to `xc.NetCfg`, engine equal to `c6_lib.mtfw_dense` with and without cap (6 of 6), Huber and GAIN limits equal to least squares | valid; all correct (`out/c8_validation.txt`); convergence of GAIN-FREE not validated |
| `explore_contrib/c8_run.py`, `c8_analyze.py` | Y.25 (C8) | Dense and standard windows in the two worlds with dispersion 0.5, sensitivity to the dispersion (0.25 and 0.8) and later sweep of the Huber threshold; recall with hierarchical bootstrap of 1000 replicates, paired differences and fraction of the margin up to CEIL-GD | exploratory; site terms do not replace the cap (negative); HUBER-N matches cap 6 with dispersion up to 0.5, is worse with 0.8 and costs 9 times more; best of 11 candidates (selection bias); E39 arises from here |
| `task2i_run.py` | Task 2I, Y.26 | Resumable runner (one part per network): controls on seed 999 and run of 2I on seeds 5151 and 6262 (primary arm j = 0..119, 5 dense windows with TAG 5201 and 5 standard ones with TAG 5202 per network), with BOMP-PK+SC, MT-FW+SC with cap 6, MT-FW-NOCAP+SC and HUBER-N+SC (`c8_lib.mtfw_general`, threshold 2.7955); secondary arms `thr` (thresholds 2.5003, 3.3682 and 5) and `s25`/`s80` (site deviation 0.25 and 0.8), j = 0..19; Prescott kernel checked per process | executed and verified; secondary arms analysed (descriptive, without independent verification; parts written between 22:04:55 UTC on 2 October and 03:33:41 UTC on 3 October). Preregistered (`PREREGISTRO_1j_1k.txt`, section Task 2I, 2 October 2026, 18:00:05 UTC, hash a983beac…6f97); frozen at 18:14:44 UTC (hashes of 27 files, `task2i_code_hashes.txt`); run interrupted twice and resumed from the parts with the code unchanged (`task2i_freeze_time.txt`); results in `results/task2i/` |
| `task2i_control.py` | Task 2I, Y.26 | Controls on seed 999: exact reproduction of the C8 pilot (4 primary methods, secondary thresholds, deviations 0.25 and 0.8), `c8_lib.mtfw_general('ls')` against `c6_lib.mtfw_dense`, FROC of the analysis against `xc.recall_caps` and the published values, bootstrap of the analysis against that of `REPORT_C8.txt`, time estimate and budget rule | valid; all controls correct (`task2i_control.txt`) |
| `analyze_2i.py` | Task 2I, Y.26 | Rules 2I-a (superiority in dense windows), 2I-b (non-inferiority in dense windows, margin 0.03) and 2I-c (TOST equivalence in standard windows, margin 0.05); exact FROC with one operating point per method; hierarchical bootstrap of 2000 replicates (strata regime × cell, networks and windows); contrasts by regime and cell, descriptive replication of 2F, hit rule with origin time and times | executed on the complete primary arm (`analyze_2i_primary_out.txt`, where the sections of the secondary arms say "not run") and, without changes to the code, on the complete task (`analyze_2i_out.txt`, paragraph of the registration of 3 October 2026, 05:54:02 UTC); frozen before data from 5151/6262 existed; 2I-a, 2I-b and 2I-c hold on both seeds; secondary arms, descriptive: thresholds 2.5003 and 3.3682 indistinguishable from the registered one, threshold 5 lower; with site deviation 0.25 the cap contributes nothing; with 0.8 HUBER-N falls below the cap (−0.055 and −0.069) |
| `explore_contrib/c9_lib.py` | Y.27 (C9) | Travel-time-free location statistics on the P and S peak vectors per station (azimuth intersection, P direction, transverse S, projection onto the moment-tensor subspace with free times, 1/r amplitudes, fit with free gain per station); vector readers (at the true arrival, blind, width of ±1.5 s); BOMP-PK and MT-FW loops with relocation and selection hooks (SELPOL, PENPOL, PEN) | valid, exploratory; seed 999 only, networks 130-149 and 199 (validation), TAG 3901 to 3913; operators equal to the registered ones (≤ 3e-16) and loops equal to `xc.bomp_pk` and to `c8_lib.mtfw_general` in 9 of 9 windows (`out/c9_validation.txt`, `c9_validation_end.txt`) |
| `explore_contrib/c9_standalone.py`, `c9_standalone_analyze.py` | Y.27 (C9) | Location of 600 isolated events (20 networks × 10 events × 3 cells, weak clutter) with each time-free statistic, under 7 paired generator configurations (polarisation deviation and site deviation); noise-free validation | exploratory; median error of 7.7-8.4 km with the best channel against 0.5-0.7 km with the times; the blind reader finds the P peak at 19-20 % of the stations |
| `explore_contrib/c9_fisher.py` | Y.27 (C9) | Bayesian Cramér-Rao bound (Van Trees) of the epicentre per channel (times, P and S polarisation, amplitudes, S/P ratio, combinations) with the benchmark prior; 200 sources per cell; Monte Carlo check on network 130 | exploratory; times 0.4-0.7 km, time-free channels together 4.4-5.6 km, joint improvement of 4-5 %; maximum-likelihood estimators 1.2-1.4 times above the bound; local bound, homogeneous rays, without the polarity information |
| `explore_contrib/c9_run.py`, `c9_extra.py`, `c9_analyze.py`, `c9_pairs.py` | Y.27 (C9) | Fusion with BOMP-PK+SC and MT-FW+SC in dense and standard windows: reordering, relocation (after and inside the loop), extra term; oracle worlds (without polarisation deviation, without amplitude factors, both) and worlds with 20 and 45 degrees; variants SELPOL and PENPOL in `c9_extra.py`; recall with hierarchical bootstrap of 1000 replicates (500 in the oracle worlds) and ceiling CEIL-GD | exploratory; relocation loses recall in all variants, oracles included (negative); SELPOL, added after seeing results and best of about 50 candidates, +0.046 [+0.007, +0.070] over MT-FW in dense windows, null at 20-45 degrees; without multiplicity correction; `c9_extra.py` did not run at 0 or at 20 degrees |
| `explore_contrib/c10_lib.py` | Y.27 (C10) | Rule for the time level from the catalog (pooled robust deviation, floor 0.05 s); penalised peak choice; MT-FW and BOMP-PK loops with Gaussian time loading (T1, T1c, T2, T2w, T4, T4h), local refinement on a 1 km grid (T3) and oracle of true statics; hit rules with and without origin time | valid, exploratory; seed 999 only, networks 160-179 and 195 (validation and choice of the level rule), TAG 4001 to 4004; in hard mode it reproduces the registered detectors in 9 of 9 windows (`out/c10_validation.txt`); the first version, with the energy clipped in the hard rule, reproduced MT-FW in 1 of 9 and was corrected before the runs |
| `explore_contrib/c10_run.py`, `c10_analyze.py` | Y.27 (C10) | Main runs in dense and standard windows (15 variants and CEIL-GD), recall with the two hit rules, hierarchical bootstrap of 1000 replicates (500 per regime), location and near-hit diagnostics, fraction of the margin | exploratory; no variant that loads the shifts gains recall with catalog corrections; widening the window to ±0.5 s loses; true statics +0.047 (dense) and +0.098 (standard), with a prior at the level of the model error +0.081 and +0.146; T3 only with the origin-time rule in standard windows; set of runs trimmed to the compute budget |
| `explore_contrib/c10_sweep.py`, `c10_extra.py`, `c10_diag.py`, `c10_diag_analyze.py` | Y.27 (C10) | Sweep of fixed levels (0.05, 0.10 and 0.20 s; networks 160-169); true statics with the catalog level; separability in the first iteration without deflation (standard windows, networks 160-169) | exploratory; `c10_extra.py` and `c10_diag.py` written after the results; no fixed level gives a gain; the prior only pays off with true statics and a small level; with the oracle 37-39 % of the isolated events still lose in the first iteration |
| `explore_contrib/c10_screloc.py`, `c10_screloc_run.py` | Y.27 (C10) | Station corrections re-estimated with the catalog events relocated by their own times (SC-R) and detection with them | exploratory, written after the results; the rule was written after seeing the oracle and its effect was checked against the true statics before the detection run; correction error from 0.145 to 0.110 s (weak clutter) and from 0.214 to 0.202 s (strong); on its own it does not change the recall |
| `task2l_run.py` | Task 2L, Y.28 | Resumable runner (one part per network): controls on seed 999 and run of 2L on seeds 7171 and 9393 (networks j = 0..159 with 5 standard windows and j = 0..99 with 5 dense ones; MT-FW with cap 6 with S0, with S0 and T3 and with S2 and T3; BOMP-PK on the standard windows of j = 0..59; oracle of true delays with T3 on j = 0..29) | executed and verified. Preregistered (`PREREGISTRO_1j_1k.txt`, section Task 2L, 4 October 2026, 21:10:33 UTC, hash 13a6e9b1…a62a); frozen at 21:24:48 UTC (hashes of 37 files, `task2l_code_hashes.txt`); 1920 parts; run interrupted when the session was suspended and resumed from the parts with the code unchanged (`task2l_freeze_time.txt`); results in `results/task2l/` |
| `task2l_control.py` | Task 2L, Y.28 | Controls on seed 999: exact reproduction of the C11 pilot in 18 parts, S0 equal to the registered correction, engine corrections equal to those saved by 2A, seed input, bootstrap of the analysis against the published one, time estimate and budget rule | valid; all controls correct (`task2l_control.txt`); estimate of 4.62 h, budget rule not triggered |
| `analyze_2l.py` | Task 2L, Y.28 | Rules 2L-a to 2L-d and secondary on BOMP-PK; exact FROC with one operating point per method; hierarchical bootstrap of 2000 replicates (strata regime × cell, networks and windows); descriptive contrasts by regime and cell, oracle arm, errors of the terms and of the hits | executed; frozen before data from 7171/9393 existed; the four rules and the secondary hold on both seeds; the verification re-ran it with byte-for-byte identical output |
| `explore_contrib/c11_lib.py`, `c11_run.py`, `c11_dev.py`, `c11_validate.py`, `c11_analyze.py`, `c11_posthoc.py`, `c11_power.py` | Y.28 (C11) | Reading of energy peaks, pick-free joint determination of hypocentres and station terms (S2), self-calibration with the detector atoms (S3), errors of the terms, development, validation, detection with each set of terms and analysis | exploratory; seed 999, networks 200-209 (development), 210-239 (evaluation) and 195 (validation); plan frozen with hash before any detection result (`C11_PLAN.txt`, 18:01:51 UTC on 4 October); primary contrast +0.070 [+0.010, +0.120]; `c11_analyze.py` was changed once after the freeze (the frozen file is kept as `c11_analyze_frozen.py.txt`) |
| `explore_contrib/c12_lib.py`, `c12_run.py`, `c12_analyze.py`, `c12_validate.py`, `c12_mech.py` | Y.29 (C12) | Full, deviatoric and double-couple statistics inside the MT-FW loop; optimiser over the double-couple manifold; constants and tail of the tube formula; null law and detection in standard and dense windows; mechanism (written after the results) | exploratory; seed 999, networks 1200-1229 and 1290; validations V0 to V6 correct (`out/c12_validation.txt`); the null law of the tube formula agrees with the simulation; null for detection (double couple − full +0.005 [−0.021, +0.033] in standard windows with origin time) |
| `explore_contrib/c13_lib.py`, `c13_run.py`, `c13_run_d3c.py`, `c13_twin.py`, `c13_analyze.py`, `c13_tables.py`, `c13_report.py`, `c13_report_text.py` | Y.29 (C13) | Decoys D1p, D1o, D2s, D2r, D3 and D3c; calibrations tdc, conformal p-value with Benjamini-Hochberg, rank competition and fpw; event-free twins; hierarchical bootstrap of 400 replicates | exploratory; seed 999, networks 1400-1411; D1p controls with the epicentre-only rule and not with origin time; D3c, fpw and the twins were added after intermediate looks, declared; first version of the twins invalid (E43); D3 and D1o with a construction defect (E44) |
| `explore_contrib/c14_lib.py`, `c14_run.py`, `c14_analyze.py`, `c14_frozen.json`, `c14_dev2.json`, `C14_FREEZE.txt` | Y.29 (C14) | MT-FW loop with remedies R1 to R4, oracle deflation and exact decomposition of the statistic; per-event parts and twins; diagnosis, development, freeze and evaluation | exploratory; seed 999, networks 1600-1602 (development) and 1610-1621 (evaluation, generated after freezing); with empty configuration it reproduces `c10_lib.mt_loop` atom by atom; negative for recall; incidents and changes around the freeze in E45 |
| `t2m_data.py`, `t2m_stations.py`, `t2m_prep.sh`, `t2m_validate.py`, `t2m_power.py` | Task 2M, Y.31 | Real-noise preparation of any PyMPA37 day with the quality rule qc1, station rule, validation against C16 (noise bit for bit, slot plans) and power analysis | valid; the noise of the C16 day regenerated bit for bit in the rebuilt environment; C16 parts not preserved (E49) |
| `t2m_lib.py`, `t2m_run.py`, `t2m_queue.sh`, `t2m_freeze.sh`, `t2m_resume.sh`, `t2m_resume_each.sh` | Task 2M, Y.31 | Injection of generator events into real noise, five methods, part files per (arm, cell, network); frozen run with SHA-256, resumption per network | valid; a declared per-network failure stopped a lane (E48), resumed with the frozen code |
| `t2m_analyze.py`, `results/task2M/posthoc/blocks.py`, `results/task2M/verif_2m/*_verif.py` | Task 2M, Y.31 | Registered analysis (primary seed 2718, Holm over two tests), post hoc breakdown by blocks of hours, independent verification | valid; the breakdown is exploratory; verification confirms with caveats of scope |
| `explore_contrib/c18_tt.py`, `c18a_timing.py`, `c18a_overlap.py`, `c18b_oracle.py`, `c18c_oracle.py`, `c18d_detect.py`, `c18d_analyze.py`, `c18e_run.py`, `c18e_analyze.py`, `c18f_run.py`, `c18f_analyze.py`, `c18g_oracle.py`, `c18h_kernel.py`, `c18h_detect.py` | Y.32 (C18) | Timing ladder with layered media on the Ridgecrest CAL readings, oracle scores at the true hypocentre for several statistics, full detector on CAL with those statistics, robustness of the statistics to real-noise transients on seed 999 | exploratory; Ridgecrest CAL only, EVAL not read; seed 999, networks 3100-3107, TAG 5181-5186 |
| `explore_contrib/c18i_widesc.py`, `c18j_layered.py` (and `c18a_timing.py` with SEARCH 2.0) | Y.32 (C18i, C18j) | Readings and station terms with a +-2 s search; detector with the layered medium of Hadley and Kanamori; the gate check on CAL | exploratory; Ridgecrest CAL only; withdraws the C18a conclusion (E50) |
| `r0/queue_close.sh`, `r0/run_close/` | Y.30 (closure) | M1 and M3 on EVAL of the four subnetworks with the frozen R0 code, re-run after the original outputs were lost (E49), and the registered secondary contrast M1 - M3 | post hoc, declared in the log before running |

---

## 3. The design errors that inverted or nullified a result

This is the central content of this manifest. The log contains 56 errors, E1 to E56. In the
first six (E1 to E6) the resulting number looked valid, and in all six it was detected **only because
it turned out to be implausible**, not because there was a visible fault in the code; the following ones were detected
by controls, preregistered rules or audits, as indicated in each one.

### E1. Exclusion radius larger than the measured separation
**Script:** `superseded/c1_certificate.py`
**Symptom:** the dual certificate was satisfied at every separation, down to 20 m, in all
geometries.
**Cause:** the exclusion radius when searching for the maximum of |η| outside the support was 0.08 km while
the tested separation went down to 0.02 km, so the region between the two peaks was excluded from the
maximum.
**Detection:** the result was physically absurd.
**Correction:** `superseded/c1_v2.py`, global-maximum criterion with tolerance.
**Note:** after correcting it, the noise-free certificate is indeed satisfied at every separation —
the original result was correct by accident, but for the wrong reason.

### E2. Metric quantised by the grid
**Scripts:** `superseded/a3b.py`, `superseded/a3c.py`
**Symptom:** all strategies gave 0.00 or exactly one cell of error, with no way to
tell them apart.
**Cause:** a 4.5 km grid measuring sub-kilometre differences.
**Correction:** `a3d.py`, continuous Nelder-Mead optimisation from the grid maximum.

### E3. Constraints turned into no-ops
**Script:** `superseded/cattrib.py`
**Symptom:** C-I and C-II gave figures **identical** to the condition with all constraints
active, and 1137 survivors in all cases.
**Cause:** two independent ones. The observability floor was at `GMIN=0.0`, so it never pruned.
And the silence test compared predicted amplitude against observed energy on incompatible scales,
so it never fired.
**Detection:** the exact identity of the figures.
**Correction:** `cattrib2.py`, with the floor at the 20th percentile and a rescaled silence test.
**Lesson:** a null result produced by an inactive constraint is not a null result. The
design-validity step of Section 13.0 comes from here.

### E4. Undersampling of the sensitive parameter
**Script:** `f1_c1.py` (kept in the root because it is cited in Appendix O)
**Symptom:** the matched filter came out catastrophically bad, rates of 0.80 to 1.00 over the whole range,
and **it would have closed the entire waveform branch**.
**Cause:** the correlation peak of a matched filter is half a period wide; the origin time was
sampled in steps of a quarter period, so the true source fell between samples while
some false candidate landed on the peak.
**Correction:** `f1b.py`, sampling at one tenth of a period. **The conclusion was completely inverted.**

### E5. Failure value coinciding with the correct answer
**Scripts:** `superseded/a_i_tuned.py`, `superseded/a_i_lateral.py`
**Symptom:** the tuned estimator gave **an error of 0.00 km at every level of model error**, even
with 20° of azimuthal bias.
**Cause:** `locate()` returned `np.array([0.,0.])` when RANSAC found no solution, and the true
source was at the origin. **Total failure was counted as a perfect hit.**
**Detection:** the 0.00 with 20° of bias was physically impossible.
**Correction:** `superseded/a_i_fixed.py`, source away from the origin and failure rate reported
explicitly. The estimator fails in 62-100% of the cases.

### E6. Averaging that destroys the structure under test
**Script:** `c1_csdm.py` (root, cited as the first attempt)
**Symptom:** the rank-based source count gave 3 for K=2 and 6 for K=3.
**Cause:** the cross-spectral matrix was averaged **over a frequency band**. The steering vector
depends on ω, so averaging over a band destroys the rank-one structure that was
meant to be measured.
**Correction:** `c1c.py`, one matrix per frequency bin. The count still fails, but now
for a real physical reason and not because of the averaging.

---

### E7. The true source falls exactly on a grid node
**Affected scripts:** eighteen, among them `a2_minimal.py`, `a2b.py`, `a2c.py`, `a_i_b.py`,
`cattrib2.py`, `cattrib3.py`, `d1_coupling.py`, `e1_critical.py`, `f1b.py`.
**Symptom:** location errors of exactly 0.00 km.
**Cause:** `SRC=(0,0,z)` with symmetric grids of the type `np.arange(-21,21.1,3.5)`, where zero
**is** a node. Measured with a step of 3.5 km: on a node 0.00 km, at half a step 1.75 km, on the diagonal
2.47 km, that is, the quantisation itself.
**Detection:** later systematic audit, not during the experiments.
**Actual scope:** the results whose endpoint is a **rate or a margin** are not affected,
because they do not use the distance between estimate and truth. That covers Appendices H, J, M, O, R, S and T,
that is, all the results that support the document. Affected: the location-error columns
of L, N and U. L and N already use continuous optimisation and are sound; in U the zeros must
be read as "within one cell".
**Correction:** documented as Appendix V.5, with the registered rule that synthetic sources
are drawn from a continuous distribution and never placed on search nodes.

### Negative check: undersampling of tau in `b1_moment.py`
The E4 pattern was suspected because the script used a step of 0.2 s with a wavelet of period 0.5 s, and its
least-squares fit is of matched-filter type. **It was verified and it does not change the conclusion**:
with a step of 0.05 s there are still 0/491 false events by residual and the isotropic test still fails
(142/491 against 134/491). Only the explained-energy ratio goes down, from 1.42 to 1.28, because the
false hypotheses also benefit from a finer search. The script was corrected to 0.05 s and
the numbers of Appendix V.3 were updated with the fine version, which is the less favourable one.

This check is recorded even though it came out negative, because a known error pattern that is searched for and
does not appear is information, and omitting it would leave the doubt open.

---

### E8. Absolute threshold on scores of different scale
**Scripts:** `task0_run.py`, `task0_common.py`, protocols P1 and P2 of Task 0.
**Symptom:** the constrained model gave a false-event rate of 1.00 in nine of ten cells, **including the
easiest one**: twelve stations with clean data, where the ablation and the two published associators
gave 0.00.
**Cause:** the ablation is a fraction of hits, so its maximum is exactly 1.000 in every
realisation. The constrained model multiplies it by a continuous polarisation term and its
maximum varies between 0.51 and 0.68 depending on the realisation. **There is no absolute threshold that works for
both**: the one that gives an operating point to one denies it to the other.
**Detection:** the 1.00 in the trivial cell was physically implausible.
**Correction:** the threshold becomes a **fraction of the maximum of each realisation**. It is the standard in
detection and the argument does not depend on any result.

### E9. Averaging over the wrong factor
**Scripts:** `task0_run.py`, protocols P1 and P2.
**Symptom:** two seeds gave irreconcilable results. The ablation went from 0.13 to 0.70
false events in the same cell; the constrained model from 0.17 to 0.97 with four stations.
**Cause:** each cell fixed **one** network and drew thirty clutters. Since `make_network` is drawn
with the seed, each seed was a different geometry, and **geometry dominates over the method**.
The averaging was over clutter when the large variance was in the geometry.
**Detection:** the discrepancy between seeds, not an internal check.
**Correction:** the network is drawn **per realisation**. With that the two seeds replicate cell by
cell, which in turn is the evidence that this was the cause.
**Note:** Appendix S had already measured that the observability constant varies 86× between
geometries. The error was foreseeable from within the document itself and was not foreseen.

### E10. Registered threshold below the chance hit rate
**Scripts:** `task1_decide.py`, registered criterion of Task 1.
**Symptom:** none visible: the criterion looked reasonable and would have run without error.
**Cause:** the criterion required a coherent fraction of the clutter below 0.20 to declare "LOW".
A direction unrelated to the event falls within 20° by chance with probability 20/90 = 0.22
if the azimuths are uniform, and 0.24-0.33 with the azimuthally clustered events that were measured. The
threshold was **below the expected value of clutter with no coherence at all**: an observed
fraction could cross it through fluctuation alone, and passing or not passing measured the geometry.
**Detection:** when comparing the clutter fraction with its permutation null in the diagnostic.
**Correction:** the criterion is applied to $\eta=(f-f_0)/(f_{\rm ctrl}-f_0)$, with $f_0$ exact by
permutation within station. $\eta=0$ for unrelated clutter and 1 for clutter indistinguishable
from real arrivals, whatever the geometry. Validated on synthetic mixtures with known
proportion.

### E11. Clutter reference chosen by proximity in origin time
**Scripts:** `task1_analyze.py`, `label_association`.
**Symptom:** at ZAR only 50 of 753 false picks remained usable, and the coherence of the clutter
systematically approached chance.
**Cause:** the reference event of a false pick was the nearest one in origin time. With the
four-station catalog (4602 events in 90 days, one every 28 minutes) that event was usually
another one, distant, or even later than the pick. Coherence was measured against an event with no physical
relation to that stretch of trace, which biases the result towards "random clutter": precisely the
conclusion that the task was testing.
**Detection:** the disproportion between unassociated picks and usable picks.
**Correction:** reference = most recent event within 200 km whose predicted P precedes the pick.
The associations do not change (550/550 identical in the test); only the clutter reference changes.

### E12. Gate without a minimum sample size
**Scripts:** `task1_decide.py`.
**Symptom:** CHI "passed" the control with 0.636 and n=11; with n=78 it gives 0.256. A clutter of n=1 gave CI95
[1.000, 1.000] and would have produced "HIGH coherence", the strongest possible conclusion from the
weakest possible data.
**Cause:** a bootstrap with very few points only resamples the same values, and a fraction
over eleven picks has a standard error of 0.15.
**Detection:** the degenerate interval.
**Correction:** minimum of 30 in the control and in the clutter, and an "inconclusive" label distinct from
"fail" and from "intermediate".

### E13. Control that mixes phases
**Scripts:** `task1_decide.py`.
**Symptom:** the control mixed two populations with different behaviour: at CHI the P picks
give 0.256 and the S picks 0.677.
**Cause:** the control required the S picks to point at the event. The motion of an S wave is not
obliged to do so (SH is transverse to the propagation direction), so an S pick does not test the
estimator in any sense. Why the S picks of CHI align better than its P picks remains unexplained.
**Correction:** control with P picks only.

### E14. Best of many configurations on few data
**Scripts:** `diagnose_polarization_config.py`, sweep A5.
**Symptom:** at SJC, with 15 picks, a fixed rotation of 165° raised the control from 0.467 to 0.667 and the
script declared it "fixed rotation detected".
**Cause:** 72 transformations tested on 15 points: one of them hits by chance. With 19 picks the
same search gives 0.579 and the "rotation" disappears.
**Detection:** the warning was made before having more data, and the additional data confirmed it.
**Correction:** every selection among configurations is validated on two disjoint halves of events;
the sweep of fixed transformations is read only with sufficient n and never as an isolated finding.

### E15. Ties of a discrete score read as a physical effect
**Scripts:** `task0_common.py` (TOL ablation), Appendices M, O, X and T.2.
**Symptom:** in one-sided networks the ablation produced false events in 30-43% of the
realisations and the prototype with polarisation in none, with $p\le0.0039$ on two seeds.
**Cause:** the ablation scores a fraction of hits; several nodes reach exactly the maximum
and with a threshold of 1.00 of the maximum all of them are emitted. The prototype multiplies by a continuous term and its
maximum is unique. The difference was tie-breaking, not polarising.
**Detection:** audit of the catalog emission code (`catalog_from_score_field`), later
confirmed with three ablations that break ties without polarisation (Appendix Y.1): all give zero
false events on both seeds.
**Correction:** the reference ablation of G1 becomes a continuous kinematic score (residual
kernel), and every contrast is made against the best score without the channel, not against the simplest one.

### E16. Weak comparator because of an overly favourable generator
**Scripts:** the entire synthetic generator (`task0_common.generate`).
**Symptom:** a Gaussian kernel of time residuals, without polarisation, equals or outperforms the
prototype with polarisation in all cells (Appendix Y.1).
**Cause:** the synthetic picks have no time error, so the true source has
zero residual at all arrivals; any continuous residual score separates it almost
perfectly. The hard tolerance of 0.6 s of the ablation threw that information away.
**Correction:** pick-time error and velocity error as design factors (Appendix Y.2); the
confirmatory generator must include pick errors with the distribution measured on real data.

### E17. Threshold relative to the maximum, incompatible with empty windows
**Scripts:** `task0_common.catalog_from_score_field(relative=True)`, protocol P3.
**Symptom:** none on synthetic data, because each realisation had exactly one source.
**Cause:** normalising by the maximum of the window itself makes a window without an event always emit
at least one. The protocol could not run on real continuous data, which are mostly
empty windows, and with one source per window the catalog endpoint reduced to that of the
argmax.
**Correction:** windows with 0, 1 or 2 sources and an absolute threshold calibrated separately (Appendix Y.3).

### E18. Asymmetry in the calibration of the baselines
**Scripts:** `task0_run.py`.
**Symptom:** oscillations of PyOcto and GaMMA between seeds in the same cell.
**Cause:** the operating point of the baselines was chosen on the evaluation realisations
(favours the baselines) and only one knob was swept, with the rest of the hyperparameters at their defaults
(disfavours them).
**Correction:** declared grid of hyperparameters evaluated only in calibration, frozen before
the evaluation (`task1jk_controls.py`).

### E19. A generator that gives the rival channel information that no instrument has
**Scripts:** `task0_common.generate` and everything that uses it (Appendices G to Y.3).
**Symptom:** in Appendix Y.1 to Y.3 polarisation contributed nothing over the best kinematic
score, and this was written up as a result.
**Cause:** the generator detected all true arrivals with exact time, even with amplitude
over noise of 1.1 to 2.5, while polarisation was measured on that same noisy waveform
(AUC 0.78). The kinematic channel was simulated at infinite signal-to-noise ratio and the polarisation channel at
the real one.
**Detection:** the author's question of whether the negative result was an implementation artefact; the
channel and the generator were audited with the same suspicion previously applied to the positives.
**Correction:** detection and time error derived from the same signal-to-noise (Appendix Y.5). The
negative stops being conclusive; what remains in its place is measured in Y.6 and Y.7 (see E21).

### E20. Combining channels with an ad hoc product instead of a likelihood
**Scripts:** `task0_common.score_grid` (multiplicative C-IV), PHYS-K of `task1j2_noise.py`.
**Cause:** multiplying a count or a kernel by a cosine is not the likelihood ratio; it weights
the evidence badly and does not calibrate it with the data.
**Correction:** per-arrival likelihood (9d), with the polarisation ratio f/g estimated in
calibration and summed in logarithmic form (PHYS-LR).

### E21. A positive read at the nominal budget and against a weak reference
**Scripts:** `task1o_consistent.py`, `task1o_baselines.py`, `analyze_1o.py` (Appendix Y.5).
**Symptom:** PHYS-LR beat KIN-LR by 0.18 to 0.25 of recall in five of six cell-seeds and
GaMMA and PyOcto in the same ones, and this was written up as evidence for G1.
**Cause:** five defects, found by an independent verification. (1) The gain was read at the
nominal false-event budget, which PHYS-LR exceeded in evaluation; (2) there was no negative control,
so it was not known whether the improvement came from polarisation or from having one more knob; (3) KIN-LR,
the reference, was weaker than GaMMA; (4) GaMMA and PyOcto searched in a wider domain than the
in-house methods and with a minimum of picks from 5; (5) the change of generator and the change of form of the
likelihood were introduced together.
**Correction:** Task 1r (Appendix Y.6): recall at realised false events, control with permuted cosines,
exact condition, same domain, CI by bootstrap. With 60 windows the result remained indeterminate;
with 200 windows and the exact curve the advantage over the kinematic reference and over GaMMA reappears,
smaller than the one read in 1o, and the second verification showed that the shuffled control was not null
(Task 1w, Appendix Y.10).
**Lesson:** the same scrutiny that turned the negative of Y.2 and Y.3 into E19 should have been applied to the
positive that replaced it, before writing it up. It was written up first.

### E22. A non-deterministic baseline
**Scripts:** `task1o_baselines.py`, `task1r_baselines.py` (version without seed).
**Symptom:** the recall of GaMMA on the same windows and with the same configuration differed
between two scripts (0.36 to 0.43 against 0.40).
**Cause:** GaMMA initialises the mixture components with the global numpy random generator
without fixing a seed when it runs in a single process (`gamma/utils.py`, line 638); its catalog depends
on everything the process executed before. In 80 calls, 15 changed their hits or false events when the
state of the generator changed (`check_gamma_det.py`).
**Correction:** fixed seed before each call, the convention that GaMMA uses in its
multiprocess mode. With it the two code paths give the same catalog in the 80 calls.

### E23. A FROC curve sampled with jumps at its strict end
**Scripts:** `task1r_fair.sweep` and everything that uses it (Tasks 1r, 1t, 1u).
**Symptom:** in 1t, WAVE-E at 4 stations had points of its curve at 0.00 and 0.14 false events per window
and none in between; its recall at cap 0.1 was read as 0.01 although at 0.14 false events it was already 0.12.
**Cause:** the thresholds were 80 quantiles (0.85 to 1.0) of all the pooled grid scores;
almost all of them fall below the peaks that decide the strict end of the curve. The reading at cap
0.1 depended on where the thresholds fell and could create or erase a difference between methods.
**Correction:** `froc_exact.py`. The catalog at a threshold is the prefix of the complete greedy list of
each window, so evaluating each score of those lists as a threshold gives the exact curve; it was verified
equal to the original emission and count at each threshold. The rules are evaluated on the exact curve
and the reading with the registered sweep is also reported. The amendment was written after seeing
that T-b failed with the registered sweep; with it R-d, U-a and U-b also fail, and with the exact curve
they hold. Both readings and the order of events are in Appendix Y.6.

### E24. A generator that radiates S in the direction the method predicts
**Scripts:** `task0_common.generate` (all tasks up to 1u).
**Symptom:** WAVE-P, which projects the S wave onto the SV direction, outperformed the same migration without
direction by 0.13 to 0.21 of recall, and PHYS-LR its control by 0.08 to 0.17.
**Cause:** the generator gave every S wave pure SV polarisation, without a source mechanism. A real
earthquake radiates S with an SH part that depends on its mechanism, unknown to the method. The
evaluated channel received the exact information that its model assumes.
**Detection:** second independent verification, with an exploratory run with double couples.
**Correction:** Task 1w (Appendix Y.10): random double-couple sources. The per-pick advantage no longer
holds against a control that breaks the geometry; the mechanism-independent pick-free form
(P along the ray, S in the transverse plane) retains it.
**Lesson:** it is the same as E19 in the opposite direction: a generator can give the own channel,
and not only the rival one, information that real data do not have.

### E25. A negative control that retained part of the signal
**Scripts:** `task1r_fair.py` (PHYS-SH).
**Symptom:** the control with shuffled cosines outperformed KIN-LR, and this was interpreted as flexibility
of the likelihood.
**Cause:** the permutation reordered the polarisation ratios among the station-phase slots of
the same hypothesis; the set of ratios of the hypothesis survived and only their pairing with
the time residuals was lost. Moreover it had exactly the same variants as KIN-LR, so there was no
additional flexibility to explain the gain.
**Correction:** control with the directions predicted by another grid node more than 9 km away
(PHYS-NODE, WAVE-PT-NODE).

### E26. Clutter picks that no picker would detect
**Scripts:** `task1o_consistent.degrade`, `task1w_mech.degrade_dc`, `task1x_mtglrt.degrade_x` (Tasks 1o
to 1x).
**Symptom:** GaMMA recovered very few events even with three or more true phases detected, and the
pick-free methods outperformed it by margins of 0.3 to 0.7 of recall.
**Cause:** the consistent generator of 1o applied the detection rule to the true arrivals but
not to the clutter. The clutter pulses have a median amplitude over noise of 0.52 and a picker would
detect them with probability 0.003, but all of them were delivered, with exact time: about 90 per
cent of the picks that GaMMA received were undetectable. The waveform methods saw that same
clutter as pulses below the noise. It is E19 in the opposite direction: the generator took away from the
rival what it gave to the own channel.
**Detection:** third independent verification (Task 1x), which repeated GaMMA with the clutter thinned
by the same rule: its recall in the one-sided cell went from 0.18 to 0.68.
**Correction:** Task 1y: the clutter follows the same detection and time-error rule as the
true arrivals, in two consistent regimes (weak clutter that is almost never detected and strong
clutter, visible and detected). The margins over GaMMA of Appendices Y.5 to Y.10 are inflated by
this defect; the comparisons between waveform methods and their controls do not depend on the
picks and do not change.

### E27. A synthetic pick model more optimistic than real pickers
**Scripts:** `task1y_clutter.degrade_y` and the same law in `task1o_consistent.degrade`,
`task1w_mech.degrade_dc` and `task1x_mtglrt.degrade_x`: the picks that GaMMA and PyOcto received in
Tasks 1o to 1y.
**Symptom:** none visible. The detection law $\Phi((a-2)/0.5)$ with a time error of 3 times the
Cramér-Rao bound was adopted in 1o and had never been compared with a picker running on the
waveform.
**Cause:** measured with `gen_realistic.py` and `real_picker.py` on 160 windows of seed 999 that did not
enter the calibration of the pickers (1608 true arrivals), with `degrade_y` applied to the
same arrivals and amplitudes. degrade_y picks 42.0 % of the true arrivals, the calibrated STA/LTA
14.3 % and PhaseNet 4.1 %; the median time error is 0.033 s against 0.078 and 0.080 s. degrade_y
detects half of the arrivals with amplitude over noise of 2; the STA/LTA, around 5; PhaseNet, around
8 to 10. The false picks per window, 14.27, 6.93 and 2.08, mix two regimes. With weak clutter
degrade_y delivers 0.05 to 0.07 false picks per window and the real pickers 1.9 to 3.9 (noise, coda and swapped
phase): there it is optimistic for GaMMA in everything. With strong clutter it delivers about 28 clutter picks per
window, against 7 to 8 from the STA/LTA and fewer than one from PhaseNet, because an isolated Ricker pulse does not
look like an earthquake: there it is optimistic in recall and in time but delivers more clutter, so it is not
uniformly optimistic. The collapse of PhaseNet is above all domain mismatch: without noise it gives
probability 0.8 to 0.9 to almost all arrivals, but the white noise of the generator occupies 0 to 8 Hz,
above the signal band, the reverse of the seismic noise it was trained with. Its figure is a
pessimistic bound, not its performance on real data.
**Detection:** when replacing the parametric model by real pickers on a mismatched generator
(Task 1p in part, preparation of Task 2A).
**Scope:** with weak clutter, GaMMA received in 1y more and more precise picks than a
real picker would give it, so this defect does not inflate the margins of MT-FW over GaMMA in that regime; with strong
clutter the direction of the bias is undetermined.
**Correction:** in Task 2A the associators receive the STA/LTA picks with the configuration
frozen on seed 999, and degrade_y remains as a secondary optimistic bound (GaMMA-DY, rule 2A-c).

### E28. A false-alarm rate guarantee written before measuring it
**Text and scripts:** Section 7.5 of the proposal (point 4) and Section 3 of `method_design.md`;
measured with `rft_threshold.py`.
**Symptom:** none. The claim, that the expected Euler characteristic of a chi-squared field with
6 degrees fixes the per-window false-alarm rate before seeing data, was analytical and the predicted
rate had never been compared with the realised one.
**Cause:** three, independent. (1) The implemented statistic is not chi-squared with 6 degrees, not
even at a single hypothesis: taking the energy peak of a window of 9 samples of filtered noise raises
the mean of $2T$ to 12.6-13.5 (6 for $\chi^2_6$), with a Kolmogorov-Smirnov distance of 0.52 to 0.58.
(2) The field decorrelates within one step of the spatial grid (correlation at one step of 0.09 to 0.21),
so the continuous EEC counts 330,000 to 640,000 effective tests where the data give 8,700 to 19,000;
with curvatures estimated by roughness it overestimates $P(M>u)$ by 20 to 60 times and gives a realised rate of
0.002 when 0.25 is requested, and the variants that respect the grid remain conservative by about 4 times.
(3) Strong clutter shifts the maximum by +8 to +9 and every threshold calibrated on Gaussian noise fails,
including those from random fields: the conformal one gives rates of 0.60, 0.78, 0.90 and 0.96 where 0.01,
0.05, 0.10 and 0.25 were requested, and the most conservative EEC 0.14, 0.27, 0.35 and 0.47.
**Detection:** Task 1z, component A: 7,200 noise-only and clutter windows in three cells, with a fixed
network and with a network per window, disjoint calibration and test.
**Correction:** the primary threshold is empirical: split conformal (finite-sample guarantee under
exchangeability of windows, Proposition AA.5(c)) or GEV fitted to the maxima of about 300 noise-only
windows from the same network and the same regime. On Gaussian noise both remain within the predictive
band in 24 of 24 combinations; recalibrated on strong-clutter windows of the same law
they again remain within it in 12 of 12. It is recalibrated every time the regime changes. The EEC and Bonferroni
remain as a prior conservative bound on Gaussian noise, not as a guarantee.

### E29. Theoretical statements stronger than the code
**Text and scripts:** Sections 7.5 and Y.11 of the proposal versus `task1x_mtglrt.fw_record`; found while
writing Appendix AA of the proposal (self-verification in `results/components/appendix_AA_selfcheck.txt`).
**Symptom:** none in the figures. The text described an idealised algorithm and attributed to it guarantees that
the implemented one does not have.
**Cause:** five discrepancies. (1) The text said that deflation removes the side lobes without a
suppression radius; `fw_record` excludes, from the second atom onwards, the nodes within 5 km epicentral distance of
the accepted atoms (`min_sep = 5.0`). (2) The $O(1/k)$ rate of Frank-Wolfe with an exact oracle was invoked;
MT-FW is a group matching pursuit with a physical dictionary:
least-squares fit of the chosen atom without refitting the previous ones, re-selection of the peak and cap,
without penalisation. It has no rate, and with the cap active the residual may fail to decrease (Proposition
AA.4(e)-(h)). (3) It was said that $\lambda^2\eta^2$ is a chi-squared field with 6 degrees: it holds at each
hypothesis for the idealised certificate at $\mu=0$; the Euler characteristic formula additionally requires
hypothesis (G), which the moment-tensor subspaces do not satisfy in general, and the implemented
statistic is not chi-squared even at a single hypothesis (E28). (4) It was said that $\eta=1$ at each declared
event: it holds for the minimisers of the penalised problem (Proposition AA.1), not for the output of
the deflation. (5) The curvature of the non-degenerate source condition was identified with the local Gram
matrix that the observability certificate computes; the curvature is the efficient Fisher
information with the moment tensor as nuisance, bounded below by the leading term of
Section 6.1 with weights equal to the energies radiated by the real mechanism (a station on a nodal
plane contributes nothing) and above by the information with known mechanism (Proposition AA.6(b)-(d)); the
global part of the condition is not a local statement.
**Detection:** while writing the proofs of Appendix AA and rereading each step against the code.
**Scope:** no measured number changes; it affects the guarantees and the description of the method, not
the recalls of Y.11 and Y.12.
**Correction:** the proposal states the 5 km radius, describes MT-FW as a group matching
pursuit and restricts the rate, the chi-squared field, $\eta=1$ and the non-degenerate source condition to the
hypotheses under which Appendix AA proves them.

### E30. A novelty claimed without the closest prior work
**Text:** Section 7 of the proposal and Sections 2 and 7 of `method_design.md`.
**Symptom:** none in the figures. The proposal presented as new the joint, sparse and moment-tensor-aware
detection and location, and cited as prior work only migration with
moment tensor in dense arrays and GRiD MT at long periods.
**Cause:** the search for prior work covered moment-tensor migration and the long-period
scan, but not the sparse inversion of simultaneous sources. Vera Rodriguez, Sacchi and Gu (2012, *GJI*
188(3)) already recover jointly
and sparsely the origin time, hypocentre and moment tensor of simultaneous sources with a
block orthogonal matching pursuit over a dictionary of Green's functions, and GRiD MT (Tsuruoka et al.,
2009) scans location, time and moment tensor operationally since 2003. Also already
published were mechanism-corrected stacks (Anikiev et al., 2014; Chambers et al., 2014; Zhebel and
Eisner, 2015), master-event corrections (Grigoli et al., 2016) and extreme-value
thresholds (Hirano et al., 2022).
**Detection:** verified novelty audit.
**Scope:** no measured number changes. What the thesis claims and the list of pending
comparators change: Task 2A does not compare with the block sparse inversion or with the mechanism-corrected
stack.
**Correction:** repositioning. What remains as a contribution is the off-grid formulation as a group
BLASSO with vector moment-tensor amplitudes, the identity between its certificate and the generalised
likelihood ratio of the adaptive subspace detector (a modest unifying result), the
certified search, the conformal thresholds and the benchmark against GaMMA, PyOcto and QuakeMigrate with
overlapping events. The block sparse inversion and the mechanism-corrected stack move to the
list of pending comparators, and Task 2B measures the contribution of the radiation operator inside the
same solver.

### E31. A different BLAS kernel between resumed sessions
**Scripts:** `task2a_confirm.py` (resumed sessions of seed 7331) versus `task2b_ablation.py`;
outputs in `results/task2b/log2b_attempt1.txt`, `task2b_control_failed_attempt1.txt` and
`task2b_kernel_check.txt`.
**Symptom:** the first attempt of Task 2B stopped at its exact control after passing 32 networks.
In the 6-station networks of seed 7331 written in the resumed sessions of 2A, the recomputed records
of MT-FW+SC differed from the saved ones in the last 1 or 2 digits (relative ≤ 6e-15, same
hits and same order), and in one network a station correction differed by 8.9e-16 s.
**Cause:** numpy's OpenBLAS (compiled with DYNAMIC_ARCH) chooses the compute kernel when the
process starts. On this machine it chooses Prescott by default, but the 18 networks that 2A wrote in its resumed
sessions of 7331 (4 stations in both regimes; 6 stations, the 5 networks with weak clutter and
networks 2 to 4 with strong clutter) were computed with the AVX-512 SkylakeX kernel. The two kernels round
matrix products differently. With `OPENBLAS_CORETYPE=SkylakeX` those 18 networks are reproduced bit for
bit; with Prescott, Haswell, Sandybridge or Zen they are not; the other 42 are reproduced bit for bit with the default
kernel.
**Detection:** the preregistered exact control of 2B (regenerated MT-FW+SC records identical to those
saved in `task2a_parts`), which requires bit-for-bit equality and not a tolerance.
**Scope:** no reported number changes. With the 18 networks recomputed with the default kernel
(exploratory), MT-FW+SC differs by at most 6e-15 relative without changing hits or order; MT-FW-FREE+SC
diverges in later Frank-Wolfe steps in some windows (up to 17 % in score, because of nearly
tied candidates) and changes hits or order in 2 of 360 windows; the recall at 0.1 and 0.25 is identical for
both methods in the 4 affected cells. The per-network results of 2A mix two kernels, with an effect
on MT-FW only at the level of rounding.
**Correction:** the final run set `OPENBLAS_CORETYPE` per network (SkylakeX in those 18, the default one in
the rest) and recorded the kernel in each part file it wrote; all networks passed the same
exact control, without relaxing any threshold or tolerance. The 32 parts of the first attempt do not carry the
field, but they ran with the default kernel.

### E32. Inherited code that could not be run as written
**Scripts:** `task2c_methods.py` and `task2c_run.py`, inherited from a previous session; review in
`results/task2c/task2c_report.txt`, section 3.
**Symptom:** none in the figures, because no evaluation run was made with that code. On
reviewing it, three defects appeared that prevented running it as it was or would have made it run incorrectly.
**Cause:** (1) `task2c_run.py validate` imported `task2c_validate`, which did not exist: the methods had no
validation. (2) The OpenBLAS kernel check only compared the environment variable with "SkylakeX".
The default numpy kernel on this machine had changed from Prescott, the one of 2A and 2B (E31), to SkylakeX,
so with the variable unset the 42 networks of 2A written with Prescott would have run with a kernel
different from that of 2A, which the exact control does not tolerate. (3) The restricted solver of GL-FW required a relative
KKT violation of 1e-11, with gradient stopping at 1e-12 λ, which floating point does not reach: the Armijo
test failed at the level of rounding and the solver looped up to 200 outer iterations, 39 s for
6 values of λ of one window and minutes for longer paths. In addition, an active block with a
certificate larger than λ raised an error instead of re-solving.
**Detection:** review of the inherited code, with the derivations redone, before any run of
seeds 2024 and 7331. The derivations (dictionary, Gram with the autocorrelation, selection, refit
and score of BOMP, MTSTACK, certificate of the group lasso and formula of the gap) were correct.
**Scope:** no reported number changes; the defective code did not produce any result.
**Correction:** `task2c_validate.py` (V1 to V4) was written. `OPENBLAS_CORETYPE` is set per process
(SkylakeX in the 18 networks of E31, Prescott in the others, in all of part (ii) and on seed 999); each
network checks with threadpoolctl the kernel actually loaded and saves it in its file, and each process only
launches the networks of its kernel. The restricted tolerance became 1e-9, a thousand times stricter than the
Frank-Wolfe admission tolerance (1e-6); the Armijo test received a rounding slack of 1e-14 |F|; the
dense O(n³) whitening was replaced by a block product O(n²), and the error, by a bounded counter
of re-solves that was never triggered (0 in 1200 windows). Afterwards all exact controls passed.

### E33. A method presented as the one that performs before measuring it
**Text:** Section 7.5 and Conclusion of the proposal in English (G1 is "the coherent moment-tensor detector
that certifies the group BLASSO of Section 7.5"); in their versions prior to Task 2C, §2, §3 and §8
of `executive_summary.md` and Sections 2, 3 and 7 of `method_design.md`.
**Symptom:** none in the figures until Task 2C. The documents presented the group BLASSO as the
method of the thesis and MT-FW as its greedy approximation, which implied that solving the convex
problem would perform at least as well as its approximation. The independent audit pointed out that the defended method
was not implemented and that all the evidence came from another algorithm.
**Cause:** the claim rested on theory (the certificate of the convex problem is the statistic that
MT-FW maximises) and not on a measurement. The convex problem on the grid reads each station at the exact predicted
time: it does not contain the ±0.25 s peak search, the cap or the greedy deflation of MT-FW, which
are what absorb the discretisation of the 3 km and 0.1 s grid and the residual time errors.
**Detection:** the descriptive rule 2C-b, preregistered in response to the audit. The group BLASSO on
the grid, solved with a duality gap ≤ 1.4e-8, falls below MT-FW in 23 of 24 cases, by
−0.064 to −0.605, and at the level of the block pursuit or somewhat below. The exploratory diagnostic, with
the sources on the nodes of the ideal generator, attributes the difference to the grid.
**Scope:** no measured number changes. What the thesis claims changes: the convex problem is the framework
and the certificate, the detector is MT-FW, and the empirical claim is robustness to grid and
time errors on a coarse grid. 2C-b tests only the restriction to the grid: the off-grid version
is not refuted, it is unimplemented.
**Correction:** in the companion documents in Spanish, the convex problem becomes framework and certificate,
MT-FW becomes the detector, and the sliding Frank-Wolfe, or a time-tolerant formulation, becomes an open
question of the mathematical line (H5 of `method_design.md`). In the proposal in English the same
change is in Appendix Y.20a and in Sections 1.6, 4.3, 7.5, 20.0 and 23, and the title and G1 describe the
family of time-tolerant detectors and not the certificate of the convex problem (E35).

### E34. A preregistration time written without consulting a clock
**Text:** entry "Task 2B" of `PREREGISTRO_1j_1k.txt`, with the stamp "14:30 UTC", and its copy in
`tasks_and_real_data_plan.md`; Appendix Y.19 and Section 4.3 of the proposal in English cite that time only
together with the erratum.
**Symptom:** the time of the entry is later than the creation of the first file of the task:
`task2b_ablation.py` was born at 14:16:49 UTC and its first output, `task2b_selfcheck.txt`, at 14:17:01
UTC, so the registration appeared to have been written after the implementation started.
**Cause:** the stamp was typed from memory, without reading a clock, and no hash fixed the content of the
preregistration at that moment.
**Detection:** the second independent audit, comparing the stamp with the birth dates of the
files of 2B.
**Scope:** the file system cannot date the text of 2B, because it only records the last
modification of the preregistration. The local transcript of the working session records the addition at
14:14:04 UTC, before `task2b_ablation.py` (14:16:49 UTC), but that transcript is also local and is not
an external stamp, and no hash fixes the text at that moment. According to that transcript, the typed time,
14:30, was written sixteen minutes before `task2b_report.txt` was born (14:30:37 UTC), so its
coincidence with it, which the third independent audit pointed out, is fortuitous. The preregistered status of 2B, and that of 2C,
is therefore classified as "registered in a local append-only file, with a date not externally
verifiable" (E36). No number changes.
**Correction:** dated erratum in the preregistration (26 September, 21:13 UTC) and in the documents.
From Task 2D onwards each entry carries the time read from the clock (`date -u`) and a SHA-256 of the
preregistration file in `PREREGISTRO_hashes.txt` (21:13:04, 02:41:50 and 03:54:55 UTC for the entry of 2D, its budget
note and its result), and the hashes of the analysis code are recorded before running
(`task2d_code_hashes.txt`).

### E35. An advantage attributed to the method that belonged to a transferable heuristic
**Text:** registered interpretation of rule 2C-a and its reading in Appendix Y.20 and Y.20a of the
proposal, in `tasks_and_real_data_plan.md` (Task 2C) and in `executive_summary.md` §3 and §5b: since 2C-a
held, "the advantage of MT-FW is not attributed to the moment tensor with deflation, which is already in Vera
Rodriguez et al. (2012), but to the elements specific to MT-FW".
**Symptom:** in 2C, MT-FW+SC outperformed BOMP+SC by +0.08 to +0.52 at 0.25 in all cells, and this was read
as an advantage of the MT-FW method over the block inversion.
**Cause:** the comparison of 2C had unequal temporal tolerance: MT-FW reads each station and phase at the energy
peak within ±0.25 s, and the comparators at the exact predicted sample. The peak search was counted among
the "elements specific to MT-FW" when it is a heuristic that any method of the family can use.
The exploratory diagnostic of 2C itself pointed to the peak search, but the conclusion about the method
was written before giving the comparators the same tolerance.
**Detection:** rule 2D-a, preregistered after the second independent audit. With the same peak
search, BOMP-PK+SC matches MT-FW+SC: differences at 0.25 from −0.155 to +0.091, all hierarchical CIs contain
0, and the rule fails in both regimes. Without its cap, MT-FW has almost the same recall as BOMP-PK in all
cells, and the cap has no significant effect (−0.091 to +0.131).
**Scope:** no measured number changes. The narrow part of the reading of 2C survives: what separated
MT-FW from the exact-time comparators was the peak search. What does not survive is that this makes
MT-FW distinctive: MT-FW is not a better detector than the block inversion with the same tolerance, and the
empirical advantage over the pick-based associators and QuakeMigrate belongs to the family of
time-error-tolerant moment-tensor detectors with deflation (BOMP-PK+SC outperforms QuakeMigrate by +0.17
to +0.49).
**Correction:** the proposal (Sections 1.3, 1.6, 3.6, 4.3, 4.5, 7.5, 13.0c, 17, 19, 20.0 and 23, and
Appendix Y.21) and the companion documents state that the advantage is attributed only to the peak search,
that BOMP-PK is the comparator to beat and that the distinctive methodological contribution becomes the principled
treatment of temporal tolerance: the group BLASSO with per-station and per-phase shifts as
nuisance parameters, off the grid, which only counts if it outperforms the peak search and BOMP-PK.
General rule: an advantage is not attributed to a method as long as its comparators do not have the same
transferable heuristics (tolerances, corrections, search windows).

### E36. A preregistration without a verifiable stamp and an undeclared edit of the analysis
**Text:** Section 4.3 of the proposal in English ("For Tasks 2A and 2C the registered times precede every
file of the task in the file-system record" and "the preregistered status of Task 2B rests on file-system
times only"), Appendices Y.18 to Y.20 of the supplementary volume, E34 of this log, and
`results/task2b/task2b_report.txt`, which did not mention the edit of `analyze_2b.py`.
**Symptom:** the third independent audit checked with `stat` that the preregistration was last modified
on 27 September at 03:54:55 UTC, so the file system cannot date
any text within it; that for 2C the claim about the file-system timestamps is not
verifiable; and that `analyze_2b.py`, created at 14:19:07 UTC, was modified at 14:26:33 UTC, 18 s
after its first output (`analyze_2b_out.txt`, 14:26:15 UTC), with the output rewritten at 14:26:43
UTC and without any declaration in the task report.
**Cause:** the entries of 2A to 2C were added to a local file with the typed time and without a hash, and the
proposal read as date evidence some file-system timestamps that only date whole
files. The edit of `analyze_2b.py` added an exploratory analysis (the fraction of bootstrap replicates
in which rule 2B holds) after seeing the first output, and it was not noted.
**Detection:** the third independent audit (`audits/independent_audit_3.md`, D6 and
inconsistencies 8, 9 and 16).
**Scope:** no rule and no number changes. The local transcript of the session shows the complete
edit: three replacements that create a container of per-replicate indicators, fill it with the rule
applied in each already computed replicate and write one line "Exploratory: fraction of replicates in which
rule 2B holds" per regime (0.59 and 0.62 with weak clutter; 0.98 and 0.97 with strong clutter). The
previous script was reconstructed by reverting those replacements and it was run on the same records: its
text output is identical line by line to `analyze_2b_out.txt` except for those two added lines, and its
JSON is identical to `analyze_2b_out.json`; the edit does not consume random numbers. Tasks 2B and 2C
are now classified as "registered in a local append-only file, with a date not externally
verifiable"; Task 2A, the same, with a consistent file-system record (first file of the
task at 23:33:22 UTC on 25 September; operating points of 2A-e frozen at 01:12:25 UTC,
before the first evaluation output, at 01:12:53 UTC); and Task 2D is the first one with SHA-256
records. None of the four has an external stamp.
**Correction:** Section 4.3 of the proposal and Appendices Y.18, Y.19 and Y.20 of the supplementary volume, with the
above classification and the edit of `analyze_2b.py` declared in Y.19; E34 of this log. General
rule: every edit of an analysis script after its first output is declared with time, reason and
output difference, even if it is exploratory; and from the next task onwards each registration carries an external
stamp (signed commit in a public repository, registration on OSF or an RFC 3161 timestamp), because a hash
saved on the same disk proves integrity of the content, not the date. Tasks 2E and 2F carry a
record of their content outside the container (the hash delivered to the user through the chat before the code or
the run); it is neither public nor third-party, and the public stamp remains pending for the test with real
data.

### E37. A mechanism read in the pilot more broadly than the confirmed effect supports
**Text:** `results/explorations/REPORT.txt`, section 2 ("per-station clipping is what makes it resolve dense
overlapping sequences"; "a few near-station or site-amplified arrivals dominate the unclipped statistic and
pull atoms to wrong nodes"), and the interpretation of `results/task2f/task2f_report.txt`, section 7, which
took that reading without bounding it: the gain of 2F was read as a better detection of the
overlapping events within each window.
**Symptom:** the effect of the cap was described as resolution of dense sequences, when the ablation had only
shown which element carries the gain (the cap), not why, and the level of the cap had never been
varied below the value fixed in 1x.
**Cause:** the mechanism was inferred from the ablation and from a plausible conjecture about amplified arrivals,
without an experiment that removed the proposed cause. An ablation identifies the element, not the mechanism.
**Detection:** the independent verification of 2F (`results/task2f/verif_2f/task2f_verification.txt`,
section 5, exploratory, seed 999). The gain peaks at the fixed cap (CAPV = 6: +0.070 at
0.25) and reverses with caps of 3 or less (−0.056) and disappears with 24 (−0.008); without clutter the gain
rises (+0.137), so it is not the clipping of clutter peaks; without clutter and without lognormal site factors
it falls to +0.023 [−0.040, +0.069]; within the window the distinct hits among the first atoms rise
only +0.001 to +0.012, while among the top 20 % of first-atom scores MT-FW hits 0.90-0.96
and BOMP-PK 0.79-0.84.
**Scope:** no registered rule or number changes; 2F-a, 2F-b and 2F-c still hold. The reading
changes: the gain is robustness to unmodelled per-station site amplification and better calibration
of the score across windows under a global threshold, tied to the level of the cap, and not a better separation of the
events within the window.
**Correction:** the proposal (Abstract, Sections 1.3, 1.6, 4.3 with hypothesis G1-R, 7.5, 17, 20.0 and 23) and
Appendix Y.24 and Y.24a state the mechanism with the scope of the verification, and the principled
generalisation of the cap (per-station site terms inside the certificate) becomes the next test.
General rule: the mechanism of a confirmed effect is stated only with the scope of an experiment that
removes the proposed cause, and the level of every fixed parameter is swept before attributing an effect to it. The
reading "per-station site amplification" of this scope and the site-term plan of this
correction turned out in turn to be an erroneous attribution, recorded in E39.

### E38. A hit rule without origin time in dense windows
**Script:** the hit rule of Tasks 2A to 2F (`task2c_methods.hit_of`: epicentre within ≤ 5 km of a
true event, with no condition on the origin time), inherited in the registration of 2F and applied by
`analyze_2f.py`.
**Symptom:** in the dense windows of 2F all events are less than 9 km from one another, so an atom
can "hit" an event with an origin error of 2.5 s (first saved window).
**Cause:** the rule was designed for windows with 0 to 2 events, where matching by location suffices, and it was
inherited without review into a design with 3 to 5 spatially clustered events.
**Detection:** the independent verification of 2F, section 2e. About half of the location-only
hits have an origin error greater than 0.5 s. With 5 km and 0.5 s, MT-FW+SC − BOMP-PK+SC at 0.25 is +0.045,
+0.024, +0.043 and +0.027 (weak 4242, strong 4242, weak 8888, strong 8888) instead of +0.072, +0.051,
+0.074 and +0.043, and the absolute recall drops to 0.09-0.13 (MT-FW) and 0.07-0.08 (BOMP-PK).
**Scope:** it is not a false positive: the rule was the registered one and it is applied equally to both methods, and the
direction of the gain and its relative size survive. The absolute gain of +0.04 to +0.07 and the absolute
recall in dense windows are inflated by the rule.
**Correction:** Appendix Y.24 (Table Y.56) and Sections 1.6, 13.0c, 14 and Appendix D of the proposal: the
gain is reported with both rules, and every future test with overlapping events also registers a rule
with origin time. General rule: the matching rule is reviewed every time the event density of
the design changes.

### E39. A per-window amplitude perturbation read as a station term, and an unsupported plan
**Text:** the mechanism reading of the verification of 2F as it stood in the proposal (Abstract,
Sections 1.3, 1.6, 4.3, 7.5, 17, 19, 20.0 and 23), in Appendices Y.24 and Y.24a and in the correction of E37:
"robustness to unmodelled per-station site amplification"; and the plan derived from it, a
certificate with per-station amplitude or site terms estimated from the catalog as a "physical generalisation
of the cap", which appeared as the next registered test, first mathematical deliverable (T2) and objective 1.
**Symptom:** the amplitude factor of the generator was described as a station term, and it was
planned to estimate it with the catalog, when in the generator there is nothing persistent to estimate.
**Cause:** `gen_realistic` draws the lognormal factors (sigma 0.5, one per station and phase) in each
call, that is, in each window and in each catalog event. The name "site factors" was taken at
its physical meaning without reading how they are drawn, and the experiment of the verification of 2F, which removed at the same time
the clutter and those factors, showed that the gain of the cap depends on them, not that the cap corrects a
station term. From a dependence it moved to a mechanism and from there to a plan.
**Detection:** the exploration C8 (`results/explorations/REPORT_C8.txt`, seed 999, exploratory), when writing
the generator with the factors returned: the deviation of the log-factor across the windows of one network is
0.42. With an added world of factors persistent per network, the catalog terms fall below the
cap (−0.040 [−0.078, +0.007]; −0.066 [−0.108, −0.028] in the unchanged generator), the true factors
without cap barely match it (−0.011 and −0.008) and with the cap in place they add +0.045 and +0.032; with dispersion 0.25
the oracle gains nothing over BOMP-PK (−0.001) and the cap still gains +0.043 [−0.026, +0.099]; no
amplitude treatment captures more than a quarter of the margin up to the realistic ceiling.
**Scope:** no registered rule or number changes; 2F-a, 2F-b and 2F-c still hold, and the
numbers of the verification of 2F are the same. The reading and the plan change: the mechanism is bounded
influence of each station-and-phase block against the amplitude error of the model, not correction of a site
term, and the certificate with site terms has no support. The error did not produce any false
result; it would have devoted the first year of the mathematical line to a direction that the generator itself
rules out. The evidence of C8 is exploratory and from a single seed.
**Correction:** the proposal (Abstract, Sections 1.3, 1.6, 3.6, 4.3 to 4.5, 7.5, 13.0c, 17, 19, 20.0, 23 and
Appendix D), Appendices Y.24, Y.24a, Y.25 and AA.8 and the documents in Spanish state the mechanism as
per-block bounded influence and replace the plan by the supported direction: a robust bounded-influence
certificate whose level is set from the noise distribution, preregistered as Task 2I, and
the open problem of the margin, of which about three quarters remain unexplained. General rule: before
giving physical meaning to a generator parameter, read how it is drawn and what persists across windows,
and a research plan does not rest on a mechanism that only has a measured dependence.

### E40. "Certificate" used for a statistic that is not the certificate
**Text:** the title of the proposal (*A Robust Subspace Certificate and a Preregistered Benchmark*), the
name of G1-R ("robust-certificate hypothesis"), the Abstract, Sections 1.3, 1.6, 3.6, 4.3 to 4.5, 7.5,
13.0c, 17, 19, 20.0 and 23, Appendix D and Appendices Y.24 to Y.26.
**Symptom:** the clipped statistic and the Huber statistic were called "certificate" while the same text
stated that the statistic with peak search and cap "is not that certificate" and Appendix AA.8
treated it as a deviation without theory.
**Cause:** the term was inherited from the convex problem, where the detection statistic is indeed the dual
certificate, and it was carried over to the object that performs without there being an optimisation problem of which that object is the
certificate, nor a false-alarm guarantee, nor a recovery guarantee.
**Detection:** fourth independent audit (weakness 4 and inconsistency 2).
**Scope:** no number changes. The title and the hypothesis promised a guarantee that has not been proven.
**Correction:** the deployed object is called "bounded-influence moment-tensor subspace statistic" in the
proposal and in the supplementary volume; "certificate" is reserved for the dual certificate of the convex problem
and for the observability certificate; the title becomes *Timing-Tolerant, Bounded-Influence
Moment-Tensor Detection of Overlapping Earthquakes: A Preregistered Benchmark and a Real-Data Test Plan*,
and G1-R is named "bounded-influence hypothesis". The one-sentence thesis of Section 1.3 describes the
greedy search with temporal tolerance and bounded influence and leaves the constrained optimisation as the
framework. Formulating the problem of which the statistic would be the certificate is work for T2.

### E41. A dependence declared removed while the arms that would measure it were unreported
**Text:** Sections 1.6, 4.3, 17, 20.0 and 23 of the proposal and Appendix Y.26: Task 2I "removes the
dependence on a tuned constant".
**Symptom:** the sentence claimed more than the registered primary arm supports. The level and
dispersion arms of 2I are not reported, and the previous exploration showed zero gain with a Huber level
of 8 and worse performance than the cap with dispersion 0.8.
**Cause:** "the level was not chosen by any detection result" was read as "the result does not depend on the
level".
**Detection:** fourth independent audit (weakness 9 and inconsistency 8).
**Scope:** rules 2I-a to 2I-c and their numbers do not change.
**Correction:** the text says that Task 2I replaces the calibrated cap by a level set by rule from
the noise distribution in the registered primary arm, and that the level and
dispersion arms are not yet reported. General rule: a dependence is not declared removed until the
registered arms that measure it are reported.

### E42. A real-data test described as model-independent without an arm for the assumption that the real data contradict
**Text:** Section 16 and Table 8 of the proposal: G1' "needs no waveform-grade velocity model".
**Symptom:** the radiation operator of the primary detector uses amplitudes and polarities of P and S on three
components and imposes P along the ray and transverse S at all stations, and the only real measurement
(Appendix T.1) gives reliable polarisation at 1 of 5 stations of the RSNC. The risk table treated that
fact only for C-IV, not for MT-FW, and there was no arm of G1' that addressed it.
**Cause:** the measurement of T.1 was made for the per-pick polarisation channel and was not carried over to the operator
of the pick-free detector, which makes the same physical assumption.
**Detection:** fourth independent audit (weakness 6 and inconsistency 12).
**Scope:** it does not affect any synthetic result; it affects the interpretation of the decisive test, which
could have failed because of an unaudited assumption without the design making it possible to know.
**Correction:** Section 16 declares the dependence and registers two arms with G1': one whose
moment-tensor subspace uses only the stations that pass the orientation and polarisation audit, which is not
run where fewer than two pass, and a sign-blind fallback one, which uses only the norm of the peak vector
of each station and phase. Table 7 adds the risk and Appendix D the arms; the decision of G1' is still
taken on the complete operator. General rule: every unfavourable real measurement is carried over to every
component that shares the assumption, not only to the one that motivated it.

### E43. An event-free twin that retained remnants of the events
**Script:** first version of `explore_contrib/c13_twin.py`, exploration C13 (seed 999).
**Symptom:** the maximum scores of the twin windows, which should have had the same clutter and the same
noise as the original window without the event part, did not match those of empty windows of the
same network.
**Cause:** the twin copied the random generators with `copy.deepcopy`, which loses the seed
sequence and its child counter; since the waveform generator spawns one child per window, the
regenerated window was not the same and the twin contained remnants of the events.
**Detection:** on 5 October at 18:18 (−05), during the analysis of C13, when comparing the maximum scores
of the twin with those of empty windows of the same network; declared in `REPORT_C13.txt`, section 2.
**Scope:** exploratory. The invalid outputs (`out/c13_twin_parts_BUGGY_v1` and the previous analysis
`out/c13_analysis_v1_5decoys.txt`) are kept as history and do not enter any number; no registered
result depends on them.
**Correction:** the twin is built with a clone that preserves the seed sequence, each window is
regenerated in full from clones and compared sample by sample with the original (assert in each part), and the
computation of the twins was repeated in full. General rule: a control that regenerates random data is validated
sample by sample against the original before being used, not only by its distribution.

### E44. Decoys by circular shift of the whole trace
**Script:** decoys D3 and D1o of `explore_contrib/c13_lib.py`, fixed before the main run of C13.
**Symptom:** in event-free windows with strong clutter, D3 and D1o were not exchangeable with the target: D3
with permutation p from 0.001 to 0.011 and between 1.5 and 1.75 times more targets than decoys in the tail; D1o
rejected in two strata of strong clutter (p 0.001-0.003) and in one of weak clutter (p 0.009).
**Cause:** the generator draws the clutter between 2 s and 36 s of the trace and the hypotheses read samples
95-451 (5.9-28.2 s); an arbitrary circular delay brings the clutter-free ends into the read zone, and the
decoy sees less clutter than the target.
**Detection:** the exchangeability test in event-free windows foreseen in the plan of C13.
**Scope:** exploratory. The numbers of D3 and D1o with strong clutter carry that bias and do not support any
conclusion; the conclusion that incoherent decoys underestimate the false events rests on D3c.
**Correction:** D3c, added after seeing the defect and declared as posterior, restricts the delays to
[−7.5 s, +3.5 s], separated by at least 1 s, and passes the test (p 0.35-0.82); D1o was not corrected for lack of
time. General rule: a data surrogate is built on the same time window that the
hypotheses read.

### E45. Two processes at once, a wrong denominator and a runner extended without detail, around the freeze of C14
**Script:** `explore_contrib/c14_run.py` and `c14_analyze.py`, exploration C14 (seed 999).
**Symptom:** (1) when launching the second development round two processes were left running at once and a
test of 1 + 1 windows wrote parts into `out/c14_dev2_parts`; (2) in `c14_analyze.py` the "far" ratio
of the twins used as denominator all the atoms of the twin, not those that are far from the
events; (3) the hash of `c14_run.py` at the freeze (58c0fe4d…, `C14_FREEZE.txt`) is not the one noted when
fixing the design on 5 October (7747fb9f…, `REPORT_C14.txt`, section 2).
**Cause:** (1) a relaunch without checking that no other process remained; (2) a code error; (3) the
runner was extended between the fixing of the design and the freeze, and the report does not detail that change.
**Detection:** (1) and (2) by whoever carried out the exploration, declared in `REPORT_C14.txt`, section 2; (3) when
integrating C14 into the documents, comparing the hash of the report with that of `C14_FREEZE.txt`.
**Scope:** exploratory; no number changes. (1) The two processes were stopped, the directory was deleted and
the second round was repeated in full with a single process; it affects neither the development nor the evaluation. (2) It was
corrected after the freeze and before reading any evaluation number, together with two declared
additions to the analysis (the TIGHT descriptor of close pairs and the table of atoms with shifted
origin), with the final hash in `C14_FREEZE.txt.post`. (3) The evaluation is fixed by the hashes of the
freeze (9 October, 11:36:11, −05), written before generating any window of networks
1610-1621; the content of the change to the runner cannot be reconstructed from the disk.
**Correction:** Appendix Y.29 of the supplementary volume declares the three things. General rule: a single
process per queue with a prior check, and every modification of a runner between the fixing of the design and
the freeze is noted with time and reason, even if it is the freeze that fixes the evaluation.

### E46. The synthetic benchmark underestimated the real timing error, and the registered detector did not transfer to real data
**Script:** `r0/r0_run.py` and `r0_analyze.py`, Task R0 (Ridgecrest 2019, real data).
**Symptom:** in the development on the first 15 minutes the registered detector hit 3 of 16 events
of M ≥ 2.0 at any score; in the registered evaluation the primary contrast was −0.383 [−0.481,
−0.310] against the best associator with PhaseNet.
**Cause:** it is not a code error (sample table, absolute origin and time alignment checked;
independent recomputation of the points of N8 and N6 identical). It is a design assumption: the synthetic
generator gave, after the catalog corrections, timing residuals of about 0.14 s, and the reading of
±0.25 s was set for that; in the real sequence the residuals after constant terms were 0.35 to
0.41 s and grow with the separation between sources, and the capped statistic saturates in a sequence that
fills all the windows.
**Detection:** in the development on CAL, declared in a registered amendment before the freeze
(`R0_DEV_LOG.txt`), and confirmed by the evaluation.
**Scope:** the conclusions G1 and G1-C become conditional on a timing model with errors smaller than
the reading tolerance; G1' on the RSNC remains registered.
**Correction:** Appendix Y.30; abstract and Section 1.6 of the proposal rewritten. General rule: calibrate the
timing-error level of the generator with real data before accepting a tolerance.

### E47. Interruptions and resumptions of R0 and C16, and flat padding in the noise of C16
**Script:** `r0/queue_eval.sh`, `queue_resume.sh`, `merge_parts.py`; `explore_contrib/c16_data.py`.
**Symptom:** (1) R0 was cut off by a container restart on 9 October at 21:07 (−05) with N8
complete; (2) C16 was cut off by an API limit on 9 October at 15:03, before running anything; (3) when
preparing Task 2M, the new quality control found that 394 of the 2037 noise slots of C16 (19 %)
contain flat stretches of telemetry padding (MURB, BADI and ATPC).
**Cause:** (1) and (2) limits of the environment; (3) C16 had no quality control for flat stretches or for digital
spikes.
**Detection:** (1) and (2) by the logs of the runs; (3) by the quality rules of Task 2M.
**Scope:** (1) R0 was resumed with the frozen code verified unchanged (sha256sum -c) and with outputs
per subset of methods merged without touching the analysis; the execution order changed, not the rules.
(2) None. (3) C16 is exploratory and is not reanalysed; part of the advantage of the cap over BOMP-PK with
real noise may come from instrumental defects and not from seismic transients; Task 2M applies the
quality control.
**Correction:** declared in Appendices Y.29c and Y.30 and in the registration.

### E48. A declared per-network failure stopped a whole lane of Task 2M
**Script:** `t2m_run.py`, `t2m_lib.py` (Task 2M, seed 3141).
**Symptom:** the REAL lane of seed 3141 ended after network 2909, with ladder_6 and unilateral_c40 of 2909 and all
networks of day 2010-02-11 not run, and the queue moved on to the descriptive GEN arm.
**Cause:** the registration states that a network whose catalogue-slot reserve runs out fails and is declared;
`Net2M.catalogue` raises `RuntimeError` in that case, but `t2m_run.py` does not catch it, so the first failure (ladder_4
2909, 47 draws) terminated the loop over networks. The validation of the 2M code ran seed 999 networks whose reserve
never ran out, so the path was untested.
**Detection:** by the log of the run, before any part was analysed.
**Scope:** none on the result. An incident note was appended to the log before resuming; the remaining networks were run
one at a time with the frozen code unchanged (`t2m_resume.sh`, `t2m_resume_each.sh`), a second network (ladder_4 2911,
50 draws) failed in the same way and was declared, and seed 3141 ran 49 of its 51 networks. The GEN arm of 3141 was
stopped and resumed after the REAL arm, within the registered budget rule.
**Correction:** declared in the log and in Appendix Y.31. General rule: a failure that the registration declares local
must be caught locally, and the failure path must be exercised in the validation.

### E49. Raw run outputs were not preserved outside the working container
**Script:** `r0/` (Task R0) and `explore_contrib/c16_run.py` (C16).
**Symptom:** after the original container was lost, the detection outputs of R0 (`r0/run/*.json`) and the part files of
C16 were not available: the pending R0 contrast M1 - M3 could not be computed from stored outputs, and the
bit-for-bit check of the Task 2M code against the C16 parts could not be repeated.
**Cause:** the repository kept the analysis outputs and the freeze files but not the raw outputs from which the analyses
are computed, and the container that held them was ephemeral.
**Detection:** when rebuilding the environment for Task 2M and for the closure of R0.
**Scope:** R0: M1 and M3 were re-run with the frozen code on the same data, declared post hoc before running, and the
re-run of M1 serves as a reproducibility check of the registered values. C16: the environment check of Task 2M was
reduced to the catalogue draws and window counts of six regenerated parts, declared as such in the registration.
**Correction:** the raw outputs of Task 2M (part files) and of the R0 closure are archived in the repository; the data
instructions in `data/README.md` were corrected (the station inventory, the catalogue and the auxiliary files of the
Alto Tiberina data are not in the folder the earlier text named).

### E50. An exploration inherited a truncated search and drew a false conclusion about the medium
**Script:** `explore_contrib/c18a_timing.py` (exploration C18a), through `station_corrections.peak_delays`.
**Symptom:** C18a concluded that no layered or fitted medium predicts the energy-peak readings of the Ridgecrest
development events better than the registered homogeneous medium, and the documents stated that the timing error was
a reading problem, not a velocity-model problem.
**Cause:** the readings were searched within +-1.0 s of the homogeneous prediction, the default of `peak_delays`,
while real S delays at several stations are 1.7 to 1.9 s; 42 per cent of the S readings and 24 per cent of the P
readings on N8 lay within 0.1 s of the edge of the search, so every model was compared on truncated readings.
**Detection:** the independent verification of R0 (defect D1, the same truncation in the station-term estimators);
then a check of the C18a readings against the search edge.
**Scope:** the C18a conclusion is withdrawn. With a +-2 s search the layered model predicts the readings better on
every subnetwork (C18i), and with terms read over 2 s, the layered medium and no cap the detector meets the registered
gate on development data (C18j). C18b to C18h remain valid as measurements of the configuration with the truncated
terms. No registered result is affected.
**Correction:** Appendix Y.32 rewritten; the documents corrected. General rule: before comparing models on readings,
check what fraction of the readings sits at the edge of the search window.

### E51. A matched filter changed to 3 Hz without changing the deflation kernel
**Script:** `explore_contrib/c18h_kernel.py` (`set_kernel`), used by C18h, C18i and C18j at 3 Hz.
**Symptom:** the 3 Hz detector readings of C18h, C18i and C18j moved by up to one event of sixteen when only the order
of two imports was changed.
**Cause:** `set_kernel` replaces the kernel of the matched filter, but the pulse shape used to deflate each accepted atom
(`task1x_mtglrt.AC`) is computed once when that module is first imported. C18h set the kernel before the import (a 2 Hz
pulse correlated with the 3 Hz kernel); C18i and C18j imported first (the 2 Hz shape throughout). None of the three
deflated with the response of the 3 Hz kernel.
**Detection:** the independent code audit of 17 October 2026 (finding F1, `results/code_audit/ac_order_bh.txt`).
**Scope:** the 3 Hz readings of Appendix Y.32 are not comparable with each other; the 2 Hz runs, including the
configuration that met the gate, and the single-atom oracle maps are not affected. No registered result is affected.
**Correction:** Appendix Y.32 marks the 3 Hz readings as affected; R1 uses the 2 Hz filter only.

### E52. Exact-time readers read the floor sample, not the nearest one
**Script:** `certified_fw.py` (`_idx_node`), `explore_contrib/c10_lib.py`, `c11_lib.py`, `c18j_layered.py`.
**Symptom:** a reader at the exact predicted time captured on average 0.78 of the energy of an arrival with no timing
error (worst case 0.42); the nearest sample gives 0.94 (worst 0.82).
**Cause:** the predicted sample is `int(t / DT)`, which truncates, so the read sample is on average half a sample
(31 ms) early, and the +-4-sample peak window is centred half a sample early.
**Detection:** the independent code audit of 17 October 2026 (finding F2, `results/code_audit/floor_bias_bh.txt`).
**Scope:** the comparators of Task 2C read at the exact sample (BOMP+SC, MTSTACK+SC, GL-FW+SC, the convex problem of
2C-b) and the MT-FW-EX ablation of 2D are handicapped, which inflates the measured value of the peak search; the
decisions are very probably unchanged (the timing error of the grid is already 0.1 to 0.3 s) but were not re-run.
Comparisons in which every method has the same +-4-sample search (2D BOMP-PK against MT-FW, 2F, 2I, 2M, R0, R1) are
not affected beyond a common shift of at most one sample.
**Correction:** declared in Appendix Y.30 next to the 2C and 2D results; the engines are left unchanged so that the
registered runs remain reproducible.

### E53. Station terms estimated on the events they are then used to detect
**Script:** `r0/r0_run.py` (calibration windows), `explore_contrib/c18j_layered.py` and the C18d to C18j readings.
**Symptom:** the development readings of C18 on the Ridgecrest calibration segment were optimistic by up to two events
of sixteen.
**Cause:** the station terms are medians over the same sixteen events that are then counted as detected, so each event
contributes to the predicted times used to find it.
**Detection:** the independent code audit of 17 October 2026 (finding F5, `results/code_audit/c18j_loeo_bh.txt`).
**Scope:** re-estimating the terms without the events each window can hit, the gate configuration still recovers 5 of
16 on U8 (gate met exactly at its threshold) and 2 instead of 4 on N8; removing every event inside the window gives 1
of 16 on U8. R0 and R1 estimate their terms on the calibration segment and evaluate on a separate segment and are not
affected; the R1 development readings on its own calibration segment are in-sample in the same way and are reported
as development only.
**Correction:** Appendix Y.32 and the R1 registration report the held-out readings and the sensitivity.

### E54. Station terms measured on energy peaks applied to onset picks
**Script:** `r0/r0_picks.py` (the +SC arms of GaMMA and PyOcto).
**Symptom:** the corrected S picks were pulled early by about 0.4 to 0.6 s.
**Cause:** the +SC terms are delays of the matched-filter energy peak relative to the catalogue prediction; PhaseNet
picks onsets, which precede the energy peak by a phase-dependent lag, so subtracting the term also subtracts that lag.
**Detection:** the independent code audit of 17 October 2026 (finding F4).
**Scope:** the descriptive +SC readings of R0; BESTPICK, which takes the best of the raw and corrected arms, and the
registered R0 decision are not affected. The hybrid development and R1 estimate the pick terms from pick residuals
(`hyb_dev2.pick_terms`), which is the correction; the R0-style +SC arms of R1 remain descriptive.
**Correction:** declared in the R0 section of Appendix Y.30 and in R1.

### E55. Calibration windows assume an origin at 6.5 s on a rounded window start
**Script:** `r0/r0_run.py`, `explore_contrib/c18a_timing.py`.
**Symptom:** none measurable; the catalogue origin sits up to 31 ms from 6.5 s of window time.
**Cause:** the window starts at the nearest sample of `t - 6.5`, while the origin passed to the estimators is 6.5 s.
**Detection:** the independent code audit of 17 October 2026 (finding F3).
**Scope:** an error common to all readings of an event, of standard deviation about 18 ms, against timing errors of
0.35 to 0.6 s. No result is affected.
**Correction:** declared.

### E56. A development run read forty seconds of the evaluation segment
**Script:** `code/r1/r1_lib.py` (`t_end`), Task R1.
**Symptom:** the amendment of R1 stated that no detector, picker or associator had read evaluation samples before the
freeze, while the development runs read the waveforms up to 880 s and the evaluation segment starts at 840 s.
**Cause:** the development read the calibration segment plus the 40 s window margin used everywhere in R0 and R1,
and the margin was not cut at the start of the evaluation segment; the docstring of `t_end` claims the opposite.
**Detection:** the independent verification of R1 (defect D1, `results/taskR1/verif_r1/r1_verification.txt`).
**Scope:** the development picks reach 853.5 s and the development lists contain the first evaluation target (843.8
s); two evaluation targets lie in [840, 880) s. No constant was set with those samples (the station terms use
calibration events only, the associator configuration was chosen on calibration origins only). With the evaluation
restricted to origins from 880 s the primary rule gives 22 against 21 of 52 targets, discordant 1 and 0, the same
conclusion. The registered result stands; the statement of the amendment was false.
**Correction:** declared in the result entry of the log, in Appendix Y.33 and in the thesis. Rule for future tests:
the development reads must end at the start of the evaluation segment, margin included.

### On the order of the three errors of Task 0
Three complete protocols, three different conclusions, and **in none of them was the cause the physics**. By
effect size: the fixed geometry (E9), the scale of the threshold (E8), and the small-sample bias
of the polarisation estimator, which was the one that received the most attention and turned out to be the smallest of the three.
Diagnosing first the most visible instead of the largest is a failure mode in itself.

---

## 4. Implementations deliberately left broken

### `a_i_tuned_attempt.py`
Tuned azimuth estimator. **Fails in 62-100% of the realisations.** Three causes
identified and not corrected, documented in Appendix U and in Task 2:

1. The sign of a covariance eigenvector is arbitrary, so one has to intersect **lines**, not
   directed rays. The `ts>0` test of the code is invalid.
2. The `len(inl)<3` gates and the support threshold are too strict at the clutter level
   tested.
3. The Lipschitz test of C-III holds between arrivals of the **same phase**, and it is applied to a set
   of inliers selected by azimuth that mixes P and S, rejecting valid configurations.

It is delivered broken because Task 2 consists precisely in repairing it, and the diagnosis is worth more than
code that half works.

### `superseded/coda_picker.py`, `superseded/coda2.py`
Test of H7b with an STA/LTA detector. **Inconclusive**: the detector caps at 62% detection on an
isolated event, so the variance between realisations swallows any coda effect. No
number from these scripts was used. Conceptually replaced by `d1_coupling.py`, which
parameterises the loss rate instead of simulating the detector.

### `superseded/b3_local.py`
Exceeded the execution time with 49 locations × 5 replicates. Replaced by `b3b.py` with a reduced grid
and fewer replicates.

### `superseded/octree.py` to `octree4.py`
Iterations of the adaptive search. `octree.py` and `octree2.py` hung because the purely
geometric criterion refines $8^{\text{nivel}}$ cells with nothing to stop it. `octree3.py` vectorises
the precomputation. `octree4.py` introduces score-based branching. The final result is in
`octree_hypothesis_space.py`, and the stopping rule **never fires**.

### `superseded/overlap.py` to `overlap3.py`
Iterations of the overlap test. `overlap.py` uses a metric that cannot see the failure,
because it counts false hypotheses that *exceed* the true one when the real failure mode is
**ties**. `overlap2.py` counts ties. `overlap3.py` brings the sources closer, because in the
previous versions they were 24 km apart and their arrivals never overlapped.

---

## 5. The checks that would have caught the errors

1. Verify that the constraint or the channel **can act** on a case built for it to act.
   It would have caught E3.
2. Verify that the sampling resolves the scale of the effect. E2 and E4.
3. Verify that the failure value of the estimator **does not coincide** with the correct answer. E5.
4. Apply every physical test to **observations**, not to quantities predicted by the model.
   E3, and the same pattern in `a3_degradation.py`.
5. Check that the endpoint has variance in the regime tested, and report which one saturated. E1, and the
   saturations of Appendices N and S.

A sixth: **distrust any number that repeats exactly across different
configurations.** The 1993 identical evaluations on two different networks revealed that the stopping
rule of the octree was never triggered.

An eighth, from E8: **never threshold two scores of different scales with an absolute value.**
Ask first whether the maximum of the statistic is the same in all realisations; if it is not,
normalise.

A ninth, from E9, and it is the most general of all: **identify which is the largest source of variance and
average over that one.** Averaging over what is convenient instead of over what is dominant produces results
that do not replicate, and the non-replication will be mistaken for method noise.

A tenth, from E10: **before setting a threshold on a fraction, compute how much that fraction is worth
by chance in the real geometry.** A threshold below chance cannot be crossed, and one barely
above it measures the geometry and not the phenomenon.

An eleventh, from E11 and E12: **every gate needs a minimum sample and three outcomes, not two**:
pass, fail, inconclusive. And the reference object of a measurement must be the one the physics
points to, not the most convenient to compute.

A twelfth, from E15 and E16: **always compare against the best counterfactual without the channel, not against
the simplest.** If the counterfactual can tie, break the tie; if the generator makes trivial the
information that the channel competes to replace, degrade it to the real level before measuring the
contribution of the channel.

A fourteenth, from E19: **before concluding that a channel contributes nothing, verify that the generator does not
give the rival channel impossible information.** Negative results deserve the same audit as
positive ones; a result that confirms what one fears is not more reliable for that reason.

A thirteenth, from E17: **a threshold that cannot operate on an empty window is not a
detection threshold.** Before adopting a normalisation, ask what it emits in a window without signal.

A fifteenth, from E21: **every positive needs a negative control with the same degrees of
freedom and none of the physics**, and it is read at the realised false-event rate, not at the nominal
budget. If the control also wins, the gain is flexibility.

A sixteenth, from E22: **run the baseline twice on the same inputs from two different
processes before comparing yourself with it.** A published method may depend on a global state.

A seventeenth, from E23: **plot the complete curve before reading it at a point.** If between two
consecutive points there is a jump larger than the difference one wants to detect, the reading at that
point is manufactured by the sampling.

An eighteenth, from E24: **ask what the method assumes that the generator gives it for free.** If the
simulator produces exactly the direction, the time or the shape that the method predicts, and the real data
do not, the positive belongs to the simulator. Replace the assumption by a realistic one before believing the advantage.

A nineteenth, from E25: **a negative control must destroy the specific information that is
attributed to the channel, not merely shuffle it.** If the control still wins over the reference, before
calling it flexibility verify that it does not retain part of the signal.

A twentieth, from E26: **every rule that the generator applies to the signal must also be applied to the
noise.** If the true arrivals go through a detector and the false ones do not, the comparison measures the
generator. Count how many picks of each type each method receives before comparing.

A twenty-first, from E27: **measure the pick model against a real picker on the same
waveform before using it to feed a rival.** A parametric law can be optimistic in recall and
in time and at the same time deliver more clutter; it is measured per regime, not on average.

A twenty-second, from E28: **every false-alarm rate guarantee is contrasted with the realised rate
in noise-only windows before being written**, and it is contrasted again when the noise or clutter
regime changes.

A twenty-third, from E29: **write the proof against the code, not against the idea.** Each
hypothesis of a cited theorem must be traceable to a line of the script; if it cannot, the guarantee
belongs to the idealised algorithm.

A twenty-fourth, from E30: **before claiming novelty, search for the antecedent by the formulation
(which unknown, which penalty, which algorithm), not only by the family of methods**, and count the
closest antecedent among the pending comparators until it has been run.

A twenty-fifth, from E31: **fix and record the BLAS kernel of each process (`OPENBLAS_CORETYPE`)
when a result must be reproduced bit for bit across sessions, and compare against what was saved with exact
equality.** The automatic choice of kernel can change when resuming a run, and only an exact
control reveals a rounding difference that is later amplified in a greedy search with
nearly tied candidates.

A twenty-sixth, from E32: **do not run inherited code without reviewing it against the preregistration and without running
its validation.** Check that every module it imports exists, that the environment checks read what
the process actually loaded and not the variable that was requested, and that the tolerances are reachable in
floating point.

A twenty-seventh, from E33: **do not present as the method that performs a formulation that you have not solved
on the same data as its approximation.** Sharing the certificate does not imply sharing the
performance: the non-convex ingredients of the approximation may be precisely the ones that perform.

A twenty-eighth, from E37: **an ablation says which element carries an effect, not why.** Before
writing the mechanism, remove the proposed cause in the generator and sweep the level of every fixed parameter;
if the effect survives that ablation, the mechanism was a different one.

A twenty-ninth, from E38: **review the matching rule every time the event density changes.**
A location-only rule suffices with isolated events and is lax when several events fall within the
tolerance; in that case a rule with origin time is also registered.

A thirtieth, from E39: **read how each parameter of the generator is drawn before giving it a physical name.**
A factor that is drawn in each window is a perturbation, not a property of the station; check what
persists between the catalog and the windows before planning to estimate it, and test the plan with an oracle that
knows the true value: if the oracle does not reach it, no estimate will.

And a seventh, which comes from E7: **never place the true answer on a node of the search
grid, nor at the origin of coordinates.** A zero in an error column almost always means
alignment, not accuracy. The practical corollary is to draw the truth from a continuous distribution and
measure the endpoint by continuous optimisation, not by grid maximum. The 1993 identical evaluations on two different networks were what
revealed that the stopping rule of the octree was never triggered.

### Verification checks of the components (Tasks 1p, 1v, 1z, Appendix AA, 2A to 2F, 2I and the exploration)

Each component of the method and each comparator was verified against an independent reference before being used in a
comparison. The complete reports are in `results/components/`; those of Task 2A, verified by an
independent verifier after the runs, in `results/task2a/task2a_report.txt`, with its scripts and
outputs in `results/task2a/verif_2a/` (suffix `_verif`); those of Task 2B, in
`results/task2b/task2b_report.txt` and `results/task2b/task2b_selfcheck.txt`; those of Task 2C, in
`results/task2c/task2c_report.txt`, `results/task2c/task2c_validation.txt` and `results/task2c/task2c_smoke.txt`;
those of Task 2D, in `results/task2d/task2d_report.txt`, `results/task2d/task2d_validation.txt` and
`results/task2d/task2d_smoke.txt`; those of Task 2E, in `results/task2e/task2e_report.txt` and
`results/task2e/task2e_validation.txt`, with its independent verification in
`results/task2e/verif_2e/task2e_verification.txt` and its scripts in `results/task2e/verif_2e/`; those of the
exploration, in `results/explorations/out/`; and those of Task 2F, in `results/task2f/task2f_report.txt` and
`results/task2f/task2f_control.txt`, with its independent verification in
`results/task2f/verif_2f/task2f_verification.txt` and its scripts in `results/task2f/verif_2f/`; and those of
Task 2I, in `results/task2i/analyze_2i_primary_out.txt` (primary arm), `results/task2i/analyze_2i_out.txt`
(complete task, with the secondary arms) and `results/task2i/task2i_control.txt`, with its independent verification in
`results/task2i/verif_2i/task2i_verification.txt` and its scripts in `results/task2i/verif_2i/`.

**Task 2L (hypothesis G1-C, Appendix Y.28).** The results are in `results/task2l/analyze_2l_out.txt` and
`results/task2l/task2l_control.txt`, with the code freeze in `results/task2l/task2l_code_hashes.txt`
(37 files) and `results/task2l/task2l_freeze_time.txt`, and the independent verification in
`results/task2l/verif_2l/task2l_verification.txt` with its scripts (suffix `_verif`); the exploration C11 that
originated it is in `results/explorations/REPORT_C11.txt` and `explore_contrib/C11_PLAN.txt`. The task was registered on 4
October 2026 at 21:10:33 UTC (hash 13a6e9b1...a62a) on seeds 7171 and 9393, never used, and
compares MT-FW with cap 6 and refinement T3 with the station terms S2 (joint determination of
hypocentres and terms, without picks, on the 20 catalog events) against the registered corrections.
Rules 2L-a to 2L-d and the secondary on BOMP-PK hold in both seeds: +0.057 and +0.062 in standard
windows and +0.038 and +0.030 in dense ones with the hit rule that requires the origin time. The verification
reproduced the rules with its own code, recomputed S2 without access to the true delays or catalog with an
identical result in 300 of 300 units, regenerated 6 identical parts and checked the 37 hashes. Its scope:
the gain is essentially that of the one-sided 8-station cell (+0.108 and +0.113), there is no
demonstrated gain with 4 stations, S2 worsens the terms in about 22 % of the networks, the effect is conditional
on T3 and assumes catalog errors of known size. Neither the report of C11, nor the registration, nor the
verification records an implementation error, and this log adds no entry. Two minor
things remain noted: the analysis script of C11 was changed once after the freeze (three lines that
omit a ratio whose methods were not run; the frozen file is kept as
`explore_contrib/c11_analyze_frozen.py.txt`), and the gauge of S2 is "median origin shift of the
relocated events equal to zero", declared in the plan before the freeze, not zero mean per phase.

**Explorations C12, C13 and C14 (Appendix Y.29).** The reports are in `results/explorations/REPORT_C12.txt`,
`REPORT_C13.txt` and `REPORT_C14.txt`, and the output summaries in `results/explorations/out/`. C12: with the
Frobenius metric the constants of the tube formula come out exact (L3 = 4π², L1 = 0); the constrained
statistic matches the explicit fit with rounding errors of 4e-15 to 7e-15; the optimiser
reaches a dense search of 300000 rotations (differences between −0.003 and +0.008); an exact double couple
without noise gives a ratio of 1.000000; the loop with the complete statistic reproduces `c10_lib.mt_loop` in 6
of 6 standard windows and 6 of 6 dense ones; and the pruned search of the maximum matches the exhaustive one in 8 of 8
windows (`out/c12_validation.txt`). C13: the decoy engines preserve orthonormal Q and the rank of G; the
sample table matches that of the C10 engine; the twins and D3c exactly reproduce the saved target
atoms; and in noise-only windows the six decoys give the same distribution as the target.
C14: with an empty configuration the loop reproduces `c10_lib.mt_loop` with position, origin and statistic
identical in bits in the 360 development windows and in the first two windows of each type of each
part of the second round and of the evaluation, and the decomposition of the statistic sums exactly (assert).
The incidents of C13 and C14 are recorded in E43, E44 and E45; C12 has none: it looked once at the
analysis of two networks to test the code, without taking decisions with it, and shared the machine with C13,
which only affects its times.

| Check | Script | Result |
|:--|:--|:--|
| Self-test of the realistic generator | `gen_realistic.py` | With g → 0 the time matches r/v0 (relative error ≤ 2.3e-7) and the ray matches the straight one (1 − cos ≤ 3.5e-13); asinh form against arccosh 1.9e-15; numerical ray tracing (DOP853, 3 pairs, P and S) ≤ 2.0e-13 in time and ≤ 1.2e-6 degrees in angle; S transport: norm 5.6e-16, transversality 1.7e-16, radiation with g = 0 equal to `generate_y` 1.8e-15; Hilbert transform via Dawson 1.8e-9; coda/direct ratio 0.301 (requested 0.3); with g = 0 and VS0 = 3 km/s times and amplitudes of all picks equal to `generate_y` (differences 0 and 3e-16) |
| External control of the generator | eikonal by fast marching (scikit-fmm), outside the package | Relative time error of 4.1e-5, 1.9e-5 and 2.3e-5 with a step of 0.02 km, convergent with the step; criterion < 1 % |
| Calibration of the pickers on seed 999 only | `real_picker.py --mode calib` / `--mode report` | Calibration on 80 windows (evbase 1000, 840 arrivals), report on another 160 (evbase 25000, 1608 arrivals); seeds 2024 and 7331 were never used; optimal thresholds in the interior of their grids; in-sample and out-of-sample recall of the same order (0.061/0.148/0.404 against 0.041/0.143/0.420 for PhaseNet, STA/LTA and degrade_y) |
| PhaseNet chain | `real_picker.py` | On labelled samples from the repository the maxima fall 0-10 samples from the arrivals (P 0.94-0.98); without noise, 0.8-0.9 at almost all arrivals of `generate_real` |
| Wiring of QuakeMigrate | `qm_wrapper.py` | Table times equal to d/V; grid coordinates equal to the local km within 1e-12 km; in a window without clutter, table and data times exact (S00 HHN with a peak at 13.00 s, expected P at 13.01 s); the in-memory trigger and the file-based `Trigger.trigger` give identical events in 192 window-configuration pairs |
| Positive controls of QuakeMigrate | `qm_wrapper.py` with `task0_common.generate` | One source, without clutter, 8 windows per noise level: default channels 3/8 and 4/8 events found (median epicentral error 0.8 and 15.4 km); three components 8/8 and 7/8 (0.5 and 0.5 km). The failure of the default channels is due to the physics of the generator (incidence of 60 to 80 degrees from the vertical), not to the wiring. `locate` on 5 windows does not change any hit decision |
| Identity of the certified oracle | `certified_fw.py` | Records identical bit for bit to `fw_record` in 72 of 72 windows and in the three modes (Lipschitz, box, envelope); 0 violations of the bound in 48 residual fields × 300 cells, raw and deflated, and in 1,004,400 hypotheses; maximum slack of the envelope 3.5e-9 |
| Assertions of the station corrections | `station_corrections.py` | The certified oracle with corrections matches the dense loop with the same corrections (`dense_fw_record_c`) and, without correction, `fw_record`; the vectorised reading of the recall matches `exact_HF` and `recall_at` in each network; ablation of independent P and S search with the same final medians |
| EEC machinery | `rft_threshold.py` (`selftest_eec`) | Simulated smooth chi-square field on a 10^4 grid: EEC 7.007/4.664/1.874 against mean Euler characteristic 6.98/4.62/1.70 at u = 12/16/20; `fast_field` against `fields_x` 5e-14 |
| Cell bound of Appendix AA | `chk_aa_bound.py` | Proposition AA.3 on the implemented statistic in 480 random cells of four windows (one-sided 8-station and random 4-station), raw and after one deflation: 0 violations, minimum slack 28.4 |
| EC densities of Appendix AA | `chk_aa_ec.py` | Symbolic ρ1 to ρ4 equal to Worsley's general formula for ν = 2 to 10 and to the Gaussian identity with ν = 1 |
| 2A: conformity with the preregistration | reading of `task2a_confirm.py` and `analyze_2a.py` against `PREREGISTRO_1j_1k.txt` | All elements of the design present; declared deviations: catalog per network and regime instead of per network (cause of the failure of 2A-d with strong clutter, neutral between methods), four parts of seed 999 from an earlier version of the code, `gen_realistic.py` edited after the preregistration to add the delays argument that the design requires, ambiguous reading of 2A-e (both programmed before the evaluation data) |
| 2A: leakage of the truth | reading of code; regenerated networks | The truth enters MT-FW, WAVE-PT and QuakeMigrate only when labelling hits; the SC come only from the catalog windows (random streams independent of those of the evaluation); the points of 2A-e depend only on seed 999, written before the first evaluation file; no network of seed 999 coincides with an evaluation one (minimum distance between stations 0.34 km, no identical networks or delays) |
| 2A: identical hit rule | `score_events`, `certified_fw_record`, `greedy_record`, `qm_records`; `verif_2a/qm_trigwin_verif.py` | Nearest true source within 5 epicentral km, without time tolerance, a second event on an already hit source counts neither as a hit nor as a false event, for all; GaMMA and PyOcto without cap or separation, that is, scored more leniently; widening the trigger window of QM-3C+SC by 0.5 s per side changes its recall by between −0.03 and +0.01 |
| 2A: sign of the corrections | `verif_2a/starve_sc_verif.py`, 80 regenerated windows | With the true delays as oracle, recall at 0.1/0.25 of MT-FW 0.25/0.36 → 0.51/0.54, QM-3C 0.17/0.24 → 0.36/0.36, GaMMA-DY 0.29/0.36 → 0.41/0.41; the oracle correction never worsens; regenerated outputs equal to the saved ones in 396 of 400 checks (the remaining 4, QuakeMigrate scores in the fifth decimal with the same hits); self-test: corrected table equal to the table plus rint(corr × 50) in all cells |
| 2A: starvation of GaMMA and PyOcto | `verif_2a/starve_sc_verif.py`, `analyze_starve_verif.py` | The frozen STA/LTA gives 1.53 correct picks within 0.3 s out of 11.5 arrivals per event; only 20 % of the events have the 3 that GaMMA needs; median pick error −0.009 s, deviation 0.13 s; positive controls with true-time picks and the same wiring: GaMMA 0.88 (0.90 with oracle SC) and PyOcto 0.69 at 0.1; GaMMA with only the correct picks of the STA/LTA, 0.08 |
| 2A: independent rederivation | `verif_2a/rederive_verif.py`, `rederive_2Ae_verif.py` | Own FROC curve on the JSON files equal to Table 1 of `analyze_2a_out.txt` in the checked cells; frozen thresholds and realised rates of 2A-e equal to Table 4; the aggregated JSON is the concatenation of the per-network files |
| 2A: internal identities | `task2a_confirm.py` | Certified oracle, with and without SC, identical to the dense loop and to `fw_record` in 2 windows of each network (the 90 per-network files report 2/2/2/1 checks); the control MT-FW-NODE+SC without violations of the 9 km separation in 810,000 node pairs |
| 2A: parts of seed 999 from an earlier version | `verif_2a/rerun_part_verif.py` | Network 1 repeated with the final version: SC different by 9e-16 s; MT-FW (three variants), WAVE-PT, GaMMA, PyOcto and picks identical; QuakeMigrate scores different in the fifth decimal and one extra false event of QM-3C+SC (score 1.39989) below the frozen threshold at 0.1 (1.40042) |
| 2A: clustering by network (exploratory) | `verif_2a/cluster_boot_verif.py` | Bootstrap by network, 2000 replicates: the 12 cells above 0 at cap 0.25 (minimum lower bound +0.03); four include 0 at cap 0.1 |
| 2A: corrections with strong clutter (exploratory) | `verif_2a/sc_strong_verif.py` | Seed 7331: rms residual 0.22/0.27/0.26 s with the strong-clutter catalog against 0.12/0.14/0.15 s with a weak-clutter one; with the latter or with the true delays the gain of the SC returns |
| 2A: windows with 0 or 1 source (exploratory) | `verif_2a/subset_k01_verif.py` | MT-FW+SC outperforms QM-3C, QM-3C+SC and GaMMA-STALTA+SC in the 12 cells at both caps (minimum margins +0.06 to +0.07 at 0.1); GaMMA-DY+SC outperforms MT-FW+SC at 0.1 in 5 of 6 weak-clutter cells |
| 2A: self-test and determinism | `task2a_confirm.py selfcheck` | Passes; floating-point non-determinism of 1e-14 in MT-FW and 1e-5 in the coalescence of QuakeMigrate, without changing any hit label or any point of the FROC curve |
| 2B: identity operator | `task2b_ablation.py selfcheck` (seed 2024, weak, 6 stations, network 0, window 0) | Shapes v (31, 2, 6, 3), G (36, 6), I (36, 36) and (I @ vf) equal to v; score of FREE equal to ‖v‖²/2 in the 675 × 31 hypotheses (relative 6.8e-16); first score of the record equal to the maximum (relative 1.9e-16); I lstsq(I, vf) − vf = 0 exactly |
| 2B: sign invariance | `task2b_ablation.py selfcheck` | Changing the sign of the vector of one station leaves the score equal bit for bit; changing the sign of the waveform of any of the 6 stations leaves the records of FREE identical (6 of 6) and changes those of MT-FW+SC (6 of 6); in that window MT-FW+SC loses the event and FREE detects it |
| 2B: exact regeneration of 2A | `task2b_ablation.py run`, reaffirmed in `analyze_2b.py` | Stations, delays, number of draws, catalog locations and times and corrections identical bit for bit to `task2a_parts` (60 of 60 networks); truth and origin times identical (1200 of 1200 windows); MT-FW+SC recomputed in windows 0 and 1 of each network EXACTLY equal to what was saved (120 of 120), with the BLAS kernel of each network fixed (E31) |
| 2B: certified oracle without operator | `task2b_ablation.py run` | Certified FREE equal to the dense loop (`dense_fw_record_c`) exactly in 120 of 120 windows |
| 2B: pairing with 2A | `analyze_2b.py` | k, network, truth, origin times and records of MT-FW+SC identical between the JSON files of 2A and 2B; the intervals of 2A-b recomputed with the same resamples reproduce `analyze_2a_out.txt` |
| 2B: sensitivity to the BLAS kernel (exploratory) | `task2b_kernel_check.py` | 18 networks, 360 windows, with the default kernel: MT-FW+SC differs by at most 6e-15 relative without change of hits or order; FREE changes hits or order in 2 of 360 windows; recall at 0.1 and 0.25 identical in the 4 affected cells |
| 2C: review of the inherited code against the preregistration | reading of `task2c_methods.py` and `task2c_run.py` against `PREREGISTRO_1j_1k.txt` | Dictionary, Gram with the autocorrelation, selection, refit and score of BOMP, MTSTACK, certificate of the group lasso and gap formula derived anew and conforming; three defects that prevented running, corrected before any run of 2024 or 7331 (E32); deviations declared in `task2c_report.txt`, section 10 |
| 2C: dictionary identities | `task2c_validate.py` (V1), 20 blocks per cell with pairs overlapping in time | Against explicit atoms in the data domain: max \|Aᵀy − h\|/max\|h\| ≤ 8.7e-16, max \|AᵀA − K\|/max\|K\| ≤ 7.3e-16; block Gram equal to the cross Gram exactly; rank of each K_bb 6 of 6 |
| 2C: recovery of a noise-free source | `task2c_validate.py` (V2) | In the model: block recovered 12 of 12, \|u − m\|/\|m\| ≤ 6.5e-13, second score / first ≤ 3e-30, first block of GL-FW true 12 of 12; physical waveform: node and origin 12 of 12, cosine(u, m) ≥ 0.978, 0.982 and 0.874 (4, 6 and one-sided 8); with the arrivals at the predicted sample, 1 − cosine ≤ 2.2e-16 |
| 2C: certified optimality of GL-FW | `task2c_validate.py` (V3) and the 1200 windows of part (i) | 24 windows with 2 λ beyond the stop: certificate off the support / λ = 1.000000000 (tolerance 1 + 1e-6), \|c/λ − 1\| on the support ≤ 5.1e-10, relative gap ≤ 1.5e-9, record with early stopping equal to that of the longer trajectory 24 of 24, λ_max²/2 equal to the first score of BOMP to 2e-16; in the 1200 windows, stop by complete record in all, maximum relative gap 1.4e-8, deviation on the support ≤ 8.7e-10, 0 re-solves |
| 2C: MTSTACK | `task2c_validate.py` (V4) | Equal to ‖Qᵀv‖²/2 at the exact times in 3 windows × 675 nodes × 31 origins per regime and cell (relative ≤ 1.5e-13); first score of the record equal to the maximum of the field |
| 2C: exact regeneration of 2A | `task2c_run.py run --part 1`, reaffirmed in `analyze_2c.py` | Stations, delays, draws and catalog locations and corrections identical bit for bit (60 of 60 networks); truth, origin times and k identical (1200 of 1200 windows); MT-FW+SC recomputed EXACTLY equal to what was saved in windows 0 and 1 of each network (120 of 120); 18 networks with SkylakeX kernel and 42 with Prescott, checked in each network |
| 2C: part (ii) | `task2c_run.py run --part 2` and `smoke` | Certified MT-FW+SC equal to the dense one exactly in 600 of 600 windows; in the smoke test, the pipeline of part (ii) reproduces the saved records of MT-FW+SC and QM-3C+SC of 2A in windows 0 and 1 (6 of 6 cells); `task2c_CONTROL_FAILED.txt` was not written |
| 2C: collapse of the exact-time methods (exploratory) | `task2c_explore.py` | Ideal generator, 4 stations, at 0.25: MT-FW 0.712, BOMP 0.458, GL-FW 0.237; with the sources on the nodes: MT-FW 0.763, BOMP 0.966, GL-FW 0.729; with oracle corrections in the windows of 2A, BOMP goes from 0.230, 0.222 and 0.380 to 0.264, 0.264 and 0.465: the dominant term is the grid |
| 2C: ties of GL-FW (exploratory) | `analyze_2c.py` | 1694 of 3000 (weak) and 1654 of 3000 (strong) events share a score with another one in their window; breaking ties by the group norm at entry changes the recall by at most 0.015 |
| 2D: order of the preregistration and of the code | `PREREGISTRO_hashes.txt`, `task2d_code_hashes.txt`, file dates | Entry of 2D with hash at 21:13:04 UTC on the 26th; code written from 02:25 to 02:29 UTC on the 27th; hashes of the three scripts, of `analyze_2d.py` and of the launcher at 02:41:59 UTC; budget note with hash at 02:41:50 UTC; first part file at 02:42:47 UTC and first of 2024/7331 at 02:53:04 UTC; result with hash at 03:54:55 UTC |
| 2D: reading at the peak | `task2d_validate.py` (V1) | Tabulated peaks equal to the `tpk` of `node_vectors_c` in all readings (0 discrepancies in 669 600 to 1 339 200 per cell), on the filter output and on a residual after one step of MT-FW+SC |
| 2D: exact reduction to 2C | `task2d_validate.py` (V2) | With the reading at the peak disabled, BOMP-PK, MTSTACK-PK and GL-FW-PK return exactly the records of BOMP+SC, MTSTACK+SC and GL-FW+SC (12 of 12 each); the loop with peak and cap is equal to `certified_fw_record` (12 of 12) |
| 2D: atoms with shifted times | `task2d_validate.py` (V3) | Explicit atoms with random shifts of ±4 samples, with overlapping pairs: max \|Aᵀy − h\|/max\|h\| ≤ 1.2e-15, max \|AᵀA − K\|/max\|K\| ≤ 6.4e-16, Gram of one atom to 3e-17, pseudo-inverse by subsets to 9e-13 |
| 2D: independent paths | `task2d_validate.py` (V4) | First score of MT-FW-NOCAP = maximum of the field of MTSTACK-PK (3.3e-15); of MT-FW-EX = vectorised maximum at the exact samples with cap (4e-16); of BOMP-PK = λ_max²/2 of GL-FW-PK (2.2e-16); orthogonality of the refit of BOMP-PK ≤ 7.9e-15 |
| 2D: optimality of GL-FW-PK | `task2d_validate.py` (V5) and the 1200 windows of part (i) | With 2 more λ: oracle certificate / λ = 1.000000000, \|c/λ − 1\| on the support ≤ 9.6e-10, relative gap ≤ 3.6e-9, record with stop = long record 12 of 12; in the 1200 windows, maximum gap 1.7e-8, 0 re-solves, 0 step caps; certificate valid only for the realised dictionary |
| 2D: exact regeneration of 2A | `task2d_run.py run --part 1`, reaffirmed in `analyze_2d.py` | Network, catalog, corrections and truth identical bit for bit (60 of 60 networks, 1200 of 1200 windows, against 2A and 2C); MT-FW+SC recomputed EXACTLY equal to what was saved (120 of 120); kernel fixed per process (SkylakeX in 18 networks, Prescott in 42) because the default one changed to Cooperlake, checked with threadpoolctl |
| 2D: calibration of part (ii) | `task2d_run.py run --part 2` | 30 calibration networks: stations and delays equal to 2A, identical corrections (maximum difference 0); the pipeline reproduces the records of MT-FW+SC and QM-3C+SC of 2A in windows 0 and 1 (60 of 60 each); `task2d_CONTROL_FAILED.txt` was not written |
| 2D: injected time error (exploratory) | `task2d_validate.py` (V6) and `task2d_explore.py e2` | ±0.2 s per station and phase, sources on nodes: the exact-time methods lose 0.16 to 0.48 and the PK variants −0.03 to +0.15, and these outperform their exact version in 12 of 12 cases (+0.04 to +0.38); without error the peak search costs 0.10 to 0.14 at 6 and 8 stations |
| 2D: false events in windows with events (exploratory) | `task2d_explore.py e1`, `task2d_explore_sum.py` | MT-FW+SC recomputed equal to 2A in 713 of 713 windows with events; with the conformal thresholds, 1.21 to 2.19 false events per window with events in weak clutter, 78 to 92 % of them 5-20 km from a true source, against 0.18 to 0.36 for QM-3C+SC |
| 2E: validation of the continuous atom | `task2e_validate.py` (V1 to V5) | Continuous identities ≤ 8.9e-16; gradients against central differences 1.4e-9 in 18 cases; convex step against an explicit dictionary (≤ 4.5e-11) and FISTA (3e-16), λ = 0 equal to least squares 5.0e-11; without refinement or sliding the first atom is that of BOMP-PK (24 of 24); V4 from the grid, 12 of 90 sources at < 0.5 km and < 0.02 s |
| 2E: order, hashes and controls | `PREREGISTRO_hashes.txt`, `task2e_code_hashes.txt`, `analyze_2e.py` | Entry of 2E with hash at 14:44:14 UTC; calibration note with hash at 15:38:20 UTC; hashes of the code at 15:38:31 UTC; first evaluation part at 15:39:40 UTC; `sha256sum -c` 5 of 5 after the run; 4800 of 4800 windows equal to what was saved; BOMP-PK+SC equal to 2D in 1200 of 1200 |
| 2E: independent verification, code | `verif_2e/check1_verif.py` | At a node with the shifts of the peaks, the continuous atom reproduces the block of BOMP-PK in 120 atoms (prediction 1.4e-15, Gram 9.0e-16, correlation 5.7e-16, objective 9.9e-16); objective against an independent implementation 2.3e-13; gradient 4.5e-8; injected shifts recovered to 1e-4 s |
| 2E: independent verification, numbers | `verif_2e/froc_rederive_verif.py`, `verif_2e/boot_verif.py` | Own FROC without importing the project: reproduces exactly Table 1 and the 24 points of 2E-a and 2E-b; own bootstrap with another random generator, intervals equal or almost equal |
| 2E: independent verification, budget | `verif_2e/check3_budget_verif.py` | 10 times the iterations and 9 starts in 270 windows (208 events): no call reaches the cap, the objective decreases by a median of 1.9 %, the recall at 0.25 changes by −0.014 [−0.052, +0.038] (6 events gained, 4 lost) |
| 2E: independent verification, ceiling | `verif_2e/check4_ceiling_verif.py`, `analyze_verif.py`, `anatomy_verif.py`, `reloc_verif.py` | CEIL-G, CEIL-GD and CEIL-C on 600 windows of seed 999 with the exact clutter and noise part (regeneration 600 of 600): BOMP-PK 0.16-0.37 under CEIL-G and 0.08-0.24 under CEIL-GD at 0.25; missed events absorbed by atoms at 5-10 km |
| 2E: independent verification, hash chain | `verif_2e/hashchain_verif.py` | The 6 hashes are increasing prefixes of the registration; the text of 2E is identical to the one hashed at 14:44:14 UTC; the calibration note precedes the first evaluation part |
| Exploration: harness and validations | `explore_contrib/t_check.py`, `c6_validate.py`, `c7_validate.py`, `c7_resolution.py` | Corrections equal to `setup_net`; `c6_lib` with the homogeneous model equal to the registered detectors (≤ 3e-13 in BOMP-PK, 8e-14 in MT-FW) and the ceiling equal to that of the verification of 2E (7.5e-5); `c7_lib` with q = 1 equal to the Gram and the correlations of `PKDict` of 2D exactly, and without a P/S cross term with the hits and locations of BOMP-PK in 2400 of 2400 windows (scores ≤ 1e-12) |
| 2F: controls on seed 999 | `task2f_control.py` | Dense pilot reproduced bit for bit (180 of 180 windows per regime, 5 methods, recall equal to the published one); MT-FW-NOCAP = BOMP-PK-SEQ in 360 of 360; standard records of C6 equal (60 of 60 per regime); fast FROC equal to the reference one (difference 0.0); pilot DE reproduced |
| 2F: order and freeze | `PREREGISTRO_hashes.txt`, `task2f_code_hashes.txt`, `task2f_freeze_time.txt`, `log2f.txt` | Registration with hash at 03:34:51 UTC; hashes of 26 modules at 03:49:08 UTC; single launch at 03:49:14 UTC; parts from 03:49:29 to 06:52:31 UTC; `sha256sum -c` 26 of 26 at 06:53:11 UTC; analysis once at 06:54:13 UTC |
| 2F: independent verification, provenance | `verif_2f/repro_verif.py`, `provenance_verif_out.txt` | Chain of 8 increasing prefixes; 3d8e4c9b ends at the end of the Task 2F section; no prior data of 4242/8888; 3 parts regenerated with BOMP-PK before MT-FW identical in 30 of 30 windows; MT-FW reimplemented bit for bit (216 of 216) |
| 2F: independent verification, scoring and statistics | `verif_2f/diag_verif.py`, `hitrule_verif.py`, `stats_verif.py` | Duplicates as false events, truncation to kmax 5/3/2 and k = 3, 4, 5 separately keep the gain; with origin time within 0.5 s, +0.024 to +0.045 (E38); own FROC and bootstrap reproduce 2F-a, 2F-b and 2F-c (8888: p = 0.047 against 0.0405 registered) |
| 2F: independent verification, mechanism (exploratory) | `verif_2f/mech_verif.py`, `mech_analyze_verif.py`, `capfine_verif.py`, `capfine_analyze_verif.py` | Sweep of the cap (peak at CAPV = 6, negative with 3 or less, null with 24); without clutter the gain rises to +0.137; without clutter or site factors it falls to +0.023 [−0.040, +0.069]; improvement mostly across windows (E37) |
| Exploration C8: validations and fact of the generator | `explore_contrib/c8_validate.py`, `c8_sitedev.py` | Copy of the generator equal to `generate_real` bit for bit (5 dense windows, difference 0); catalog and corrections equal to `xc.NetCfg`; engine equal to `c6_lib.mtfw_dense` with and without cap (6 of 6 and 6 of 6); Huber with infinite threshold equal to least squares (6e-14); bounded search equal to the maximum of the complete lattice (8 of 8); the amplitude factors are drawn per window (deviation of the log-factor across the windows of one network 0.42); catalog terms against the true factors in PERS, median correlation 0.72 (E39) |
| Explorations C9 and C10: validations | `explore_contrib/c9_standalone.py` and `c9_run.py` (validation mode), `c9_fisher.py`, `c10_run.py` (validation mode) | C9: vectorised operators equal to `task1x_mtglrt.radiation_operator` (3e-16) and rays equal to `gen_realistic._ray` (3e-16); noise-free recovery from 0 to 1 km with the consistent statistics (the 1/r fit is not consistent, median 6.6-15.6 km); `bomp_c9` equal to `xc.bomp_pk` and `mtfw_c9` equal to `c8_lib.mtfw_general` in 9 of 9 windows, repeated after the later edits; bound against Monte Carlo, maximum likelihood 1.2-1.4 times above; amplitude calibration 74.4 noise deviations per unit (analytical value 77). C10: loops in hard mode equal to `xc.bomp_pk` and to `c8_lib.mtfw_general` in 9 of 9 windows; with exact model and statics, source on a node located at 0.00 km by all variants; each later script regenerates the windows and checks that its MT-FW atoms are those of the main run |
| 2I: controls on seed 999 | `task2i_control.py` | Pilot C8 reproduced exactly (90 of 90 dense windows and 72 of 72 standard, 4 methods; 60 of 60 with the secondary thresholds; 24 of 24 with deviation 0.25 and 0.8); general engine with least squares equal to `c6_lib.mtfw_dense` with and without cap (60 of 60, relative difference of scores 2.2e-15); FROC equal to `xc.recall_caps` and to the published values; bootstrap on the pilot +0.064 [+0.027, +0.104] against +0.064 [+0.028, +0.105] published; estimate of 5.27 h on 2 cores, under the budget of 8 h |
| 2I: order, freeze and interruptions | `PREREGISTRO_hashes.txt`, `task2i_code_hashes.txt`, `task2i_freeze_time.txt`, `log2i.txt` | Registration with hash at 18:00:05 UTC; freeze at 18:14:44 UTC; first job at 18:14:45 UTC; parts from 18:14:58 to 22:04:25 UTC; interrupted at 20:33:25 (resumed 20:34:09) and at 20:46:15 UTC (resumed 21:29:10), each time from the parts and with a correct `sha256sum -c`; 1101, 83 and 256 parts per session; result in the registration at 22:08:07 UTC |
| 2I: independent verification, provenance and regenerated parts | `verif_2i/regen_verif.py`, `regen_all_verif.sh`, `regen_out_verif.txt` | Hashes of 27 files without failures; chain of 11 prefixes of the registration, a983beac ends at the end of the Task 2I section; no part prior to the freeze; internal and file times equal within 1.0 s in the 1440 parts; exactly two gaps of more than 60 s (82 s and 2599 s); no prior data of 5151/6262; 6 parts regenerated from scratch, 2 per session, identical; 260 windows of the `thr` arm with HUBER-N identical to that of the primary. Not verifiable from the disk: that the hash was sent before the run |
| 2I: independent verification, symmetry and statistics | `verif_2i/stats_verif.py`, `run_stats_verif.py`, `stats_verif_out.txt` | Same windows, grid, origins, corrections, kmax, exclusion and hit function for the 4 methods; hit labels recomputed without differences; point recall equal to 4 decimals; 2I-a, 2I-b and 2I-c reproduced with own loader, FROC and bootstrap (for example +0.0582 [+0.0402, +0.0753] against +0.058 [+0.039, +0.076]) |
| 2I: independent verification, scope (exploratory) | `verif_2i/run_explore_verif.py`, `explore_verif_out.txt` | With origin time within 0.5 s, +0.031 and +0.019 (the point criterion of 0.03 would fail in 6262); equivalence not met in standard windows with weak clutter of 6262 (−0.032 [−0.064, −0.003]); recall with all atoms accepted without increase (0.623 and 0.622; 0.628 and 0.616), gain from ranking and choice of atoms; −0.018 against the cap at 0.1 in 5151; 0.93 s against 0.10 s per dense window, 9 times MT-FW; low absolute recall. Secondary arms not verified |

---

## 6. How to reproduce

All scripts are self-contained except those that `exec` another one, and they depend only on `numpy` and
`scipy`. None needs a GPU, a network connection or external data. Typical times range between seconds and a few
minutes; `b3_local.py` and the first versions of `octree` do not finish.

The seeds are fixed with `np.random.default_rng(n)` in each file, so the numbers are
exactly reproducible. If when reproducing you obtain different numbers, first check the version of
numpy: `np.linalg.eigh` can differ in the sign of the eigenvectors between versions, and several
scripts depend on polarisation directions.

The scripts of the components and of Task 2A have additional dependencies: `real_picker.py` uses
ObsPy and, for PhaseNet, TensorFlow 2.15 with the original repository in `~/.cache/phasenet`;
`qm_wrapper.py` uses QuakeMigrate (`pip install quakemigrate`); `chk_aa_ec.py` uses SymPy; Task 2A
also uses GaMMA from `gamma_src/` and PyOcto. The commands and times are in
`results/components/README.txt` and, for Tasks 2A to 2F and the exploration, in `results/task2a/README.txt`,
`results/task2b/README.txt`, `results/task2c/README.txt`, `results/task2d/README.txt`,
`results/task2e/README.txt`, `results/task2f/README.txt`, `results/task2i/README.txt`, `results/task2l/README.txt` and `results/explorations/README.txt`. To reproduce bit for
bit the records of Task 2A, and the exact controls of Tasks 2B to 2F and 2I, one has to set
`OPENBLAS_CORETYPE` to the kernel with which each network was written (E31); Tasks 2C to 2F and 2I set it and
check it in each process (E32) and additionally need threadpoolctl.
