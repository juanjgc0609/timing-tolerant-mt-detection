# Timing-Tolerant, Bounded-Influence Moment-Tensor Detection of Overlapping Earthquakes

Undergraduate thesis (trabajo de grado) by Juan José Gordillo Córdoba. It studies whether detecting and locating
overlapping earthquakes on sparse seismic networks directly from three-component waveforms, before any phase picking,
with the physics of the source (a moment-tensor subspace), a timing tolerance and greedy deflation, does better than
pick-based association, and whether bounding the influence of each station helps. The main document is
`docs/thesis/Thesis.md`. The work was first planned in the format of a doctoral proposal; that document is kept in
`docs/proposal/` as the record of the original, broader plan, together with the evidence volume that holds every task
in full.

## Status in one paragraph

Eleven tests were registered before their run; which of them also had frozen code, a single run and a verification
by a separate AI session is listed in Table 3 of the thesis, and Section 10 of the thesis states how AI tools were
used. On a synthetic benchmark, timing-tolerant moment-tensor detection exceeded GaMMA, PyOcto and QuakeMigrate fed
by an STA/LTA picker in Task 2A (with operating points chosen on the evaluation data), but block orthogonal matching
pursuit with the same timing tolerance matched the thesis detector (2D) and an off-grid convex variant was refuted
(2E). Bounding the influence of each station (a cap or a Huber fit) added 0.04 to 0.07 in recall in dense synthetic
sequences (2F, 2I) and +0.228 [+0.148, +0.298] (epicentre rule; +0.085 with origin time) on one day of real noise with transients and
injected events (2M); the gain belongs to the
per-station bound, which block pursuit can also use. On real seismicity the pick-free detector failed: at most 0.022
of the targets on one hour of the 2019 Ridgecrest sequence against up to 0.667 for PhaseNet with PyOcto (R0, -0.383
[-0.481, -0.310]), and 1 of 54 on a 2012 Hi-net swarm (R1). A hybrid that adds pick-free candidates to a PhaseNet
associator and scores them by the picks they explain found 23 of 54 against 22 for the associator alone (R1, exact
McNemar p = 0.50): no gain demonstrated. A short theory of the deployed statistic (Appendix AA.9) shows why the cap
helps with transients and hurts in continuous seismicity. The contribution is the registered benchmark with its
positive and negative results, the bounded-influence result and its theory, and a measured account of why the
pick-free detector does not transfer to real sequences.

## Hypotheses and confirmatory results

| Task | Hypothesis | Result | Where |
|:--|:--|:--|:--|
| 2A-2D | G1: timing-tolerant MT detectors beat published associators and migration | Ahead of GaMMA, PyOcto and QuakeMigrate with an STA/LTA picker in all 12 cells of 2A; BOMP-PK with the same peak search equals MT-FW (2D-a not met); operating points chosen on evaluation | `results/task2a` to `task2d`; Y.18-Y.21 |
| 2E | G1-M: off-grid group BLASSO beats block OMP with equal tolerance | Refuted (+0.000 to +0.039, all intervals include 0) | `results/task2e`; Y.22 |
| 2F | G1-R: per-station cap in dense sequences | Confirmed, +0.043 to +0.074 | `results/task2f`; Y.24 |
| 2I | G1-R: Huber statistic with noise-fixed level | Confirmed, +0.058 and +0.053 | `results/task2i`; Y.26 |
| 2L | G1-C: pick-free station terms + local refinement | Confirmed, +0.057 and +0.062 (standard), +0.038 and +0.030 (dense), mainly 8-station one-sided networks | `results/task2l`; Y.28 |
| R0 | G1 on real data (Ridgecrest 2019) | Negative, -0.383 [-0.481, -0.310]; independently re-run | `results/taskR0`; Y.30 |
| 2M | G1-R with real noise of independent days (Alto Tiberina) | Confirmed on one primary day with transients, +0.228 (cap) and +0.241 (Huber) with the epicentre rule; +0.085 and +0.104 with origin time | `results/task2M`; Y.31 |
| R1 | Hybrid: pick-free candidates scored by PhaseNet picks vs the associator (Hi-net 2012) | Not met: 23 vs 22 of 54, p = 0.50; vs best configuration 23 vs 26 | `results/taskR1`; Y.33 |

Exploratory studies C1 to C16 (calibration seed 999 only) are in `results/explorations/` and Appendices Y.23, Y.25,
Y.27, Y.29 to Y.29c; exploration C18 (Ridgecrest development segment and seed 999) in `results/explorations/c18/` and
Appendix Y.32. An independent audit of the code is in `results/code_audit/`.

## How the project is run

The process (stages, gates, data partitions, registration with an external timestamp, verification) is in
`docs/protocol/RESEARCH_PROTOCOL.md`; the plan in `docs/ROADMAP.md`; the map from the original plan to what was
implemented, refuted or is pending in `docs/TRACEABILITY.md`; decisions in `docs/decisions/`; data access in
`docs/data_access.md`; changes in `CHANGELOG.md`.

## Repository layout

```
docs/
  thesis/Thesis.md                                    the thesis (main document); evidence_table.* its Table 1
  proposal/Doctoral_Proposal.md                       the original plan, in the format of a doctoral proposal
  proposal/Supplementary_Evidence_Appendices_G_to_AA.md   evidence volume: every task and exploration (Appendix Y)
  executive_summary.md                                decision document: what is measured, discarded and pending
  method_design.md                                    design of the method and its hypotheses
  tasks_and_real_data_plan.md                         every task: design, rules, results, next steps
  artifact_manifest_and_error_log.md                  artifacts and the register of errors E1-E58
  audit_status.md, audits/                            independent audits and the status after them
  task1_step_by_step.md
preregistration/                                      append-only registration log, hash chain, drafts
code/                                                 all scripts, in the original flat layout (see code/README.md)
  explore_contrib/                                    exploration code C1-C18 and their freeze files
  r0/, r1/                                            real-data tasks R0 and R1
  hyb/                                                hybrid development (Ridgecrest CAL) and the R1 rule simulation
results/                                              per task: README, analysis outputs, raw results
  explorations/                                       reports C7-C16, outputs of C18 (c18/)
data/                                                 how to obtain the real data, with checksums
```

## Reproducing

The pinned environment is listed in `requirements-frozen.txt`. The code was run with Python 3.11, NumPy 1.26.4, SciPy 1.17.1, ObsPy 1.5.1, scikit-learn 1.8.0, PyOcto 0.2.0 and
OpenBLAS pinned with `OPENBLAS_NUM_THREADS=1` and `OPENBLAS_CORETYPE=Prescott` (bit-for-bit reproducibility of the
frozen runs depends on the BLAS kernel). GaMMA must be cloned into `code/gamma_src/` and PhaseNet installed as
described in `data/README.md`. Every task script states its command line in its docstring; the analysis of a task is
run from its frozen outputs, for example `python code/r0/r0_analyze.py eval <dir with the run outputs> 2000`. Paths in
the scripts and freeze files refer to the original working directory `/home/claude/work`; place `code/` there or
adjust the path constants at the top of the scripts.

## Language

All documents are in English. Two kinds of text remain in Spanish on purpose: the sealed preregistration log in
`preregistration/original_sealed_es/`, whose hashes would break if it were translated (an English translation sits
next to it), and comments, docstrings and printed labels inside the code. Analysis outputs in `results/` are English
translations of what the scripts printed; numbers are unchanged, and running a script prints its original labels.

## Licensing

Code: MIT License (`LICENSE`). Documents, figures and written results: CC BY 4.0 (`LICENSE-docs`). Third-party data
and software keep their own licenses: SCSN data (via benchflow-ai/skillsbench, Apache-2.0), INGV data (via
avuan/PyMPA37, GPL-3.0), NIED Hi-net data (via Geolab-USTC/MatchLocate2, GPL-3.0 for the repository, NIED terms of use
for the waveforms), PhaseNet, GaMMA, PyOcto and QuakeMigrate.
