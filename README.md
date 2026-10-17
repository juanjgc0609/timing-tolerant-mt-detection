# Timing-Tolerant, Bounded-Influence Moment-Tensor Detection of Overlapping Earthquakes

Doctoral research project by Juan José Gordillo Córdoba. It studies whether detecting and locating overlapping
earthquakes on sparse seismic networks directly from three-component waveforms, before any phase picking, with the
physics of the source (a moment-tensor subspace) and greedy deflation, does better than the published pick-based
associators and migration locators. The repository holds the proposal, the preregistration log, all code, the
confirmatory and exploratory analyses with their outputs, and the first test on real data.

## Status in one paragraph

On synthetic data, in eight preregistered tasks (five of them, 2A, 2E, 2F, 2I and 2L, independently verified), the
timing-tolerant moment-tensor detectors beat GaMMA, PyOcto and QuakeMigrate by 0.09 to 0.49 in event recall at equal
realized false-event rates in Task 2A (Tasks 2C and 2D compared them with QuakeMigrate only);
bounding the influence of each station (a per-station cap or a noise-scaled Huber fit) added 0.04 to 0.07 in dense
sequences (G1-R); and station terms estimated without picks added about 0.06 in origin-time-matched recall (G1-C). An
exploration that injected synthetic events into real noise strengthened G1-R. The first registered test on real data,
one hour of the 2019 Ridgecrest aftershock sequence, was **negative**: the detector recalled at most 0.022 of the 45
catalogued events of magnitude 2 or more at 30 false detections per hour, against up to 0.667 for PhaseNet with
PyOcto (pooled difference -0.383 [-0.481, -0.310]), because a homogeneous medium with constant station terms leaves
timing errors of 0.35 to 0.41 s and the capped score saturates in continuous seismicity. A registered test on real
noise of independent days then confirmed G1-R (+0.228 and +0.241 over BOMP-PK on a day with transients, Task 2M),
and an independent re-run of R0 reproduced every registered number and found that its station terms were read within
1 s while real S delays reach 1.7 to 1.9 s; on the development segment, with terms read over 2 s, a published layered
medium and an uncapped score, the pick-free detector recalls 4 and 5 of 16 events with origin time, which meets a
registered gate for a new real test (drafted for a Hi-net swarm, not yet registered). The contribution is the
preregistered benchmark, its positive and negative results, and the measured conditions under which pick-free
physics helps.

## Hypotheses and confirmatory results

| Task | Hypothesis | Result | Where |
|:--|:--|:--|:--|
| 2A-2D | G1: timing-tolerant MT detectors beat published associators and migration | +0.09 to +0.49 in all 12 cells (synthetic; operating points chosen on evaluation) | `results/task2a` to `task2d`; Appendices Y.18-Y.21 |
| 2E | G1-M: off-grid group BLASSO beats block OMP with equal tolerance | Refuted (+0.000 to +0.039, all intervals include 0) | `results/task2e`; Y.22 |
| 2F | G1-R: per-station cap in dense sequences | Confirmed, +0.04 to +0.07 | `results/task2f`; Y.24 |
| 2I | G1-R: Huber statistic with noise-fixed level | Confirmed, +0.058 and +0.053 | `results/task2i`; Y.26 |
| 2L | G1-C: pick-free station terms + local refinement | Confirmed, +0.057 and +0.062 (isolated), +0.038 and +0.030 (dense), mainly 8-station one-sided networks | `results/task2l`; Y.28 |
| R0 | G1 transfer to real data (Ridgecrest 2019) | Negative, -0.383 [-0.481, -0.310]; independently re-run and confirmed; closure M1 - M3 +0.000 [-0.023, +0.021] | `results/taskR0`; Y.30 |
| 2M | G1-R with real noise on independent days (Alto Tiberina) | Confirmed on the primary day with transients, +0.228 [+0.148, +0.298] (cap) and +0.241 [+0.179, +0.293] (Huber); +0.057 and +0.120 on near-Gaussian days | `results/task2M`; Y.31 |

Exploratory studies C1 to C16 (calibration seed 999 only) are reported in `results/explorations/` and Appendices Y.23,
Y.25, Y.27, Y.29 to Y.29c; exploration C18, on the development segment of the Ridgecrest hour and seed 999, in
`results/explorations/c18/` and Appendix Y.32.

## Repository layout

```
docs/
  proposal/Doctoral_Proposal.md                       full proposal (sections 1-23, references)
  proposal/Supplementary_Evidence_Appendices_G_to_AA.md   evidence volume: every task and exploration (Appendix Y)
  executive_summary.md                                decision document: what is measured, discarded and pending
  method_design.md                                    design of the method and its hypotheses
  tasks_and_real_data_plan.md                         every task: design, rules, results, next steps
  artifact_manifest_and_error_log.md                  artifacts and the register of errors E1-E50
  audit_status.md, audits/                            four independent audits and the status after them
  task1_step_by_step.md
preregistration/                                      append-only registration log, hash chain, drafts
code/                                                 all scripts, in the original flat layout (see code/README.md)
  explore_contrib/                                    exploration code C1-C18 and their freeze files
  r0/                                                 real-data task R0
results/                                              per task: README, analysis outputs, raw results
  explorations/                                       reports C7-C16, outputs of C18 (c18/)
data/                                                 how to obtain the real data, with checksums
```

## Reproducing

The code was run with Python 3.11, NumPy 1.26.4, SciPy 1.17.1, ObsPy 1.5.1, scikit-learn 1.8.0, PyOcto 0.2.0 and
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

No license has been chosen yet for the code and documents of this project. Third-party data and software keep their
own licenses: SCSN data (via benchflow-ai/skillsbench, Apache-2.0), INGV data (via avuan/PyMPA37, GPL-3.0), PhaseNet,
GaMMA, PyOcto and QuakeMigrate.
