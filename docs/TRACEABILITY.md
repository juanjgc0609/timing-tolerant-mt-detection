# Traceability: original plan, what was done, what was refuted, what is pending

Status as of 18 October 2026. The original plan is `docs/proposal/Doctoral_Proposal.md`; the state of each
hypothesis there is its master table (Section 1.6). This file maps every component of that plan to its code, its
evidence and its status, and adds the components of the next model (Roadmap, phases 3 and 4). Status values:
**implemented** (code exists and was evaluated), **refuted** (a registered test failed it), **partial** (evaluated in a
reduced form), **not implemented**, **planned** (in the roadmap with a gate).

## 1. Hypotheses of the original plan

| ID | Component | Status | Evidence | Code |
|:--|:--|:--|:--|:--|
| G1 | Timing-tolerant moment-tensor detection with deflation adds detections over the best associator | refuted on real data; supported on the synthetic benchmark with an STA/LTA picker | Tasks 2A-2D (Y.18-Y.21), R0 (Y.30), R1 (Y.33) | `code/certified_fw.py`, `code/explore_contrib/c10_lib.py`, `code/r0/`, `code/r1/` |
| G1-M | Off-grid group BLASSO with per-station shifts beats block pursuit | refuted (synthetic) | Task 2E (Y.22) | `code/task2e_methods.py` |
| G1-R | Bounded influence per station adds detections in dense sequences | implemented; confirmed on synthetic data (2F, 2I) and on one day of real noise with injected events (2M) | Y.24, Y.26, Y.31; theory AA.9 | `code/explore_contrib/c10_lib.py`, `code/task2f_*`, `code/task2i_*` |
| G1-C | Pick-free station terms | implemented; confirmed once on synthetic data | Task 2L (Y.28) | `code/explore_contrib/c11_lib.py` |
| A1-A4 | Ordering and minimal set of constraint families; model independence | partial (preliminary, synthetic) | Appendices M, N | see artifact manifest |
| B1-B3 | Observability certificate; network design | partial (B1 numerical evidence; B2 not run; B3 preliminary) | Appendices N, Q, S | see artifact manifest |
| C1-C3 | Station-count crossover; geometric failure class; single station | partial (C1 numerical; C2 not run; C3 preliminary) | Appendices O, P | see artifact manifest |
| D1 | Pick-free value grows with picker miss rate | partial (numerical evidence) | Appendices R, W | see artifact manifest |
| E1-E5 | Registered nulls (resolution, VP/VS, overlap, degradation, rank) | partial (numerical evidence for the nulls) | Appendices G, H, K, L, V | see artifact manifest |
| F1 | Waveform gate beyond the minimal set | partial (oracle matched filter only) | Appendix O | |
| F2 | Optimal-transport misfit beats time-domain misfit | not implemented (conditional on F1) | | |
| F3-F5 | Learned wavefield representation, exact symmetries, compute break-even | not implemented (conditional on F1) | | |
| G2-G4 | Clutter scaling; calibrated uncertainty; silence channel | partial (G2 numerical; G3 not run; G4 preliminary) | Appendices H, J | |

## 2. Physical modelling of the original plan

| Component | Status | Note |
|:--|:--|:--|
| 3-D heterogeneous elastodynamic forward model, Green tensors | not implemented | The detector uses analytic far-field radiation and kinematic travel times |
| Adjoint sensitivities and constrained full-wave inversion | not implemented | |
| Differentiable or learned wavefield representation (PINN branch) | not implemented | Replaced in the roadmap by learned travel-time fields (eikonal), which address the measured timing failures; full-wave modelling kept behind a gate (ADR-0003) |
| Layered 1-D travel times (direct and head waves) | implemented | `code/explore_contrib/c18_tt.py`; used in C18j and R1 |
| Station terms: catalogue (+SC), pick-free (S2), pick residuals, anchored terms | implemented | E53, E54 mark their limits |

## 3. Components of the next model (planned)

| Component | Failure it addresses | Phase | Gate |
|:--|:--|:--|:--|
| Benchmark harness and failure atlas: core set PhaseNet and PhaseNO with GaMMA and PyOcto, EQNet, template matching, BOMP-PK, MT-FW; extended set EQTransformer, REAL, QuakeMigrate, MALMI | Unknown failure regimes of the comparators | 3 | Atlas reproduces published comparator numbers on their own example data; regimes and counting rule fixed beforehand (`experiments/exploratory/ATLAS-0/card.md`) |
| Evidence from picker probability fields instead of energy peaks | Readings on the wrong peak (R1: 36 % of P, 20 % of S off by more than 0.5 s) | 4 | Fewer wrong-peak readings and more origin-time hits on development data |
| Multi-station coupling with physics (travel time and moment-tensor consistency) | Single-station picking misses small and overlapping arrivals | 4 | Recall gain over PhaseNO, EQNet and PhaseNet in the atlas regimes. Travel-time coupling of learned features already exists (EQNet, Zhu et al. 2022); the untested part is the moment-tensor consistency of amplitudes and polarisation and the joint decomposition of overlapping events |
| Learned or eikonal 3-D travel times with source-dependent terms | Timing error of 0.35 to 0.41 s (R0, homogeneous medium, readings truncated at +-1 s; E50). With a +-2 s search the layered medium already predicts better (C18i) and meets the gate on development (C18j), but only exactly at threshold and not with every event of the window held out (E53) | 4 | Held-out reading fraction within tolerance improves over the layered medium of C18j, with terms estimated without the scored events |
| Adaptive bounded-influence statistic | Saturation of the cap in continuous seismicity | 4 | Keeps the 2M gain without the C18 loss |
| Joint multi-event refinement; deflation with empirical Green's functions | Coda ghosts and under-subtraction | 4 | Ghost detections reduced by half or more |
| False-alarm calibration with a deployable threshold | Oracle operating points | 4 | Realized rate within budget on unseen data |
| Full-wave modelling | Path effects not captured by the above | 4 (conditional) | Atlas shows residual path errors in a band where 3-D models predict waveforms |
