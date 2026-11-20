# Glossary: one definition per term (rule 27)

Started 20 November 2026. Every script, card, ADR, thesis passage and the preregistration uses these terms with these
meanings, and cites the module that computes them. A term with more than one meaning in use is listed with each
meaning and a name that keeps them apart; new text must use the qualified name. Changing a definition changes this
file and every use in the same PR.

## Signal-to-noise ratio (three meanings in use; always qualify)

| Name to use | Definition | Computed by | Used in |
|:--|:--|:--|:--|
| **SNR_blk** (per-block peak) | median over the strong (station, phase) blocks of the peak of the matched-filter output of the source, over the robust sd of the filtered noise of that block | `M6-JOINT/run.py` (`noise_mf`, source scaling) | M6-JOINT, M6-CAL, M6-THEORY, M6-BIAS ("SNR 12", "SNR 8") |
| **SNR_net** (network energy) | sqrt(E), E the whitened energy of the noise-free event over all blocks (the matched-filter SNR of the whole network) | Appendix AB, A.7; `src/ttmt/detectability.py` | theory (resolvable separation, A.5, A.7) |
| **SNR_med** (atlas cell) | median over the subset stations of the per-station SNR of the data card | `ttmt.atlas.cells.network_snr` | ATLAS-0 cells; FE-FIX population (`snr_med >= 5`) |

For the same scene SNR_net is about 2.7 times SNR_blk (supervision X2, 14 November 2026; M6-JOINT card, note on
"SNR 12"). Never set a number of one beside a number of another without converting and saying so.

## Regimes, endpoint and scoring (ADR-0008)

| Term | Definition | Computed by |
|:--|:--|:--|
| Regime **T** (target) | arrival separation to the nearest other event at the median station in [0.1, 1) s | `ttmt.atlas.cells` (separation), ADR-0008 item 1 |
| Regime **C** (complement) | every other event (separation < 0.1 s or >= 1 s) | same |
| **Recall T, recall C** | fraction of the regime's targets hit under LOC+T-V at 30 false detections per hour, exact FROC, operating point per method on the evaluated data | `ttmt.stats.label_voronoi`, `operating_point`, `hits_at` |
| **Objective** (comparator tuning) | mean of recall T and recall C on the tune partition | tuning ledgers (`st_retune_a6`, `harpa-a5`) |
| **LOC+T** | ADR-0007 one-to-one match in decreasing score: epicentre within 5 km and origin within 0.5 s | `ttmt.stats.cost_matrix` (TOL_KM 5.0, TOL_S 0.5) |
| **LOC+T-V** | LOC+T plus Voronoi eligibility: cost to the event <= 0.5 x cost to every other catalogued event | `ttmt.stats.label_voronoi` (DOMINANCE 0.5) |
| **Duplicator** | negative control: every detection of the best comparator emitted twice at origin +- half the cell separation | `ttmt.stats.duplicator` |
| **False detection** | a detection with no eligible event under LOC+T-V | `ttmt.stats.FALSE`, `false_count` |

Two scoring modules exist. `ttmt.stats` holds the ADR-0008 rule (LOC+T-V, duplicator, operating point, tests) and
is imported by `ttmt.v2.m6`, `ttmt.atlas`, M6-CAL, M6A-REALNOISE, STAGE-LOSS and P6-SAP. `ttmt.scoring` holds the
ports of the legacy FROC (`c10_lib.label`, `froc_exact`; `match="legacy"` has no Voronoi condition) and the one-to-one
matcher that `ttmt.stats` itself imports; it is used by `ttmt.atlas.froc`, `ttmt.atlas.instrument` and
`code/froc_strict.py`. A new script scores with `ttmt.stats`; a legacy number is reproduced with `ttmt.scoring` and
says so. [Open: a test that the two agree on isolated events, or a note in each module of where they differ; review
thread.]

## Model statistics

| Term | Definition | Computed by |
|:--|:--|:--|
| **G** | radiation operator of a node: (2 ns 3) x 6, P and S displacement per moment-tensor component, rows ordered phase, station, component, with 1/d | `task1x_mtglrt.radiation_operator`, `ttmt.radiation` |
| **Q** | orthonormal basis of the column space of G (rank <= 6) | `task1x_mtglrt.basis` |
| **Q_w** | basis of diag(1/sd) G, sd the noise level per station: the only correct basis when data are divided per station (ADR-0012) | ADR-0012 option; `WHITEN-REVIEW/run.py` (independent 3-line version) |
| **S_mt** | ½ ‖Qᵀ φ(v)‖², v the read vector of a hypothesis, φ the cap | `c10_lib`, `ttmt.detector` |
| **CAP6** | the registered cap φ: the norm of each (station, phase) 3-component read vector of the normalised matched-filter output limited to 6 | `c10_lib._vectors`, `CAPV = 6.0` |
| **f** | fraction of a vector's energy in the subspace: ‖Qᵀv‖² / ‖v‖² (noise-free unless said) | FE-FIX `run.py` |
| **f_emp** | f of the vector read on a real event at its catalogue position | FE-FIX step 3 |
| **margin m** | S_mt of a target at its true position over the 0.999 quantile of S_mt on the event-masked noise block at the same position | FE-FIX card, step 3 |
| **pair gain** (M6) | RSS of the best one-event fit minus RSS of the two-event joint fit of a near atom (whitened units) | `ttmt.v2.m6` (`gain = rss_one - rss_t`) |
| **G-M6a score** | min(s, pair gain) for a tested atom (owner, 15 November 2026) | `M6A-REALNOISE/run.py` (`m6-min`) |

## Data partitions

| Term | Definition |
|:--|:--|
| **tune partition** | the comparator tuning set b60d99b3d0bbbbfa: semisynthetic EGF sets `ss_ridgecrest_2019_egf_v1` and `ss_hinet_2012_egf_v1`, 50 records, n_T 223, n_C 305 (staged on the owner's laptop, `C:\Users\juanj\ttmt\tune`) |
| **CAL / EVAL** (R0) | Ridgecrest hour, origins in [0, 900) s / [900, 3560] s (`code/r0/r0_lib.py`) |
| **development / evaluation** (R1) | Hi-net 2012 swarm, origins before 840 s / the rest (`data/cards/hinet_swarm_2012.md`) |
| **reserved windows** | the confirmatory windows sealed in PR #50; never read |
| **spent** | a held-out set already read for a decision; see `docs/DATA_USE_LEDGER.md` |
