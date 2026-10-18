# ADR-0003: Direction of the next model

Date: 18 October 2026. Status: accepted as the working direction; each component enters only through its gate.
Revised on 18 October 2026 (review of the plan, `docs/reviews/2026-10-18_plan_review.md`): evidence on PhaseNO and
the closest prior work.

Context. The measured causes of the real-data failures are timing (0.35 to 0.41 s at Ridgecrest), readings on the
wrong energy peak (R1), saturation of the cap in continuous seismicity, and coda ghosts after deflation. Single-station
neural pickers are the strongest component of the comparators; a multi-station neural picker (PhaseNO, Sun et al.
2023) already improves on PhaseNet by using coherence across stations, which shows that the information PhaseNet
leaves unused is the joint structure of the network. That evidence is mixed outside Ridgecrest: on four local
hydraulic-fracturing networks (Kong et al., arXiv 2510.15281, GJI) PhaseNO had higher recall and better low
signal-to-noise picking, but PhaseNet kept the higher P-pick F1 in all four tests and the higher association F1 in one,
and PhaseNO degraded with few stations; overlapping events were not analysed. The 0.35 to 0.41 s timing error was
measured with a homogeneous medium and truncated readings (E50); the layered medium of C18i and C18j removes much of
it. Travel-time coupling of learned features across stations is already published (EQNet, Zhu et al. 2022), so the
novelty to be tested is narrower: moment-tensor consistency of amplitudes and polarisation, joint decomposition of
overlapping events, and a bounded-influence statistic.

Decision. The next model combines (1) learned per-station evidence (picker probability fields, not thresholded picks),
(2) multi-station coupling constrained by physics (travel times and moment-tensor consistency), (3) an adaptive
bounded-influence statistic, (4) joint multi-event inference with deflation by empirical Green's functions, and (5)
calibrated false alarms. Full-wave modelling (elastodynamic solver, adjoints, PINN) is kept behind a gate: it is built
only if the failure atlas shows residual path errors in a frequency band where 3-D models predict waveforms.

Consequences. The comparator set includes PhaseNet and PhaseNO with the best associators, EQNet, template matching,
and, as an extended set, EQTransformer and MALMI. PhaseNet runs with its peak distance and threshold tuned on
development, not the defaults (0.5 s, 0.3) used in R0 and R1, since the overlapping-event regime depends on them. The aim is to win in most regimes of the atlas, and the claim is tested, not assumed.
