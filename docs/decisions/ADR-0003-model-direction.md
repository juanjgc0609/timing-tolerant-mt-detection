# ADR-0003: Direction of the next model

Date: 18 October 2026. Status: accepted as the working direction; each component enters only through its gate.

Context. The measured causes of the real-data failures are timing (0.35 to 0.41 s at Ridgecrest), readings on the
wrong energy peak (R1), saturation of the cap in continuous seismicity, and coda ghosts after deflation. Single-station
neural pickers are the strongest component of the comparators; a multi-station neural picker (PhaseNO, Sun et al.
2023) already improves on PhaseNet by using coherence across stations, which shows that the information PhaseNet
leaves unused is the joint structure of the network.

Decision. The next model combines (1) learned per-station evidence (picker probability fields, not thresholded picks),
(2) multi-station coupling constrained by physics (travel times and moment-tensor consistency), (3) an adaptive
bounded-influence statistic, (4) joint multi-event inference with deflation by empirical Green's functions, and (5)
calibrated false alarms. Full-wave modelling (elastodynamic solver, adjoints, PINN) is kept behind a gate: it is built
only if the failure atlas shows residual path errors in a frequency band where 3-D models predict waveforms.

Consequences. The comparator set includes PhaseNet, PhaseNO and EQTransformer with the best associators, template
matching and MALMI. The aim is to win in most regimes of the atlas, and the claim is tested, not assumed.
