DRAFT, not registered. Built from `docs/protocol/templates/registration.md`. Fields in <angle brackets> are filled
before anchoring; nothing here may be changed after any reserved waveform is read.

Task R3 (__FECHA__, wall clock), written before any method reads the evaluation partition. Anchored by tag and Release
`reg/R3/registration-1`, OpenTimestamps and OSF Registries (ADR-0004) before the run.
Origin. Phase-4 gate G4 met on development under ADR-0008 (report <path>); analysis plan `docs/protocol/PHASE6_SAP.md`
v<version>, Sections 3 to 10 part of this registration.
Question. When recorded events of networks never used in development are injected into those networks' own
continuous noise at controlled separations, does the model detect more events than every tuned comparator at arrival
separations of 0.1 to 1 s and no fewer, by more than 0.05 of recall, elsewhere, at 30 false detections per hour?
Data. Reserved networks <list> (T3.R4). Noise library: <hours> of continuous record with no catalogued event above
<M> within <km>; event library: <n> recorded events, disjoint by event and by time block from every library a method
may use for templates, deflation or priors (with the near-duplicate check of CR-06: no injected event with waveform
correlation above <c> to a library event). Development: records ending at least 1 h before the evaluation period.
Targets and truth. Injected events; design grid (separation in 5 bins, signal-to-noise ratio, station subsets 4, 6,
8, 12, distance to the calibration centroid, K_mec/E recorded as covariate); one block = one network realisation (one
noise block, one station subset, its injected events). Seeds derived from the anchor (ADR-0004 item 4): first 8 bytes
of SHA-256(<tag commit hash> || <Release created_at>). Number of blocks: <B>, from `sample_size.py solve` with rho
<r>, icc <i>, f <f>, +0.05 in T (output <path>).
Methods. As in R2 (core set, model's picker arm, 60-configuration budget, provenance table).
Endpoint. Recall under LOC+T-V at 30 false detections per hour, operating point per method on the evaluated record,
weights = real prevalence of each cell / its share in the design (<table, rule in docs/protocol/CR04_PREVALENCE.md>, ADR-0008 item 4); duplicator control.
Rules. H1 superiority in T, then H2 non-inferiority by 0.05 in C, intersection-union of pairwise tests (H1 block
sign-flip at one-sided 0.025; H2 clustered score test at the level <a2> of the design check `sample_size.py calibrate` over the range
from `sample_size.py estimate` on development <path>,
ADR-0008 item 7), fixed sequence; secondary Holm family (SAP Section 2); intervals stratified by network (SAP 7.5);
budget check stratified by parent event and the e-value check at 3 to 15 km and 1 to 4 s (SAP 7.6, CR-01). Power of the
chain <p>, size <s>.
Prediction. <declared from development>.
Interpretation. H1 and H2 met: superior in T and not inferior elsewhere under known truth on unseen networks. H1 not
met: no demonstrated gain where the theory predicts one; the cells are described. The duplicator control positive: no
claim. The e-value rate above 2p on parent-event hypotheses: the M7 guarantee does not cover event-induced false
detections, and that is reported.
Freeze and run. As in R2; no interim look.
Declared limits. Injection adds signals linearly to real noise; no wave interaction; events of one network injected
into its own noise.
