DRAFT, not registered. Built from `docs/protocol/templates/registration.md`. Fields in <angle brackets> are filled
from the data card of the reserved sequences (T3.R4) and the frozen code before anchoring; nothing here may be
changed after any reserved waveform is read.

Task R2 (__FECHA__, wall clock), written before any method reads the evaluation partition. Its SHA-256 and the
commit that carries it are anchored by tag and Release `reg/R2/registration-1`, OpenTimestamps and OSF Registries
(ADR-0004) before the run.
Origin. Phase-4 gate G4 met on development under the rule of ADR-0008 (report <path>); statistical analysis plan
`docs/protocol/PHASE6_SAP.md` v<version>, whose Sections 3 to 10 are part of this registration.
Question. On real sequences from networks never used in development, does the model detect more catalogued events
than every tuned comparator in the target regime (arrival separation 0.1 to 1 s) and no fewer, by more than 0.05 of
recall, elsewhere, at 30 false detections per hour?
Data. <sequences, networks, periods, file list with SHA-256 from the sealed data card>. Development partition: records
ending at least 1 h before the evaluation period (causal split, SAP 3.2). Evaluation partition: <exact UTC
boundaries>, margins included. Eligibility by the frozen script <path, SHA-256> (SAP 3.3); excluded blocks listed.
Targets and truth. Dense or nodal network <name> with its catalogue <name, version>; evaluation on station subsets of
4 to 12 stations of the same network that exclude the stations used to build the truth (network decimation, CR-04),
subsets drawn by <rule> before the run. Targets: catalogue events with M >= Mc(t) and at least 4 usable stations of the
subset. Excluded events and neutral detections as in SAP 3.3; enhanced catalogue for the sensitivity of SAP 3.4:
<name>. Regime T can be tested only if the truth resolves pairs at 0.1 to 1 s with origin-time standard deviation sigma_t
(measured on T3.V) at most a ninth of their separation (SAP 3.4): sigma_t = <value>, pairs qualifying <n>. Expected targets in T
and C from the catalogue: <numbers>.
Methods. Model <tag, version>; comparators: core set of ATLAS-0 (PhaseNet and PhaseNO with GaMMA and PyOcto, EQNet,
template matching where the truth is not its own, BOMP-PK and MT-FW) and the best comparator with the model's picker
if fine-tuned (CR-03), each tuned with 60 configurations under ADR-0008 item 9 (ledger <path, SHA-256>); shared
waveforms, station list, velocity model and typed station terms; pretrained-weight provenance table <path> (CR-22).
Endpoint. Recall under LOC+T-V (ADR-0008 item 2) at 30 false detections per hour, exact FROC, operating point per
method on the evaluated record, prevalence weights <table, rule in docs/protocol/CR04_PREVALENCE.md> (ADR-0008 item 4); duplicator negative control on the best
comparator (item 3); deployable threshold fixed on development (secondary).
Rules. H1 (only if T is testable) superiority in T, then H2 non-inferiority by 0.05 in C, each an intersection-union
of pairwise tests over the comparators (H1 block sign-flip at one-sided 0.025; H2 clustered score test at the level
<a2> of the design check `sample_size.py calibrate --blocks <B> --p0 <range> --m <range> --k <range> --rho <range> --icc <range>`, ends
from `sample_size.py estimate` on development <path>, ADR-0008 item 7), fixed
sequence; blocks of 1 h (or T/30) within each
sequence; secondary Holm family and descriptive readings as in SAP Section 2; intervals SAP 7.5; budget check SAP 7.6.
Simulated with `sample_size.py` for this design (rho <r>, icc <i>, <m> targets per block, f <f>, +0.10 in T): power of
the chain <p>, size <s> (output <path>).
Prediction. <declared from development>.
Interpretation. H1 and H2 met: the model is superior in T, which holds <Y> % of the events of dense sequences, and not
inferior by more than 0.05 elsewhere, on this real data. H2 met alone: not inferior on real data; the T claim rests on
R3. H2 not met: the model is worse than at least one tuned comparator outside T on real data, reported with the
discordant targets. The duplicator control positive: the analysis gives no beat claim.
Freeze and run. Code, environment locks, weights, ledgers and the eligibility script frozen (`reg/R2/freeze-1`); single
run; non-binding futility look after the first sequence (SAP Section 8); verification by reimplementation in a separate
AI session plus ten discordant targets re-derived by a person.
Declared limits. <network, region, period, magnitude range; truth from one network; decimation is not an independent
network>.
