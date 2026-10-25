# ADR-0009: ATLAS-0 amendment 1 (minimum event interval of REAL and PyOcto, tuning objective, `instance` overlap)

Date: 24 October 2026. Status: proposed by the atlas thread after the questions of the comparator thread (PR #18),
before any tuning trial is drawn. It takes effect when anchored as `reg/ATLAS-0/amendment-1` (ADR-0004 item 1). The
seeds stay those of `reg/ATLAS-0/registration-1` (`experiments/exploratory/ATLAS-0/seeds.json`); the amendment does not
re-derive them.

## Context

The ATLAS-0 card was anchored on 24 October 2026 (tag `reg/ATLAS-0/registration-1` on 7bafb7a). Preparing the first
tuning run, the comparator thread raised three points, and the generator thread a fourth:
1. The REAL search space gives "window in [1, 5] s". The parameter it maps to in REAL is `tint`, the minimum time
   between two events. With `tint` >= 1 s REAL cannot output two events closer than 1 s, so in regime T (arrival
   separation in [0.1, 1) s) it is a straw man by construction, and prediction P-RE1 (separation 0.5-3 s) would be
   forced rather than measured. PyOcto, in the core set, has the same limit: `min_interevent_time` is not in its
   search space and keeps its default of 3 s, so it cannot output two events closer than 3 s. The card's "association
   cutoff in [0.5, 2] s" is mapped to `pick_match_tolerance` (the only PyOcto association cutoff in seconds;
   `association_cutoff_distance` is in km), and "tolerance in [0.3, 2] s" to the velocity model's tolerance.
2. The card says the tuning objective is "development recall in T and in C at 30 false detections per hour under
   LOC+T-V, weighted as in the endpoint", but it does not say how the two recalls combine into one number.
3. CR-22 requires the ledger to declare overlaps of pretrained weights with an atlas network and period. The SeisBench
   PhaseNet weights `instance` are trained on INSTANCE (Italy, 2005 to 2020), which covers the Alto Tiberina 2010
   noise records.
4. The generator thread found that the 30/70 split "by a hash of the hour" cannot split the Ridgecrest set (one
   hour) or the Hi-net set (2,000 s), and asked for seeds for the shared event library, which no record owns. Its
   realisations were seeded from a hash of a commit string, which is known before the Release exists and so is
   weaker than ADR-0004 item 4.

## Decision

1. REAL: `tint` log-uniform in [0.1, 5] s (replaces "window in [1, 5] s"). Its other parameters are unchanged.
   PyOcto: `min_interevent_time` log-uniform in [0.1, 5] s is added to its search space; "tolerance" is the velocity
   model's tolerance and "association cutoff" is `pick_match_tolerance`, both with the ranges of the card. Any other
   method with a fixed minimum time between output events is listed in the ledger with its value.
2. Tuning objective, for every method, the prototype and the phase-4 model included: the unweighted mean of the
   development recall in T and the development recall in C (C holds every target outside T, the isolated events
   with separation above 3 s included, as in `ttmt.atlas.froc.gate` and ADR-0008), each at 30 false detections per
   hour under LOC+T-V and weighted within its regime as in the endpoint (ADR-0008 item 4). The ledger stores both
   recalls and the mean; the configuration with the largest mean is selected, ties broken by the lower trial number.
   If the tuning partition has no target in T or in C, the objective is undefined (no fallback to the other recall)
   and tuning stops until the partition is fixed.
3. `instance` overlap: the ledger and the atlas report declare that PhaseNet `instance` was trained on data that
   include the Alto Tiberina 2010 period. On Alto Tiberina records, results of pipelines using `instance` are reported
   separately and flagged "training-period overlap"; the comparator chosen there is the better of `original` and
   `instance` as tuned, so the overlap can only make the comparator stronger. No claim of the atlas rests on an
   `instance` result on Alto Tiberina alone. Other weights are checked the same way when the ledger lists them.
4. Partition unit: the unit is the noise block of 300 s, not the hour, because a set of one hour
   (Ridgecrest) or of 2,000 s (Hi-net) cannot be split by the hour. Blocks are the consecutive 300 s windows of a
   record, start = `rec.t0` + 300 i in seconds on the record's time axis (`rec.t0` = 0 for these sets), the last
   partial block dropped (`ttmt.synth.realdata.blocks`). Only noise-role blocks enter the split; EGF-role blocks
   (`realdata.block_roles`) are in neither partition. A noise block goes to tuning when
   `int.from_bytes(sha256(f"{dataset}|{start:.3f}|tune-read-v1".encode()).digest()[:8], "little") / 2**64 < 0.3`
   (`realdata.tune_read`). Per dataset, if no noise block falls in tuning, the one with the smallest value moves to
   tuning; if none falls in reading, the one with the largest value moves to reading. Blocks do not overlap and
   every realisation and injected event follows its block.
5. Seeds of semi-synthetic data, one rule, all from the anchored master of registration 1: realisation r on noise
   block i of set s draws its events from `ttmt.atlas.seeds.realisation_rng(s, i, r)` (stream "inject", key
   (first 4 bytes big-endian of sha256(s), i, r)) and its contamination stations from
   `realisation_rng(s, i, r, "subset")`; event library j (drawn whole: medium, event sequence, receiver effects;
   Ridgecrest j = 0, Hi-net j = 1) from `library_rng(j)` (stream "inject", key 1000000 + j). Seeds hashed from a
   commit or any other string are not used.

## Consequences

REAL and PyOcto can separate events inside regime T, so their failures there are measured, not built in. The
registered lines of the card are not edited; its Amendment 1 section carries the changes. Tuning is reproducible from
the ledger with one stated objective. The `instance` overlap is visible and conservative for any beat claim, and the
atlas does not attribute a picker success on Alto Tiberina to generalisation. The rule of ADR-0008 and the comparator
set are unchanged; PyOcto, a core comparator, becomes stronger in T, which makes H1 harder, not easier.
