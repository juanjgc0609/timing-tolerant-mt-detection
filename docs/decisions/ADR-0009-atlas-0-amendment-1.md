# ADR-0009: ATLAS-0 amendment 1 (minimum event interval of REAL and PyOcto, tuning objective, `instance` overlap)

Date: 24 October 2026. Status: proposed by the atlas thread after the questions of the comparator thread (PR #18),
before any tuning trial is drawn. It takes effect when anchored as `reg/ATLAS-0/amendment-1` (ADR-0004 item 1). The
seeds stay those of `reg/ATLAS-0/registration-1` (`experiments/exploratory/ATLAS-0/seeds.json`); the amendment does not
re-derive them.

## Context

The ATLAS-0 card was anchored on 24 October 2026 (tag `reg/ATLAS-0/registration-1` on 7bafb7a). Preparing the first
tuning run, the comparator thread raised three points:
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

## Decision

1. REAL: `tint` log-uniform in [0.1, 5] s (replaces "window in [1, 5] s"). Its other parameters are unchanged.
   PyOcto: `min_interevent_time` log-uniform in [0.1, 5] s is added to its search space; "tolerance" is the velocity
   model's tolerance and "association cutoff" is `pick_match_tolerance`, both with the ranges of the card. Any other
   method with a fixed minimum time between output events is listed in the ledger with its value.
2. Tuning objective, for every method, the prototype and the phase-4 model included: the unweighted mean of the
   development recall in T and the development recall in C, each at 30 false detections per hour under LOC+T-V and
   weighted within its regime as in the endpoint (ADR-0008 item 4). The ledger stores both recalls and the mean; the
   configuration with the largest mean is selected, ties broken by the lower trial number.
3. `instance` overlap: the ledger and the atlas report declare that PhaseNet `instance` was trained on data that
   include the Alto Tiberina 2010 period. On Alto Tiberina records, results of pipelines using `instance` are reported
   separately and flagged "training-period overlap"; the comparator chosen there is the better of `original` and
   `instance` as tuned, so the overlap can only make the comparator stronger. No claim of the atlas rests on an
   `instance` result on Alto Tiberina alone. Other weights are checked the same way when the ledger lists them.

## Consequences

REAL and PyOcto can separate events inside regime T, so their failures there are measured, not built in. Tuning is reproducible from
the ledger with one stated objective. The `instance` overlap is visible and conservative for any beat claim, and the
atlas does not attribute a picker success on Alto Tiberina to generalisation. The rule of ADR-0008 and the comparator
set are unchanged; PyOcto, a core comparator, becomes stronger in T, which makes H1 harder, not easier.
