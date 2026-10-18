# ADR-0007: Primary hit rule for future tests

Date: 18 October 2026. Status: adopted as the default of the review of the plan
(`docs/reviews/2026-10-18_plan_review.md`, Section 2.2); the author can reverse it.

Context. The thesis describes the LOC rule as "epicentre within 5 km and origin within 3 s", but the scoring code
applies no test on the origin time. In Task 2M, 40 per cent of the LOC hits have origin errors above 0.5 s. Duplicate
detections of an event that is already hit are neutral in `froc_exact.py`. The rescoring of the registered tasks is
carried out separately and will be logged in the error register with its effect on each result.

Decision. For every test registered from now on:
1. The primary hit rule is **LOC+T**: epicentre within 5 km and origin time within 0.5 s.
2. Matching is one to one: detections are taken in score order, each one is matched to the nearest unmatched event
   within tolerance, and a detection with no unmatched event left counts as a false detection.
3. Results under the epicentre-only rule are reported as secondary, and duplicates are reported in a separate column
   (protocol, Section 4.4).

Consequences. Registered results keep their registered rule, and their rescoring under this rule is reported as post
hoc. The failure atlas (phase 3) and the confirmatory tests (phase 6) use this rule for every method.
