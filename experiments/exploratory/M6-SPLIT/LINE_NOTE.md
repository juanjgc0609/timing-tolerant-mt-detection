# M6-SPLIT line: ledger row and thesis note (ready to paste, 22 November 2026)

Nothing was run for this file. `docs/DATA_USE_LEDGER.md` is not on main yet, and `docs/thesis/Thesis.md` belongs to
the thesis thread. Both texts below are ready for those owners to paste as they are.

## Row for `docs/DATA_USE_LEDGER.md`, section 2 (attempts per line of work)

| M6 S-wave splitting (constraint inside the joint fit) | 4: M6-SPLIT (window-only constraint), M6-SPLIT2 (shared with 10 isolated events; implementation defect found post hoc), M6-SPLIT3 (estimated once from isolated events, then fixed; run by M6-SPLIT2's declared condition), M6-SPLIT2b (second attempt of M6-SPLIT2 with converged isolated fits) | negative on validity in all four: 24, 126, 57 and 114 of 200 control scenes split (limit 12); G-M6a not run with any of them | `M6-SPLIT/card.md`, `M6-SPLIT2/card.md`, `M6-SPLIT3/card.md`, `M6-SPLIT2b/card.md` |

Simulated scenes only, new seeds for each attempt (20261121-23, 20261124-27, 20261131-34, 20261141-45); no real
or reserved data was read.

## Note for the thesis (results on model v2, M6; negative result)

S-wave splitting enters the M6 resolvability test as a constraint inside the joint fit, as the director asked on
20 November 2026: only S is split, with a fast axis and a delay of at most 0.2 s per station. In simulation on the
M6-CAL scenes (8 stations, mismatched pulses, SNR 12), four declared attempts fail the validity criterion. That
criterion allows at most 12 of 200 single-source control scenes to be split, and it was not relaxed.

| attempt | how the splitting is estimated | control scenes split | near-pair recall (reference without splitting) | AUC |
|:--|:--|--:|--:|--:|
| M6-SPLIT | inside the fit, pair window only | 24 | 0.744 (0.894) | 0.885 |
| M6-SPLIT2 | shared with 10 isolated events per station | 126 | 0.850 (0.912) | 0.931 |
| M6-SPLIT3 | once from the isolated events, then fixed | 57 | 0.806 (0.906) | 0.873 |
| M6-SPLIT2b | as M6-SPLIT2, isolated fits converged | 114 | 0.850 (0.931) | 0.904 |

Without splitting in the data, the same test splits 1 to 6 of 200 controls. With splitting in the data and a model
that ignores it, the test splits 34 to 48 of 200 controls. The largest spurious gain in a control is normally negative
(median about -100). With the splitting shared with isolated events, it becomes positive (median +115, and +87 after
converging the isolated fits), and the cause is not diagnosed. Splitting shared with isolated events was not
supported by these results. A second atom can fit the S-wave residual that a splitting estimated from other events
leaves behind. That explanation is a hypothesis and has not been tested. M6-SPLIT2 counts as a failed attempt even
though its defect was found afterwards. M6-SPLIT2b is reported as a second attempt, not as a replacement.
