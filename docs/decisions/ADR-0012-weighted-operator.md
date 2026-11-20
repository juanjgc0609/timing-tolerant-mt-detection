# ADR-0012: Weight the operator like the data (Q of diag(1/sd) G)

Date: 20 November 2026. Status: decided by the owner on 18 November 2026 (card "Arreglar pesos y correr el piloto",
8:56 PM Bogotá on 17 November); written by the FE-FIX thread from the supervisor's design
(`tras-reunion-director-2026-11-20.md`, section 1, in the project files). Number reserved by the master plan v2.18.

## Context

The model divides the data of each station by its noise level, but projects them on the span of an operator G
built from geometry alone. In `ttmt.r1.Net` the 16 Hz traces are divided by `station_scale` (src/ttmt/r1.py,
`self.d16 = d16 / self.scale`) while `node_operators(XN, STA)` builds Gs and Qst without weights. A noise-free event
v = G m becomes diag(1/sd) G m, which is not in span(G) when the sd differ between stations.

Evidence [measured, FE-FIX step 2, PR #68, laptop run 24b2e4a]: noise-free analytic events in the real geometries
of the tuning partition give f = ||Q^T v||^2 / ||v||^2 of 0.705 (Hi-net) and 0.771 (Ridgecrest) with Q from G, and
0.999 in both with Q from diag(1/sd) G. The max/min ratio of sd between stations is 5.70 and 7.63. No other element
of the chain (resampling, band-pass, Ricker, read window, peak rule, common shift) moves f by more than 0.02.

## Decision

1. Wherever the data of each station are divided by its own noise level, the operator is weighted the same way:
   G_w = diag(1/sd) G per node (rows ordered phase, station, component) and Q its orthonormal basis. Gs is weighted
   too, so the moment-tensor estimate and the deflation use the same operator as the projection.
2. Implementation: `ttmt.radiation.node_operators(XN, STA, weights=None)` and `station_weights`; `ttmt.r1.Net(...,
   weighted_q=False)` and `python -m ttmt.r1 p1 ... --weighted-q`. Every new path that divides by a per-station
   level must pass `weights = 1 / sd`.
3. **Default off.** The registered R1 path is unchanged: with `weights=None` the code is the old one, and the frozen
   lists reproduce bit for bit. The weighted option is the one used by model v2, pilot P1 and G-M6a.
4. **The band does not change** (owner, 18 November 2026): 0.5-6 Hz, 16 Hz, Ricker of 2 Hz. The question stays open
   for the director's seismologist adviser, with the numbers of FE-FIX step 3 (Hi-net +112 % [-9 %, +261 %],
   Ridgecrest -17 % [-48 %, +9 %] in the detection margin; the two sets point in opposite directions).

## Tests

`tests/test_adr0012.py`: with station sd in a ratio of 10 and no noise, f ≥ 0.99 for every event with the weighted
operator and a median f below 0.95 without it; `weights=None` gives Gs and Qst bit for bit equal to the per-node
`radiation_operator` and `basis` of the registered engines; unit weights give the same span; wrong-shaped, zero or
non-finite weights are rejected. The existing regression and legacy-equivalence tests cover the frozen outputs.

## Consequences

- Results computed with real noise and per-station scaling are affected (supervisor's review,
  `revision-blanqueo-2026-11-17.md`): prototype R0 and R1, C18a, G-M1a and the `--local` step of M5. Task 2M and the
  default M5 steps divide by one sd for the whole window and are not affected. They are re-run descriptively with
  the fix, on the same partitions and seeds, without changing any criterion already read; a registered negative
  stays registered. Nothing reserved or spent is read again for this.
- M5 `--local` (experiments/exploratory/M5-EVALUE/run_c.py) normalises per slot window; it is not changed here. Its
  re-run must build the weighted operator from the same per-window levels.
- Error register: E66.
