# ADR-0002: Package, experiments and archive

Date: 18 October 2026. Status: accepted, to be implemented in phase 1.

Context. About 270 scripts in a flat layout, absolute paths, module-level state changed at run time, and the
registered results depending on the exact bytes of that code.

Decision. (1) A Python package `src/ttmt/` holds every function that a confirmatory test uses, with configuration
objects instead of module globals. (2) The current code moves unchanged to `archive/legacy/`, with its hash index, so
that every registered result stays reproducible from the bytes that produced it. (3) Experiments live in
`experiments/registered/` and `experiments/exploratory/`, each with a card, a configuration and its commands.
(4) A regression suite proves that the package reproduces the registered outputs of 2F, 2M, R0 and R1 before it is
used for anything new. (5) Known defects (E51, E52, E54, E56) are fixed as explicit, versioned conventions.

Consequences. No registered number changes. Corrected conventions produce new, separately labelled results.
