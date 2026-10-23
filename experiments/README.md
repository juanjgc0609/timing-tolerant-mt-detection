# Experiments

Layout of the research protocol (`docs/protocol/RESEARCH_PROTOCOL.md`, Sections 1 and 5).

- `exploratory/<id>/`: stage S0 and S1 work. Each folder holds a card (`card.md`, from
  `docs/protocol/templates/experiment_card.md`), the code it ran and its outputs. Nothing here is confirmatory.
- `registered/<id>/`: created at stage S2 for every registered test, with its registration, freeze list, run commands,
  result entry and verification, anchored under ADR-0004.

The registered tasks 2A to R1 predate this layout. Their code stays frozen in `code/` and their results in `results/`;
the error register and the sealed log are the index to them.

| ID | Stage | Question |
|:--|:--|:--|
| `exploratory/P6-power` | S1 (done) | How many targets does a paired recall test need? |
| `exploratory/ATLAS-0` | S0 | Which regimes, budgets and counting rule define the failure atlas, fixed before it runs? |
| `exploratory/THEORY-5` | S1 (done) | Do the statements of Appendix AB hold numerically? |
