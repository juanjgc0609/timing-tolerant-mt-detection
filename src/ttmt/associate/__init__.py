"""Comparators of the failure atlas and their tuning (thread "Comparadores y ajuste"; card
experiments/exploratory/ATLAS-0/card.md, sections "Tuning" and "Methods and comparators").

Modules:
  space       search spaces of the card and the shared random search (60 trials, one seed list for every method)
  ledger      the append-only tuning ledger (experiments/exploratory/ATLAS-0/tuning_ledger.csv), freeze, sensitivity
  tune        the objective (development recall in T and in C at 30 false detections per hour, LOC+T-V) and the loop
  velocity    the one 1-D model, in the form each program takes
  base        the frame of a run and the conversions to and from ttmt.atlas.schema
  gamma_assoc, pyocto_assoc, real_assoc   associator wrappers (GaMMA always with ncpu = 1)
  pickers     probability traces once per weights, picks per trial
  pipeline    picker + associator pipelines and the status of every method of the card
  provenance  training data of the pretrained weights against the atlas data (CR-22)
The programs themselves (GaMMA, PyOcto, REAL, SeisBench) are imported only inside the wrappers; the rest needs numpy.
"""
