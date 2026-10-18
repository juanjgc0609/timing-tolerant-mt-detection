# Data card: Hi-net swarm 2012

- Source: repository `Geolab-USTC/MatchLocate2`, commit `bac5ccb00816a4ec0aaab215ea1c342c730d412f`, folder `example/`
  (GPL-3.0; the waveforms are NIED Hi-net data under the NIED terms of use).
- Record: Fukushima-Yamagata swarm of 2 September 2012, 33 min 20 s, seven Hi-net stations, three components; 14 JMA
  events (`catalog.dat`). The target list of R1 (54 events of magnitude 0.5 or more) is the Match and Locate output,
  so it cannot be used to evaluate template matching.
- Checksums: `data/checksums/hinet_matchlocate2_example.sha256`.
- Partitions used by R1: development = origins before 03:34:00 (840 s, 30 targets); evaluation = the rest of the
  record. The development reads extended to 880 s (E56); the rule now is that development reads end at the start of
  the evaluation segment, margin included.
- Uses: R1 (registered, run once). Post hoc: timing readings (`code/r1/posthoc_timing_cal.py`). Open item: two targets
  of magnitude 2 or more that no method detects, to be reviewed by hand as possible artefacts of the truth list.
- Status: development only from now on.
