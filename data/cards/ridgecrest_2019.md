# Data card: Ridgecrest 2019, one hour

- Source: repository `benchflow-ai/skillsbench`, commit `9a1f4dd5f7659f75707435da3ce854b6e48321d1`, folder
  `tasks/earthquake-phase-association/`. Original data: SCSN (Caltech/USGS); cite SCEDC.
- Record: 2019-07-04 19:00-20:00 UTC, CI network; SCSN analyst catalogue of 102 events, M 1.49 to 4.50.
- Checksums: `data/checksums/ridgecrest_2019_R0.sha256`.
- Partitions used by R0 (`code/r0/r0_lib.py`, `T_CAL`, `T_EVAL`): calibration = origins in [0, 900) s, evaluation =
  origins in [900, 3560] s, with t = 0 at 19:00:00 UTC.
- Uses: R0 (registered, evaluation partition read once); explorations C18a to C18j on the calibration partition.
  Error notes: E50 (truncated search in C18a), E53 (station terms in-sample in C18), E54 (+SC terms of R0).
- Status: development only from now on. Any later use is exploratory.
