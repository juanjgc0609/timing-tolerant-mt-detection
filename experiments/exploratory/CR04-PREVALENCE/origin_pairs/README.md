# CR-04, catalogue bound: events with another origin within 1 s (Ridgecrest 2019, SCSN)

Exploratory, 2026-11-14. Catalogue fields only; no waveform, pick or detection. This is not the registered CR-04
estimand (`docs/protocol/CR04_PREVALENCE.md` uses predicted arrivals on station subsets and the QTM catalogue). It is a
quick bound from the origin times of the network catalogue, used by the data plan of the real-data thread.

- Input: `s3://scedc-pds/earthquake_catalogs/SCEC_DC/2019.catalog`, SHA-256
  `8c8d3ec0a04f6da56c110879337aa5870366ecadc318d71bfcd67e178fa5f3b3` (fetched 2026-11-14).
- Selection: type `eq`, 2019-07-04 00:00 UTC + 30 days, 35.25-36.25 N, 118.1-117.1 W. Network CI (development).
- Command: `curl -sO https://scedc-pds.s3.amazonaws.com/earthquake_catalogs/SCEC_DC/2019.catalog && python -I
  ridgecrest_origin_pairs.py 2019.catalog`. Deterministic (no seed). Output: `ridgecrest_origin_pairs_out.txt`.

Results (measured): 28,060 events; 96 (0.34 %) have the next event within 1 s and 10 km; 1.89 % have another event
within 1 s at any distance in the box, against 3.7 % expected from a Poisson process at the local rate (events in a
centred 10-minute window).

Reading (hypothesis, not tested): the gap between 1.89 % and 3.7 % suggests that the catalogue misses roughly half of
the events with a neighbour within 1 s. Other explanations are not excluded: the Poisson reference ignores spatial
clustering and the magnitude-dependent dead time of the network, and the 10-minute rate smooths bursts. The 49 pairs
under 0.1 s are probably duplicates in large part (not checked). The registered measurement with QTM decides.
