# CR-04: prevalence of the target regime in dense real sequences

Status: rule and code fixed 2026-11-01, measurement not run. Owner of the measurement: the reference-catalogue thread.
Sources: ADR-0008 items 1, 4 and 7; ATLAS-0 card, "Regime grid" and "Endpoint"; critique register CR-04.

## 1. Why this is needed

ADR-0008 item 4 weights each stratum of a regime by its share in real sequences divided by its share among the
evaluated targets, and item 7 words the claim "superior in the regime T, which holds Y % of the events of dense
sequences". CR-04 adds that if Y is below about 1 % the gain is of little practical use. No measured share exists yet.

## 2. What is already fixed, and stays fixed

- **ATLAS-0 and the phase-4 gate G4 use weights 1** (each evaluated target counts once). The card was anchored on
  2026-10-24 (`reg/ATLAS-0/registration-1`) before any prevalence table existed, and it states this default, as
  ADR-0008 item 4 requires. Replacing it would need an amendment tag (ADR-0004); none is proposed. Tuning (ADR-0009
  item 2) also uses weights 1.
- **The registrations R2 and R3 need the table.** Their `<table>` placeholders (prevalence weights, ADR-0008 item 4)
  are filled from Section 4 below, and the table's SHA-256 enters their freeze (SAP conformance item 17). If the
  measurement is not available when a registration freezes, that registration uses weights 1, says so, and the claim
  carries no value of Y.

## 3. Estimand

For a sequence with reviewed catalogue E (origins, hypocentres) and station inventory S:

1. Strata are the ATLAS-0 axes that a catalogue and an inventory determine: separation (5 bins: < 0.1, 0.1-0.5,
   0.5-1, 1-3, > 3 s), station subset k in {4, 6, 8, 12} and distance to the catalogue centroid (< 10, 10-30,
   > 30 km). The signal-to-noise ratio and K_mec/E need waveforms and mechanisms; they are not reweighted (factor 1
   within each stratum), and the table says so.
2. For each event and each k, the subset is `ttmt.atlas.cells.station_subset` (the k stations nearest the network
   centroid that are available over [origin - 10 s, origin + 60 s], `cells.availability`, ties by station code).
   Availability comes from the gap-free data intervals of the archive; where only station epochs are available, the
   epochs are used and the table is flagged.
3. Arrival times are predicted with `ttmt.traveltime` and the same velocity model the atlas uses for that network.
   The separation is `ttmt.atlas.cells.separation`: the median over the subset stations of |t_P(e) - nearest P or S
   arrival of any other catalogued event|, the definition of the card.
4. Y_k is the fraction of the catalogued events with separation in [0.1, 1) s (regime T). The weight table is the
   share of each stratum, per sequence and pooled with equal weight per sequence.
5. Interval: 95 % percentile interval of a block bootstrap over 1-hour blocks of origin time (2,000 replicates,
   seed fixed in the script and recorded with the table).

## 4. Catalogues

The densest reviewed catalogue of each atlas sequence, as named in its data card: Ridgecrest 2019, the QTM
catalogue of Ross et al. (2019); Hi-net swarm 2012, the Match and Locate catalogue (`catalog.dat`); Alto Tiberina
2010, the TABOO template catalogue. A sequence whose catalogue is not reviewed is reported but does not enter the
pooled table.

## 5. Known bias: the measured Y is a lower bound

Analyst and template-matching catalogues rarely list two events with origins less than about 1 s apart (CR-04): a
matched-filter detection suppresses others within its window, and analysts merge overlapping picks. The catalogue Y is
therefore reported as a **lower bound**. Two checks are reported with it, and neither replaces it in the weights:
- the minimum and the 1st percentile of inter-origin times in each catalogue, which show its dead time;
- a Poisson reference: with the local rate lambda(t) of the catalogue (events in a centred 10-minute window), the
  expected share of events with another event inside the T window, computed with the same predicted arrivals by
  drawing partner origins from lambda(t).
If the Poisson reference exceeds the catalogue Y by more than a factor of 2, the claim states both values.

## 6. Where it runs

The catalogues and inventories come from SCEDC, NIED and INGV. The cloud container cannot reach them (proxy 403 on
2026-11-01), so the measurement runs on the author's machine, only after the author asks for it. The functions are in
`ttmt.atlas.prevalence` (tested in `tests/test_prevalence.py`); the command is
`experiments/exploratory/CR04-PREVALENCE/run.py` (catalogue, stations and optional gap-free data intervals as CSV;
velocity model `hk` for Ridgecrest, `jma2001` for Hi-net; Alto Tiberina has no layered model in the package yet, so
it runs with `homogeneous:VP` and the table says so; seed 20261101). Each output table and its SHA-256 are committed
together, before the first of R2 or R3 freezes. Separations are computed over events within 60 s of origin, which
gives the same bins as all events for travel times below 57 s (`separation_windowed`).
