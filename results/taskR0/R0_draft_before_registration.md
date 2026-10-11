# Task R0 (draft, not registered): pilot with public real data, Ridgecrest 2019

Status: DRAFT. No detector has been run on these data; it was only verified that the traces can be read (inventory of stations, channels, sampling and window). This draft becomes a registration (with the hash sent before running) only with the author's approval.

## Why

All the evidence for G1, G1-R and G1-C is synthetic and from a single generator family, and the four audits name the first result with real data as the main gap. The registered test G1' needs RSNC data that are not accessible from this environment. One hour of Ridgecrest aftershocks with an analyst catalogue is accessible (repository benchflow-ai/skillsbench, Apache-2.0 licence, CI network data and SCSN catalogue). R0 does not replace G1': it is a preliminary pilot, in another region and with another network, that tells whether the advantage measured on synthetic data survives on real data with overlapping events.

## Data

One continuous hour, 2019-07-04 19:00-20:00 UTC, 1.5 to 2.5 h after the Mw 6.4. Eleven three-component broadband CI stations (CLC, SRT, CCC, TOW2, LRL, WRC2, WBM, SLA, JRC2, MPM, WCS2) at 16-43 km from the centroid, 100 Hz, raw counts; scalar sensitivity per channel. SCSN catalogue of 102 events (M 1.49-4.50) with hypocentre; no picks. Minimum inter-event time 1.5 s, median 22.9 s.

## Proposed design

Temporal partition fixed before running: calibration 19:00-19:15 (S2 station terms by pick-free joint relocation with the catalogued events of that segment; noise level; thresholds), evaluation 19:15-20:00. Subnetworks fixed by rule before running: the 8 HH stations closest to the centroid, a subnetwork of 6, one of 4 and a one-sided one (the stations on a single side of the cluster in azimuth). Preprocessing: counts to velocity with the scalar sensitivity, 1-6 Hz bandpass, resampling to 16 Hz, 40 s windows with origins 5.0-8.0 s and a 3 s step, as in the benchmark. Medium: the project's homogeneous velocities plus station terms for the pick-free methods; the standard regional 1-D model of southern California for the pick-based methods.

Methods: MT-FW with cap 6, S2 terms and T3 refinement (the thesis method); MT-FW with cap 6 and catalogue corrections; BOMP-PK; HUBER-N; published comparators on the same data: PhaseNet (original model from the repository) with GaMMA and with PyOcto, and QuakeMigrate.

Ground truth and endpoint: catalogue events with M >= 2.0 (approximate completeness in the first hours of the sequence). Hit: epicentre within <= 5 km and origin within <= 0.5 s (LOC+T rule) and, as secondary, epicentre only (LOC). A detection without a match is classified into three classes, as Section 16 requires: it coincides with a catalogued event of M < 2.0, verified by correlation with templates of catalogued events, or unverified; only those verified as false count as false in the primary endpoint, and the conservative count (all non-matching ones as false) is also reported. Primary endpoint: recall of M >= 2.0 at an equal number of false detections per hour, read on each method's curve, with 60 s block bootstrap intervals. Criterion, in the form of G1': the thesis method loses no more than 0.02 of recall against the best published comparator and emits fewer false detections at equal recall in at least half of the subnetworks; any result is reported.

## Limits declared in advance

One hour, one sequence, one network; incomplete catalogue in the early sequence, which biases downward the apparent precision of all methods; scalar-only instrument response; region and network different from the RSNC. A positive result is evidence of transfer to real data, not the confirmation of G1'; a negative one delimits where it does not transfer.
