# Real seismic data downloadable in this environment (search of 9 October 2026)

From the environment, github.com (git clone), raw.githubusercontent.com and GitHub release assets work. The FDSN servers (IRIS, SGC, GEOFON), zenodo, huggingface and storage.googleapis.com respond 403 due to the proxy policy, so only repositories that version the waveforms in git are usable.

## Candidate 1 (recommended): Ridgecrest 2019, one hour of aftershocks, 14 stations

Source: benchflow-ai/skillsbench, path tasks/earthquake-phase-association/ (Apache-2.0 licence); original data from the CI network (SCSN/Caltech) and analyst-reviewed SCSN catalogue. Local copy: /home/claude/realdata/skillsbench_eqassoc/ (wave.mseed, stations.csv, catalog.csv; 28.8 MB). One continuous hour, 2019-07-04 19:00-20:00 UTC, between 1.5 and 2.5 h after the Mw 6.4. Fourteen CI stations at 16-48 km from the sequence centroid; eleven broadband HHE/HHN/HHZ (CLC, SRT, CCC, TOW2, LRL, WRC2, WBM, SLA, JRC2, MPM, WCS2) and three with only EHZ plus accelerometer (WNM, WVP2, WRV2). 100 Hz, raw counts, no gaps, 45 traces. Inter-station distance 3.8-71.6 km (median 32.6 km). stations.csv gives coordinates and scalar sensitivity per channel, without poles and zeros. Catalogue: 102 events (M 1.49-4.50; 23 with M >= 2.5), with hypocentre, without picks. Strong overlap: minimum inter-event time 1.5 s, median 22.9 s, 17 pairs less than 10 s apart and 64 less than 30 s apart. Limitations: catalogue incomplete below M 1.5-2 (precision must be evaluated with that bias, for example by restricting to M >= 2), no reference picks or mechanisms, scalar-only response. Complement: AI4EPS/PhaseNet demo/mseed (CCC and CLC, 7 h, no catalogue).

## Candidate 2: Fukushima-Yamagata swarm, 2012-09-02, Hi-net, 7 stations

Source: Geolab-USTC/MatchLocate2, folder example/ (GPL-3.0; NIED Hi-net data, whose redistribution is subject to NIED). Local copy: /home/claude/realdata/Geolab-USTC_MatchLocate2/example/ (continuous data in Trace/20120902/). 33 min 20 s continuous, 7 three-component stations at 100 Hz, at 12-30 km from the swarm, coordinates in SAC headers. Catalogues: 14 JMA events (M 2.1-3.2) and 140 Match&Locate detections (relative magnitude -0.41 to 3.2). Extreme overlap (minimum 6.3 s, median 11.9 s). Limitations: prefiltered data (energy in 1-15 Hz), no instrument response, 1 Hz borehole sensors with known horizontal orientation errors (relevant for a moment-tensor detector), ground truth derived from another detector.

## Candidate 3: Alto Tiberina, Italy, 2010-07-23, 24 h, 14 stations

Source: avuan/PyMPA37 (GPL-3.0; INGV IV network). Local copy: /home/claude/realdata/pympa37_altotiberina_100723/ (142 MB). 24 h, three components decimated to 25 Hz, stations at 1.5-39 km from the cluster; TABOO catalogue with 32 events that day (M -0.8 to 1.9). Little overlap; useful mainly as a false-alarm control over 24 hours of real noise.

## Unusable or limited

FAST (Hector Mine; 6 of 7 stations vertical only, no catalogue), EQcorrscan (90 s event windows, not continuous), MALMI (5 minutes at Hengill, no catalogue), LOKI, GPHS445, fast_matched_filter, easyQuake, seismo_live, EQTransformer (isolated stations or event windows). LOC-FLOW, QuakeFlow, GaMMA, REAL, QuakeMigrate, BackTrackBB, PyOcto, SeisBench, PAL/PALM/PAD/MFT, Seismic_BPMF and HinetPy do not version continuous waveforms and download from blocked servers.

## Recommendation

Ridgecrest is the best case for a first real test: raw, 100 Hz, broadband, analyst catalogue and real overlap; a sparse subset of 4 to 8 HH* stations can be extracted. The Hi-net swarm serves as a second case of extreme overlap with 7 stations; Alto Tiberina as a false-alarm control.
