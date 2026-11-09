# Reserved confirmatory sequences (plan T3.R4)

Sealed on 9 November 2026. Nobody reads these sequences before the phase-6 registration that uses them: no waveform
sample, pick, detection or per-event catalogue entry. The reservation exists so that phase 6 has data that no
development, tuning or atlas step has seen.

## Files

| File | What it is |
|:--|:--|
| `reservation.json` | The sealed list: selection rule, criteria, sequences (window in UTC, disc of radius `radius_km`), exclusions. `sealed_sha256` is the canonical SHA-256 of `sequences`. Written by `tools/data/seal_reservation.py`, which `tests/test_registry.py` re-runs with `--check`. |
| `<id>.listing.tsv` | Archive listing of every continuous-waveform file of the sequence inside the disc and the window: S3 key, size, ETag (content hash). Pins the bytes without downloading them. Written by `tools/data/seal_listing.py`. |

## How the seal is enforced

`ttmt.registry.check_read(network, lat, lon, t0, t1, purpose)` raises `SealedError` for any read that touches a
sealed window, whatever dataset or path it comes from; every station of every network in the disc is sealed, and a
station of unknown position counts as inside. The purposes `fetch` (store and hash bytes) and `qc-metadata` (headers,
StationXML, catalogue aggregates; `ttmt.qc` without `samples`) pass. A confirmatory run opens a sequence only with
`TTMT_UNSEAL=<id>:<registration id>`, where the registration id is listed under `unseal` for that sequence by an
amendment commit made after the registration is anchored, and appears in `preregistration/PREREGISTRATION_LOG_EN.txt`.
Every unsealed read is logged to `<work root>/unseal_log.jsonl`. Changing a sequence after sealing changes the hash,
and the registry then refuses the whole file.

## Sequences (9 November 2026)

| Id | Role | Networks | Window (UTC) | Catalogue aggregates | Three-component stations | Open item |
|:--|:--|:--|:--|:--|:--|:--|
| `conf_geysers_2024` | confirmatory | BG, NC | 2024-03-06 to 2024-04-05 | 2,348 events; Mc 1.0 (maximum curvature), 1.2 (goodness of fit 90 %); 624 at M >= 1.0; median horizontal error 0.33 km | 41 of 46 within 30 km on 2024-03-15, all with StationXML | 182 GB in the disc of 60 km (BG at 500 Hz); the run needs only the evaluated subsets |
| `conf_antelope_valley_2021` | confirmatory | NN, NC, BK | 2021-07-08 22:00 to 2021-08-07 22:00 | 2,213 events; Mc 1.0 / 2.2 (moving Mc 0.6 to 2.8: short-term incompleteness after the M 6.0); 1,405 at M >= 1.0 | NN stations only near the source; their StationXML is not in the bucket | NN metadata from EarthScope, then QC; if fewer than 8 stations pass, `conf_bodie_2020` replaces it |
| `conf_bodie_2020` | spare | NN, NC | 2020-04-11 14:00 to 2020-05-11 14:00 | 745 events; Mc 1.6 / 1.3; 477 at Mc | as above (NN) | used only if Antelope Valley fails QC |
| `conf_colombia_cm` | confirmatory slot | CM (SGC) | chosen by the rule in the file | - | - | data access requested from SGC (T3.R6); the window is added by amendment and sealed from that commit |

The aggregates come from `tools/data/scan_catalogue.py --days 30` on the NCEDC catalogues of 2019 to 2026, and the
station counts from `tools/data/s3_availability.py`; both read catalogue fields, archive listings and StationXML only.

## Why these and not others

- **Networks.** CI, IV and N. (Hi-net) are development networks; BG, NN and CM were never used. Using NC (in
  Antelope Valley) for a confirmatory sequence means NC must not become development data either.
- **Training leakage (CR-22).** `phasenet/original` and `phaseno` were trained on NCEDC recordings of 1984-2019
  (`src/ttmt/associate/provenance.py`). Every reserved NCEDC window starts after 2019: no temporal overlap, a
  regional overlap that the phase-6 report declares. The Long Valley swarm of September 2019 meets every other
  criterion and was excluded for this reason; it can be development data.
- **Truth for the target regime (CR-04).** The NCSN catalogue lists almost no pair of events less than 1 s apart
  (1 in 2,348 at the Geysers). The Geysers network is dense enough for the decimation design: truth from the full
  network, evaluation on 4 to 12 station subsets that exclude the truth stations. The truth protocol belongs to the
  reference-catalogue thread.
- **Size.** Phase 6 needs 400 to 1,200 targets per sequence; both confirmatory sequences have 600 or more events at
  or above Mc in 30 days.

## External anchor

The hash of `reservation.json` and of each listing is in the commit that adds them. The author pushes the tag
`reserve/T3.R4-1` on that commit (the cloud cannot push tags):

```bash
git tag -a reserve/T3.R4-1 -m "T3.R4 reservation, sealed_sha256 <hash>" <commit> && git push origin reserve/T3.R4-1
```
