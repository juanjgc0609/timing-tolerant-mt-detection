# Data access

## What the working environment can reach (checked 18 October 2026)

| Host | Reachable | Use |
|:--|:--|:--|
| github.com through git, raw content, release assets (objects.githubusercontent.com), LFS (media.githubusercontent.com) | yes | code, example datasets hosted in repositories, PhaseNet weights |
| pypi.org | yes | packages (ObsPy, SeisBench, PyOcto, EQcorrscan, QuakeMigrate) |
| scedc-pds.s3.amazonaws.com, s3.us-west-2.amazonaws.com | no | SCEDC continuous waveforms and catalogues of Southern California on AWS |
| ncedc-pds.s3.amazonaws.com | no | NCEDC continuous waveforms of Northern California on AWS |
| service.iris.edu | no | FDSN web services (waveforms, stations, events) of many networks |
| earthquake.usgs.gov | no | ComCat catalogues |
| huggingface.co and its LFS hosts | no | labelled datasets of AI4EPS (QuakeFlow) |
| hifis-storage.desy.de | no | SeisBench benchmark datasets and pretrained weights (EQTransformer, PhaseNet variants) |
| zenodo.org | no | published datasets and model weights |

Refused hosts can be allowed by the account owner in the claude.ai settings (Capabilities, allowed domains).

## Datasets held

| Dataset | Use so far | Status |
|:--|:--|:--|
| Ridgecrest 2019, one hour, SCSN (via benchflow-ai/skillsbench) | R0 development and evaluation | development data from now on |
| Alto Tiberina 2010, noise days (via avuan/PyMPA37) | C16, 2M | development data from now on |
| Hi-net swarm 2012, seven stations (via Geolab-USTC/MatchLocate2) | R1 | development data from now on |

No unused real dataset remains.

## What the next phases need

1. **Failure atlas (phase 3):** several sequences with dense seismicity and reviewed catalogues, from at least two
   networks, with continuous three-component data; labelled picks for picker evaluation.
2. **Confirmatory tests (phase 6):** at least two sequences never used for development, from networks different from
   the development ones, each with an independent reviewed catalogue.
3. **Training or fine-tuning of learned components (phase 4):** a labelled dataset of the scale used by PhaseNet or
   PhaseNO, or pretrained weights; this also needs compute beyond the two cores of the working environment.

## Options

(a) Allow the domains above in the account settings (recommended: scedc-pds.s3.amazonaws.com,
s3.us-west-2.amazonaws.com, ncedc-pds.s3.amazonaws.com, service.iris.edu, earthquake.usgs.gov, huggingface.co,
hifis-storage.desy.de, zenodo.org). (b) Download on a personal computer and attach the files or link the computer to
the session. (c) Use only GitHub-hosted examples, which limits the confirmatory phase to small datasets.

Every dataset enters through a data card (`data/cards/<dataset>.md`) with its source, version or commit, checksums,
catalogue, partition boundaries and the record of its uses.
