# Licences and redistribution of the real datasets (plan T3.R7)

Rule of the project: the repository holds no waveform. Every dataset is fetched by a script from its archive and
checked against committed checksums or listings. Derived products (tables, figures, QC reports, catalogue aggregates)
are published with the citation of each source. Raw waveforms go to Zenodo only where the last column says yes.

"Verified" means read on the source's own page or file on the date given; "to verify" means stated from the
repository's earlier records or common practice and not yet read at the source.

| Dataset | Source and citation | Terms | Raw waveforms to Zenodo | Status |
|:--|:--|:--|:--|:--|
| Ridgecrest 2019, one hour (R0) | SCEDC, doi:10.7909/C3WD3xH1; network CI, Caltech/USGS; copy in `benchflow-ai/skillsbench` (Apache-2.0) | SCEDC data are open; cite SCEDC and the network | no: fetch script instead | to verify |
| Alto Tiberina 2010 (C16, 2M) | INGV network IV; copy in `avuan/PyMPA37` (GPL-3.0) | INGV data policy: open with attribution | no: fetch script instead | to verify |
| Hi-net swarm 2012 (R1) | NIED Hi-net; copy in `Geolab-USTC/MatchLocate2` (GPL-3.0) | NIED terms of use: registered users; no redistribution to third parties | no | to verify |
| `conf_geysers_2024` (reserved) | NCEDC, doi:10.7932/NCEDC; networks BG (LBNL, doi:10.7914/SN/BG) and NC (USGS, doi:10.7914/SN/NC) | open access through the NCEDC AWS open dataset; acknowledge NCEDC and the networks | no: listing and fetch script | to verify |
| `conf_antelope_valley_2021`, `conf_bodie_2020` (reserved) | NCEDC; networks NN (Nevada Seismological Laboratory, doi:10.7914/SN/NN), NC, BK (doi:10.7932/BDSN) | as above; NN metadata from EarthScope | no: listing and fetch script | to verify |
| `conf_colombia_cm` (slot) | Servicio Geologico Colombiano, network CM, doi:10.7914/SN/CM | to be set by the SGC reply (T3.R6) | only if the SGC allows it | pending |

The NCEDC acknowledgement to use in every publication: "Waveform data, metadata, or data products for this study
were accessed through the Northern California Earthquake Data Center (NCEDC), doi:10.7932/NCEDC."
