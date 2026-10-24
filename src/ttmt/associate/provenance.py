"""Training provenance of the pretrained weights the comparators use, against the atlas data (critique item CR-22;
card, "Tuning": "the tuning ledger lists the training sets, regions and periods of every pretrained model, and an
overlap with an atlas network and period is declared").

Each entry states its source; "verify" marks a period taken from the dataset paper's abstract or a secondary page
rather than from the dataset files, which the real-data thread (owner of CR-22 together with this one) checks against
the dataset metadata. :func:`overlaps` gives, for every weight and atlas dataset, "period" (same region and
overlapping period: the evaluation data may be in the training set), "region" (same region, other period) or "".
"""
WEIGHTS = {
    "phasenet/original": dict(
        dataset="NCEDC catalogue recordings (Zhu and Beroza, 2019, GJI 216, 261-273, Section 2: 779,514 recordings)",
        region="northern California", period=(1984, 2019), source="Zhu and Beroza (2019); arXiv 1803.03211",
        verify="period: 'over thirty years of recordings' (abstract); exact years not stated there"),
    "phasenet/190703-214543": dict(
        dataset="same model as phasenet/original, TensorFlow checkpoint of AI4EPS/PhaseNet at 62005c6 (R0, R1)",
        region="northern California", period=(1984, 2019), source="AI4EPS/PhaseNet", verify="as original"),
    "phasenet/instance": dict(dataset="INSTANCE (Michelini et al., 2021, ESSD 13, 5509)", region="Italy",
                              period=(2005, 2020), source="SeisBench weight docs; Michelini et al. (2021)",
                              verify="period from the dataset paper"),
    "phasenet/stead": dict(dataset="STEAD (Mousavi et al., 2019, IEEE Access)", region="global",
                           period=(1984, 2018), source="SeisBench weight metadata (trained on STEAD, 100 epochs)",
                           verify="period from the dataset paper; STEAD includes Italian and Californian stations"),
    "phasenet/scedc": dict(dataset="SCEDC (SeisBench benchmark dataset)", region="southern California",
                           period=(2000, 2020), source="SeisBench benchmark datasets page ('from 2000-2020')",
                           verify=""),
    "phasenet/ethz": dict(dataset="ETHZ (SeisBench benchmark dataset)", region="Switzerland", period=None,
                          source="SeisBench benchmark datasets page", verify="period not stated on the page"),
    "eqtransformer/original": dict(dataset="STEAD (Mousavi et al., 2020, Nat. Commun.)", region="global",
                                   period=(1984, 2018), source="Mousavi et al. (2020)", verify="period"),
    "eqtransformer/instance": dict(dataset="INSTANCE", region="Italy", period=(2005, 2020),
                                   source="SeisBench weight docs", verify="period"),
    "phaseno/epoch=19-step=1140000": dict(
        dataset="NCEDC events 1984-2019 with STEAD noise (Sun et al., 2023, GRL; arXiv 2305.03269)",
        region="northern California", period=(1984, 2019), source="Sun et al. (2023), training data section",
        verify=""),
    "eqnet/none": dict(dataset="no released weights for the 2022 EQNet model (environments/DEVIATIONS.md)",
                       region="", period=None, source="AI4EPS/EQNet at 7dd24ed", verify=""),
}

# atlas data (ttmt.config.DATASETS keys): region and year of the records
ATLAS_DATA = {
    "ridgecrest_2019": dict(region="southern California", year=2019, network="CI (SCSN)"),
    "hinet_2012": dict(region="Japan", year=2012, network="Hi-net"),
    "altotiberina_c16": dict(region="Italy", year=2010, network="Alto Tiberina (TABOO) and INGV"),
    "altotiberina_days": dict(region="Italy", year=2010, network="Alto Tiberina (TABOO) and INGV"),
}

# regions contained in a weight's training region
_CONTAINS = {"global": {"southern California", "northern California", "Italy", "Japan", "Switzerland"},
             "California": {"southern California", "northern California"}}


def _same_region(train, data):
    return train == data or data in _CONTAINS.get(train, set())


def overlaps():
    """{weight: {dataset: "period" | "region" | ""}}."""
    out = {}
    for w, e in WEIGHTS.items():
        row = {}
        for d, a in ATLAS_DATA.items():
            if not e["region"] or not _same_region(e["region"], a["region"]):
                row[d] = ""
            elif e["period"] is None or e["period"][0] <= a["year"] <= e["period"][1]:
                row[d] = "period"
            else:
                row[d] = "region"
        out[w] = row
    return out


def table_markdown():
    ds = list(ATLAS_DATA)
    lines = ["| Weights | Training data | Region | Period | " + " | ".join(ds) + " |",
             "|:--|:--|:--|:--|" + ":--|" * len(ds)]
    ov = overlaps()
    for w, e in WEIGHTS.items():
        per = "-".join(map(str, e["period"])) if e["period"] else "not stated"
        lines.append(f"| {w} | {e['dataset']} | {e['region'] or '-'} | {per} | "
                     + " | ".join(ov[w][d] or "-" for d in ds) + " |")
    return "\n".join(lines)
