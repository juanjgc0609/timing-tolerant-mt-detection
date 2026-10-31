"""Reproduction gate G3a, target 2: PhaseNet in-domain on ETHZ with the SeisBench pick-benchmark code
(Münchmeyer et al. 2022, JGR), against the published metrics. Runs on the GPU machine.

    git clone https://github.com/seisbench/pick-benchmark %BENCH% && git -C %BENCH% checkout d039fd8  (jgr_revision)
    python repro_ethz.py prepare --out %W%\\weights\\ethz_phasenet     (published ETHZ weights as a checkpoint)
    (download the paper's evaluation targets to %W%\\targets\\ethz, link in the benchmark README)
    cd %BENCH%\\benchmark
    python eval.py %W%\\weights\\ethz_phasenet %W%\\targets\\ethz --sets dev,test
    python collect_results.py %W%\\pred %W%\\results.csv
    python repro_ethz.py report %W%\\results.csv published_ethz.json --out %W%\\ethz.json  (refused until verified)
    python gate.py write --puente ... --ethz %W%\\ethz.json

`prepare` wraps ``seisbench.models.PhaseNet.from_pretrained("ethz")`` (the weights the benchmark trained and exported)
in the folder layout eval.py reads (version_0 with hparams.yaml, metrics.csv and one checkpoint), so the benchmark's own
evaluation and metric code run unchanged; nothing is retrained. `report` compares the test metrics of results.csv with
published_ethz.json, {"source": "<paper, table or figure>", "metrics": {"test_det_auc": v, "test_phase_mcc": v,
"test_P_rmse_s": v, ...}}, whose values are copied from the paper with their source; it writes the gate's evidence
format ({"conditions": [{"name", "value", "published"}], "source"}).
"""
import argparse
import json
import sys
from pathlib import Path

CONFIG = {"model": "PhaseNet", "data": "ETHZ", "model_args": {"sample_boundaries": [None, None], "lr": 1e-3}}


def prepare(out, weights="ethz"):
    import seisbench
    import seisbench.models as sbm
    import torch
    import yaml
    net = sbm.PhaseNet.from_pretrained(weights)
    ver = Path(out) / "version_0"
    (ver / "checkpoints").mkdir(parents=True, exist_ok=True)
    hp = dict(lr=1e-3, sigma=20, sample_boundaries=[None, None], phases=net.labels, norm=net.norm,
              in_channels=net.in_channels, classes=net.classes)
    try:
        import lightning
        plv = lightning.__version__
    except ImportError:
        import pytorch_lightning
        plv = pytorch_lightning.__version__
    torch.save({"state_dict": {"model." + k: v for k, v in net.state_dict().items()}, "hyper_parameters": hp,
                "pytorch-lightning_version": plv, "epoch": 0, "global_step": 1},
               ver / "checkpoints" / "epoch=0-step=1.ckpt")
    (ver / "hparams.yaml").write_text(yaml.safe_dump(CONFIG))
    (ver / "metrics.csv").write_text("epoch,step,val_loss\n0,0,0.0\n")
    meta = dict(weights=weights, seisbench=seisbench.__version__, phases=net.labels, norm=net.norm)
    (Path(out) / "provenance.json").write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta))


def report(results, published, out):
    import pandas as pd
    res = pd.read_csv(results)
    row = res[(res["data"].str.lower() == "ethz") & (res["model"].str.lower() == "phasenet")]
    if len(row) != 1:
        raise SystemExit(f"{results}: expected one ETHZ PhaseNet row, found {len(row)}")
    pub = json.load(open(published))
    if not pub.get("verified"):
        raise SystemExit(f"{published}: published values not verified against the paper (see verification_note)")
    cond = [{"name": k, "value": float(row.iloc[0][k]), "published": float(v)} for k, v in pub["metrics"].items()]
    ev = {"conditions": cond, "source": f"pick-benchmark on ETHZ vs {pub['source']}"}
    Path(out).write_text(json.dumps(ev, indent=1))
    for c in cond:
        print(f"{c['name']:20s} {c['value']:.4f}  published {c['published']:.4f}  "
              f"rel {(c['value'] - c['published']) / c['published']:+.3f}")


def main(argv=None):
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--out", required=True)
    p.add_argument("--weights", default="ethz")
    r = sub.add_parser("report")
    r.add_argument("results")
    r.add_argument("published")
    r.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "prepare":
        prepare(a.out, a.weights)
    else:
        report(a.results, a.published, a.out)


if __name__ == "__main__":
    sys.exit(main())
