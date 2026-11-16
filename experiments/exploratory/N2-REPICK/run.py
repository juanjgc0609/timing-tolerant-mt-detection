"""N2-REPICK: "pick, subtract, pick again" in the ideal case (card.md, declared 15 November 2026 before running).

Question: of the T targets that frozen PhaseNet+PyOcto loses at the picker stage on the tuning partition
(b60d99b3d0bbbbfa), how many get at least min_picks of their own arrivals picked once the waveform of their pair
partner is subtracted exactly (the injected waveform itself, not an estimate)?

Two stages, so that the part without waveforms can be checked anywhere:

  baseline  needs only the tuning manifests (truth + stored PhaseNet probabilities). Reproduces the frozen objective,
            lists the T targets lost at the picker (STAGE-LOSS rule) and computes the control K: the same stored picks
            matched with the partner's arrivals removed from the truth (no subtraction).
  repick    needs the set files (sig, noise) under $TTMT_WORK_ROOT/synth/<set id>/ and the real record under
            $TTMT_DATA_ROOT (to rebuild every injected event on its own, with the anchored seed, and check that the
            events add up to the stored signal). Subtracts the partner, runs PhaseNet again with the frozen weights,
            threshold and peak distance, and matches as in K.

    python run.py baseline <tune_phasenet.json>... --out lost.json
    python run.py repick lost.json [--limit N] --out repick.json

Tuning partition only; nothing here reads a read-partition file, a reserved window or a ledger.
"""
import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

REPO = Path(__file__).resolve().parents[3]
FROZEN = REPO / "experiments/exploratory/COMPARATORS/frozen/phasenet+pyocto.json"
spec = importlib.util.spec_from_file_location("sl", REPO / "experiments/exploratory/STAGE-LOSS/stage_loss.py")
SL = importlib.util.module_from_spec(spec)
spec.loader.exec_module(SL)
TU = SL.TU
EGF_KW = {"ridgecrest_2019": {}, "hinet_2012": {"iso_pre": 4.0, "iso_post": 15.0, "iso_dm": 1.0, "pre": 4.0,
                                               "post": 20.0}}         # tools/build_semisynth.py
SIG_RTOL = 1e-3          # rebuilt events must add up to the stored signal within this share of its peak


def without(truth, drop):
    """the two fields SL.match reads, with the arrivals of event `drop` removed."""
    a = np.array(truth.arr_t_obs, float)
    a[drop] = np.nan
    return SimpleNamespace(arr_t_obs=a, n_events=truth.n_events)


def partner(truth, e):
    p = int(truth.pair[e])
    if p < 0:
        return None
    o = [int(j) for j in np.flatnonzero(truth.pair == p) if j != e]
    return o[0] if len(o) == 1 else None


def baseline(a):
    from ttmt.associate.pipeline import associator
    fz = json.load(open(FROZEN))
    cfg = dict(fz["params"])
    recs = [r for m in a.manifests for r in TU.load_manifest(m, fz["method"])]
    part = TU.partition_hash(recs)
    if part != fz["partition"]:
        sys.exit(f"partition {part} is not the frozen one {fz['partition']}")
    acfg = {k: v for k, v in cfg.items() if k.startswith("as_")}
    mp = int(cfg["as_min_picks"])
    cache, dets, pk = {}, [], {}
    for r in recs:
        lab = r.name.split("__")[0]
        if lab not in cache:
            p = SL.phasenet_picks(r, cfg)
            cache[lab] = associator("pyocto").associate(r.data["frame"], p, acfg)[0], p, SL.match(p, r.truth)
        d, p, m = cache[lab]
        dets.append(d)
        pk[r.name], pk[r.name + "#picks"] = m, p
    ev = SL.evaluate(recs, dets)
    st = SL.attribute(ev, recs, dets, pk, mp)
    print(f"partition {part}  T {ev['recall_T']:.4f} (frozen {fz['recall_T']})  C {ev['recall_C']:.4f} "
          f"(frozen {fz['recall_C']})")
    if abs(ev["recall_T"] - fz["recall_T"]) > 5e-4 or abs(ev["recall_C"] - fz["recall_C"]) > 5e-4:
        sys.exit("the frozen objective is not reproduced: stop")
    lost = []
    for j, (r, e, sep, _snr, _tot) in enumerate(ev["meta"]):
        if not ev["T"][j] or st[j] != "picker":
            continue
        f = partner(r.truth, e)
        n0 = int(pk[r.name][1][e])
        nk = int(SL.match(pk[r.name + "#picks"], without(r.truth, f))[1][e]) if f is not None else n0
        lost.append({"record": r.name, "label": r.name.split("__")[0], "set_id": r.name.split("_0")[0],
                     "event": int(e), "partner": f, "sep_s": float(sep), "n_picked": n0, "n_control": nk})
    nT = int(ev["T"].sum())
    K = sum(x["n_control"] >= mp for x in lost)
    print(f"T targets {nT}; lost at the picker {len(lost)}; without a partner in the truth {sum(x['partner'] is None for x in lost)}")
    print(f"control K (stored picks, partner's arrivals removed from the truth): {K} of {len(lost)} reach {mp} picks")
    json.dump({"partition": part, "min_picks": mp, "pk": {k: cfg[k] for k in cfg if k.startswith("pk_")},
               "manifests": [str(m) for m in a.manifests], "n_T": nT, "lost": lost}, open(a.out, "w"), indent=1)


def rebuild(set_id, cache={}):
    """callable label -> (data, [waveform of each injected event], truth, codes, fs) for one set."""
    if set_id in cache:
        return cache[set_id]
    from ttmt.synth import inject as INJ, realdata as RD, sets
    from ttmt.atlas import seeds as SEEDS
    rd = sets.SetReader(set_id, "tune", verify="fingerprint")
    man = rd.man
    rec = RD.load_record(man["dataset"])
    atoms, atom_ids = RD.template_candidates(rec)
    src = [INJ.Source.from_egf(g, rec.fs) for g in RD.egf_library(rec, **EGF_KW[man["dataset"]])]
    real = {r["label"]: r for r in man["realisations"]}

    def one(label):
        r = real[label]
        comps, orig = [], INJ.radiated_vectors

        def spy(wave, *args, **kw):                       # realise calls it once per target, in target order
            comps.append(np.array(wave, np.float32))
            return orig(wave, *args, **kw)
        INJ.radiated_vectors = spy
        try:
            sig, _n, truth = INJ.realise(rec, tuple(r["block"]), src, SEEDS.realisation_rng(set_id, r["block_index"], r["r"]),
                                         man["inject_config"], atoms, atom_ids, label)
        finally:
            INJ.radiated_vectors = orig
        s0, n0 = rd.arrays(label)
        t0 = rd.truth(label)
        err = float(np.abs(sig - s0).max() / max(np.abs(s0).max(), 1e-30))
        tot = float(np.abs(np.sum(comps, axis=0) - s0).max() / max(np.abs(s0).max(), 1e-30))
        if len(comps) != int(t0.target.sum()) or err > SIG_RTOL or tot > SIG_RTOL or not np.allclose(truth.t, t0.t):
            raise RuntimeError(f"{label}: the rebuilt realisation is not the stored one (sig {err:.2e}, sum {tot:.2e}, "
                               f"{len(comps)} events for {int(t0.target.sum())} targets)")
        return s0 + n0, comps, t0, list(man["stations"]), float(man["fs"])
    cache[set_id] = one
    return one


def repick(a):
    from ttmt.associate import pickers
    from ttmt.associate.records import stream_of
    L = json.load(open(a.lost))
    mp, pkc = L["min_picks"], L["pk"]
    lost = [x for x in L["lost"] if x["partner"] is not None]
    labels = sorted({x["label"] for x in lost})
    if a.limit:
        labels = labels[:a.limit]
    done, t00 = {}, time.time()

    def picks_of(data, codes, fs, t0):
        prob, dt, tp = pickers.annotate_seisbench(stream_of(data, codes, fs, t0), codes, "PhaseNet", pkc["pk_weights"],
                                                  t0=t0)
        return pickers.picks_from_probabilities(prob, dt, pkc["pk_threshold"], pkc["pk_peak_distance"], tp)

    for li, lab in enumerate(labels):
        set_id = next(x["set_id"] for x in lost if x["label"] == lab)
        data, comps, truth, codes, fs = rebuild(set_id)(lab)
        if li == 0:                                       # the picker of this machine against the stored one
            p = picks_of(data, codes, fs, truth.t0)
            n_new = SL.match(p, truth)[1]
            ref = next(x for x in lost if x["label"] == lab)
            print(f"check on {lab}: {len(p.t)} picks on the unmodified record; event {ref['event']} has "
                  f"{int(n_new[ref['event']])} picked arrivals here, {ref['n_picked']} with the stored probabilities",
                  flush=True)
            if abs(int(n_new[ref["event"]]) - ref["n_picked"]) > 1:
                sys.exit("the picker of this machine does not reproduce the stored picks: stop")
        for f in sorted({x["partner"] for x in lost if x["label"] == lab}):
            p = picks_of(data - comps[f].astype(float), codes, fs, truth.t0)
            n = SL.match(p, without(truth, f))[1]
            for x in lost:
                if x["label"] == lab and x["partner"] == f:
                    done[(lab, x["event"])] = int(n[x["event"]])
        print(f"{lab}: done ({li + 1}/{len(labels)}, {time.time() - t00:.0f} s)", flush=True)
    rows = [dict(x, n_repick=done[(x["label"], x["event"])]) for x in lost if (x["label"], x["event"]) in done]
    n = len(rows)
    K = sum(x["n_control"] >= mp for x in rows)
    S = sum(x["n_repick"] >= mp for x in rows)
    both = sum(x["n_control"] >= mp and x["n_repick"] >= mp for x in rows)
    uniq = {(x["label"], x["event"]): x for x in rows}
    Su = sum(x["n_repick"] >= mp for x in uniq.values())
    Ku = sum(x["n_control"] >= mp for x in uniq.values())
    print(f"\ntargets lost at the picker and re-picked: {n} (of {len(L['lost'])} lost; {len(L['lost']) - len(lost)} have no partner in the truth)")
    print(f"reach {mp} picks: control K {K}, after subtraction S {S}, both {both}; net S - K = {S - K}")
    print(f"unique targets (a target counts once whatever k): {len(uniq)}; K {Ku}, S {Su}, net {Su - Ku}")
    for lo, hi in ((0.1, 0.25), (0.25, 0.5), (0.5, 1.01)):
        rr = [x for x in rows if lo <= x["sep_s"] < hi]
        print(f"  separation {lo:.2f}-{hi:.2f} s: n {len(rr)}, K {sum(x['n_control'] >= mp for x in rr)}, "
              f"S {sum(x['n_repick'] >= mp for x in rr)}")
    d = np.array([x["n_repick"] - x["n_control"] for x in rows])
    print(f"picked arrivals per target, after minus control: median {np.median(d):+.0f}, mean {d.mean():+.2f}, "
          f"share with more {np.mean(d > 0):.2f}, with fewer {np.mean(d < 0):.2f}")
    labs = sorted({x["label"] for x in rows})              # interval: bootstrap over realisations (labels), seed fixed
    per = np.array([sum((x["n_repick"] >= mp) - (x["n_control"] >= mp) for x in rows if x["label"] == l) for l in labs])
    rng = np.random.default_rng(20261120)
    bs = per[rng.integers(0, len(per), (10000, len(per)))].sum(1)
    print(f"net S - K, 95 % bootstrap interval over {len(labs)} realisations: [{np.quantile(bs, .025):.0f}, "
          f"{np.quantile(bs, .975):.0f}]")
    if not a.limit:
        net = S - K
        zone = "PASSES the target" if net >= 40 else ("MIDDLE zone: target not met, effect reported" if net >= 15
                                                       else "LOW zone")
        print(f"\nREADING (card.md): net {net} of {len(L['lost'])} -> {zone}")
    json.dump({"rows": rows, "K": K, "S": S, "limit": a.limit}, open(a.out, "w"), indent=1)


def snr(a):
    """Step 2 of card.md: is a target that stays below min_picks after the subtraction a weak one? Per channel, the
    peak of the target's own injected waveform over the robust sd of the noise of the record."""
    L, R = json.load(open(a.lost)), json.load(open(a.repick))["rows"]
    mp, rows = L["min_picks"], []
    for lab in sorted({x["label"] for x in R}):
        set_id = next(x["set_id"] for x in R if x["label"] == lab)
        data, comps, truth, codes, fs = rebuild(set_id)(lab)
        noise = (data - np.sum(comps, axis=0)).reshape(-1, data.shape[-1])
        sd = 1.4826 * np.median(np.abs(noise - np.median(noise, axis=1, keepdims=True)), axis=1)
        for x in R:
            if x["label"] == lab:
                r = np.abs(comps[x["event"]].reshape(-1, data.shape[-1])).max(axis=1) / np.maximum(sd, 1e-30)
                rows.append(dict(x, snr_med=float(np.median(r)), n_ch5=int((r >= 5).sum()), n_ch=int(r.size)))
    rec = np.array([x["n_repick"] >= mp for x in rows])
    print(f"targets {len(rows)}; recovered after subtraction {int(rec.sum())}, not recovered {int((~rec).sum())}; "
          f"array shape of the last record {data.shape}, {len(codes)} stations")
    for k in ("snr_med", "n_ch5"):
        v = np.array([x[k] for x in rows], float)
        t, f = v[rec], v[~rec]
        auc = float(np.mean(t[:, None] > f[None, :]) + 0.5 * np.mean(t[:, None] == f[None, :]))
        labs = sorted({x["label"] for x in rows})
        rng, bs = np.random.default_rng(20261121), []
        for _ in range(2000):                              # bootstrap over realisations
            pick = rng.choice(labs, len(labs))
            vv = np.concatenate([[x[k] for x in rows if x["label"] == l] for l in pick]).astype(float)
            rr = np.concatenate([[x["n_repick"] >= mp for x in rows if x["label"] == l] for l in pick])
            if rr.any() and (~rr).any():
                tt, ff = vv[rr], vv[~rr]
                bs.append(np.mean(tt[:, None] > ff[None, :]) + 0.5 * np.mean(tt[:, None] == ff[None, :]))
        print(f"{k}: recovered median {np.median(t):.2f} (q25 {np.quantile(t, .25):.2f}, q75 {np.quantile(t, .75):.2f}); "
              f"not recovered median {np.median(f):.2f} (q25 {np.quantile(f, .25):.2f}, q75 {np.quantile(f, .75):.2f}); "
              f"AUC {auc:.3f}, 95 % interval [{np.quantile(bs, .025):.3f}, {np.quantile(bs, .975):.3f}]")
    v = np.array([x["snr_med"] for x in rows])
    for s_ in ("ridgecrest", "hinet"):
        m = np.array([s_ in x["set_id"] for x in rows])
        print(f"  {s_}: n {int(m.sum())}, recovered {int(rec[m].sum())}, median snr_med {np.median(v[m]):.2f}")
    v = np.array([x["snr_med"] for x in rows], float)
    auc = float(np.mean(v[rec][:, None] > v[~rec][None, :]) + 0.5 * np.mean(v[rec][:, None] == v[~rec][None, :]))
    zone = "SUPPORTED" if auc >= 0.75 else ("PARTLY supported" if auc >= 0.60 else "NOT supported")
    print(f"\nREADING (card.md, step 2): AUC of snr_med {auc:.3f} -> low-signal hypothesis {zone}")
    json.dump({"rows": rows}, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="stage", required=True)
    b = sub.add_parser("baseline")
    b.add_argument("manifests", nargs="+")
    b.add_argument("--out", default="lost.json")
    r = sub.add_parser("repick")
    r.add_argument("lost")
    r.add_argument("--limit", type=int, default=0)
    r.add_argument("--out", default="repick.json")
    q = sub.add_parser("snr")
    q.add_argument("lost")
    q.add_argument("repick")
    q.add_argument("--out", default="snr.json")
    a = ap.parse_args()
    {"baseline": baseline, "repick": repick, "snr": snr}[a.stage](a)
