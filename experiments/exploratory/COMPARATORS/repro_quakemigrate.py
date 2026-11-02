"""Reproduction gate G3a, target 3: QuakeMigrate on its published example (Winder et al. 2026, Seismica, Section 4.1:
Rutford Ice Stream icequakes, 48 hours), against its stated event counts. Runs on the device.

    python repro_quakemigrate.py fetch   --work %W%     (station file and published catalogues from the manuscript
                                                         repository at the pinned commit; waveforms from EarthScope)
    python repro_quakemigrate.py lut     --work %W%
    python repro_quakemigrate.py detect  --work %W% [--threads 12]
    python repro_quakemigrate.py trigger --work %W%
    python repro_quakemigrate.py locate  --work %W% [--threads 4]
    python repro_quakemigrate.py report  --work %W% published_quakemigrate.json --out %W%\\qm.json
    python gate.py write --puente ... --ethz ... --quakemigrate %W%\\qm.json

Of the paper's examples only this one states event counts (the Askja 24-hour example gives none; the Askja manual-pick
comparison is a location benchmark). Every parameter is copied from the manuscript's scripts
(figures/icequake_example/generate_results/0_get_data.py to 5_collate_filter_catalogue.py at the commit in
published_quakemigrate.json) and is not tuned: homogeneous LUT (vp 3.841, vs 1.970 km/s, 50 m nodes, decimated 2x for
detect), STA/LTA onsets at 250 sps (P 20-124 Hz on *1, S 10-124 Hz on *[2,3]; windows P 0.01/0.25 s, S 0.05/0.5 s),
detect timestep 30 s, static trigger threshold 3.0 on the normalised coalescence, marginal window 0.06 s, minimum
event interval 0.12 s, Gaussian picker, and the catalogue filter COV_Err_XYZ <= 0.15 km and COA >= 5.5. The manuscript
code itself is not run; only its input files (station list) and published outputs (catalogues) are read.

``report`` gates the three counts of Section 4.1 (triggered, located, filtered); it refuses to run before the locate
stage has finished. It also reports, not gated, the fraction of the published located catalogue matched by an event of
this run within 0.1 s.
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

PUBLISHED = Path(__file__).with_name("published_quakemigrate.json")
EXAMPLE = "figures/icequake_example/generate_results"
START, END = "2009-01-20T00:00:00.0", "2009-01-22T00:00:00.0"
BLOCKS = (("2009-01-19T23:59:30.0", "2009-01-19T23:59:59.999"), ("2009-01-20T00:00:00.0", "2009-01-20T23:59:59.999"),
          ("2009-01-21T00:00:00.0", "2009-01-21T23:59:59.999"), ("2009-01-22T00:00:00.0", "2009-01-22T00:00:30.0"))
RUN = "paper_run"
FILTER = "COV_Err_XYZ <= 0.15 and COA >= 5.5"
MATCH_S = 0.1


def _paths(work):
    w = Path(work)
    return dict(stations=w / "inputs" / "rutford_stations.txt", mseed=w / "inputs" / "mSEED",
                lut=w / "outputs" / "lut" / "icequake.LUT", runs=w / "outputs" / "runs")


def _single_thread():
    os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1", MKL_NUM_THREADS="1")


def fetch(work, published=PUBLISHED, download=True):
    """clone the manuscript repository at the pinned commit (data only), copy the station list and the published
    catalogues into work, then download the waveforms as 0_get_data.py does."""
    pub = json.load(open(published))
    m = pub["manuscript"]
    w, p = Path(work), _paths(work)
    repo = w / "manuscript"
    if not repo.exists():
        subprocess.run(["git", "clone", "-q", m["repo"], str(repo)], check=True)
    subprocess.run(["git", "-C", str(repo), "checkout", "-q", m["commit"]], check=True)
    ex = repo / EXAMPLE
    p["stations"].parent.mkdir(parents=True, exist_ok=True)
    p["stations"].write_bytes((ex / "inputs" / "rutford_stations.txt").read_bytes())
    for name in ("rutford_icequakes.csv", "rutford_icequakes_gc150_coa55.csv"):
        (w / f"published_{name}").write_bytes((ex / name).read_bytes())
    if not download:
        return
    from obspy import UTCDateTime
    from obspy.clients.fdsn.mass_downloader import GlobalDomain, MassDownloader, Restrictions
    from quakemigrate.io import read_stations

    def storage(network, station, location, channel, starttime, endtime):
        return (p["mseed"] / f"{starttime.year}" / f"{starttime.julday:03d}" / f"{station}_{channel[2]}.m").as_posix()

    names = ",".join(read_stations(str(p["stations"]))["Name"])
    for a, b in BLOCKS:
        r = Restrictions(starttime=UTCDateTime(a), endtime=UTCDateTime(b), chunklength_in_sec=86400, network="YG",
                         station=names, channel_priorities=["GL[123]"], minimum_interstation_distance_in_m=0,
                         reject_channels_with_gaps=False)
        MassDownloader(providers=["IRIS"]).download(GlobalDomain(), r, threads_per_client=3, mseed_storage=storage,
                                                    stationxml_storage=str(w / "inputs" / "DATALESS"))
    missing = missing_waveforms(p["mseed"], read_stations(str(p["stations"]))["Name"])
    if missing:
        raise SystemExit(f"{len(missing)} waveform files missing after download (rerun fetch): {missing[:5]}")


def missing_waveforms(mseed, names):
    """expected YEAR/JD/STATION_C.m files (each station, components 1-3, each day of BLOCKS) that are absent; the
    mass downloader skips a channel on a server error without failing."""
    from obspy import UTCDateTime
    days = sorted({(UTCDateTime(a).year, UTCDateTime(a).julday) for a, _b in BLOCKS})
    want = [Path(mseed) / f"{y}" / f"{d:03d}" / f"{n}_{c}.m" for y, d in days for n in names for c in "123"]
    return [str(f) for f in want if not f.exists()]


def lut(work):
    from obspy.core import AttribDict
    from pyproj import Proj
    from quakemigrate.io import read_stations
    from quakemigrate.lut import compute_traveltimes
    p = _paths(work)
    p["lut"].parent.mkdir(parents=True, exist_ok=True)
    g = AttribDict()
    g.ll_corner, g.ur_corner, g.node_spacing = [-84.1, -78.19, 1.0], [-83.7, -78.1, 2.5], [0.05, 0.05, 0.05]
    g.grid_proj = Proj(proj="lcc", units="km", lon_0=-83.9, lat_0=-78.145, lat_1=-78.16, lat_2=-78.13, datum="WGS84",
                       ellps="WGS84", no_defs=True)
    g.coord_proj = Proj(proj="longlat", datum="WGS84", ellps="WGS84", no_defs=True)
    compute_traveltimes(g, read_stations(str(p["stations"])), method="homogeneous", phases=["P", "S"], vp=3.841,
                        vs=1.970, log=True, save_file=str(p["lut"]))


def _onset(position):
    from quakemigrate.signal.onsets import STALTAOnset
    o = STALTAOnset(position=position, sampling_rate=250)
    o.phases = ["P", "S"]
    o.bandpass_filters = {"P": [20, 124, 4], "S": [10, 124, 4]}
    o.sta_lta_windows = {"P": [0.01, 0.25], "S": [0.05, 0.5]}
    o.channel_maps = {"P": "*1" if position == "classic" else "*[1]", "S": "*[2,3]"}
    return o


def _scan(work, position, decimate, picker=False):
    from quakemigrate import QuakeScan
    from quakemigrate.io import Archive, read_lut, read_stations
    p = _paths(work)
    st = read_stations(str(p["stations"]))
    arch = Archive(archive_path=str(p["mseed"]), stations=st, archive_format="YEAR/JD/STATION")
    lt = read_lut(lut_file=str(p["lut"]))
    if decimate:
        lt.decimate([2, 2, 2], inplace=True)
    onset = _onset(position)
    kw = {}
    if picker:
        from quakemigrate.signal.pickers import GaussianPicker
        kw["picker"] = GaussianPicker(onset=onset)
        kw["picker"].plot_picks = False
    return QuakeScan(arch, lt, onset=onset, run_path=str(p["runs"]), run_name=RUN, log=True, loglevel="info", **kw)


def detect(work, threads=12):
    _single_thread()
    from quakemigrate.io import read_stations
    p = _paths(work)
    missing = missing_waveforms(p["mseed"], read_stations(str(p["stations"]))["Name"])
    if missing:
        raise SystemExit(f"{len(missing)} waveform files missing; rerun fetch: {missing[:5]}")
    s = _scan(work, "classic", decimate=True)
    s.timestep = 30.0
    s.threads = threads
    s.detect(START, END)


def trigger(work):
    _single_thread()
    from quakemigrate import Trigger
    from quakemigrate.io import read_lut
    p = _paths(work)
    t = Trigger(read_lut(lut_file=str(p["lut"])), run_path=str(p["runs"]), run_name=RUN, output_all_columns=True,
                log=True, loglevel="info")
    t.marginal_window, t.min_event_interval, t.normalise_coalescence = 0.06, 0.12, True
    t.write_event_time_windows, t.threshold_method, t.static_threshold = True, "static", 3.0
    t.trigger(START, END, interactive_plot=False)


def locate(work, threads=4):
    _single_thread()
    s = _scan(work, "centred", decimate=False, picker=True)
    s.marginal_window = 0.06
    s.threads = threads
    s.plot_event_summary = False
    s.write_cut_waveforms = False
    s.locate(starttime=START, endtime=END)


def counts(work):
    """(triggered, located catalogue, filtered catalogue) of the run in work."""
    import pandas as pd
    run = _paths(work)["runs"] / RUN
    trig = sorted((run / "trigger" / "events").glob("*_TriggeredEvents.csv"))
    ev = sorted((run / "locate" / "events").glob("*.event"))
    if not trig:
        raise SystemExit(f"{run}: no trigger output; run the trigger stage")
    if not ev:
        raise SystemExit(f"{run}: no located events; the report needs the locate stage (Section 4.1 counts)")
    n_trig = sum(len(pd.read_csv(f)) for f in trig)
    cat = pd.concat([pd.read_csv(f) for f in ev])
    return n_trig, cat, cat.query(FILTER)


def matched_fraction(published_dt, run_dt, tol=MATCH_S):
    """fraction of published origin times with an origin time of the run within tol seconds."""
    import numpy as np
    a = np.sort(np.asarray(run_dt, float))
    b = np.asarray(published_dt, float)
    if not len(a) or not len(b):
        return 0.0
    i = np.clip(np.searchsorted(a, b), 1, len(a) - 1)
    d = np.minimum(np.abs(a[i] - b), np.abs(a[i - 1] - b))
    return float(np.mean(d <= tol))


def _seconds(dt):
    import pandas as pd
    return (pd.to_datetime(dt, utc=True) - pd.Timestamp("2009-01-20", tz="UTC")).dt.total_seconds()


def report(work, published, out):
    import pandas as pd
    import quakemigrate
    pub = json.load(open(published))
    if not pub.get("verified"):
        raise SystemExit(f"{published}: published counts not verified against the paper")
    n_trig, cat, filt = counts(work)
    val = {"triggered": n_trig, "located": len(cat), "filtered": len(filt)}
    cond = [{"name": k, "value": float(val[k]), "published": float(v)} for k, v in pub["metrics"].items()]
    ref = pd.read_csv(Path(work) / "published_rutford_icequakes.csv")
    frac = matched_fraction(_seconds(ref["DT"]), _seconds(cat["DT"]))
    ev = {"conditions": cond, "source": f"QuakeMigrate {quakemigrate.__version__} on the Rutford example vs "
          f"{pub['source'].split(',')[0]} et al. 2026, Section 4.1 (manuscript {pub['manuscript']['commit'][:7]})",
          "reported_only": {f"published_located_matched_within_{MATCH_S}s": frac},
          "quakemigrate": quakemigrate.__version__}
    Path(out).write_text(json.dumps(ev, indent=1))
    for c in cond:
        print(f"{c['name']:10s} {c['value']:8.0f}  published {c['published']:8.0f}  "
              f"rel {(c['value'] - c['published']) / c['published']:+.3f}")
    print(f"published located events matched within {MATCH_S} s: {frac:.3f}")


def main(argv=None):
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("fetch", "lut", "detect", "trigger", "locate", "report"):
        s = sub.add_parser(name)
        s.add_argument("--work", required=True)
        if name in ("detect", "locate"):
            s.add_argument("--threads", type=int, default=12 if name == "detect" else 4)
        if name == "fetch":
            s.add_argument("--no-download", action="store_true")
        if name == "report":
            s.add_argument("published", nargs="?", default=str(PUBLISHED))
            s.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "fetch":
        fetch(a.work, download=not a.no_download)
    elif a.cmd == "lut":
        lut(a.work)
    elif a.cmd in ("detect", "locate"):
        {"detect": detect, "locate": locate}[a.cmd](a.work, a.threads)
    elif a.cmd == "trigger":
        trigger(a.work)
    else:
        report(a.work, a.published, a.out)


if __name__ == "__main__":
    sys.exit(main())
