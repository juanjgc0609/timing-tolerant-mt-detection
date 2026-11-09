"""Quality control of real datasets: waveforms, station metadata and catalogues (data registry, protocol Section 2).

Every check here reads headers, metadata or catalogue fields only, unless ``samples=True`` is passed, so the checks
are allowed on sealed evaluation partitions (``ttmt.registry``: an evaluation partition may be inspected for QC
metadata and nothing else). With ``samples=True`` the waveform check also reads the samples to find flat lines
(dead or clipped channels); that decodes samples, so it calls ``ttmt.registry.check_read`` first and is refused on
a sealed window.

The checks and their thresholds (``LIMITS``) are the ones the data cards report:

waveforms   coverage of the requested span, gaps and overlaps, sampling-rate changes and mismatches with the
            metadata, miniSEED timing quality (blockette 1001), "clock locked" and "time tag questionable" flags,
            flat lines.
stations    StationXML present; one channel epoch covering the span; a response with an overall sensitivity
            whose units and value agree with the product of the stage gains; orientation present; a complete,
            orthogonal three-component set; response epochs that change inside the span.
catalogue   magnitude of completeness (maximum curvature plus 0.2, and goodness of fit at 90 and 95 %),
            b-value (Aki-Utsu with the binning correction), completeness in moving windows (short-term
            incompleteness after large events), location-quality quantiles, close pairs in time (overlapping
            events) and duplicates.
"""
import math
from collections import Counter, defaultdict
from dataclasses import dataclass, field, asdict

import numpy as np

LIMITS = {
    "coverage_min": 0.98,          # fraction of the span with data, per channel
    "gap_max_s": 60.0,             # longest single gap tolerated in a usable channel
    "timing_quality_min": 70.0,    # mean miniSEED timing quality (0-100) when the channel reports it
    "clock_unlocked_max": 0.05,    # fraction of records without the "clock locked" flag, when it is reported
    "sensitivity_rel_tol": 0.05,   # |overall sensitivity / product of stage gains - 1|
    "orthogonality_tol_deg": 5.0,  # horizontals at 90 deg from each other, vertical at dip -90 or 90
    "flatline_s": 10.0,            # a run of identical samples this long counts as flat
    "flatline_max": 0.01,          # fraction of samples in flat runs tolerated
}


# ================================================================================================= catalogue
def _bins(m, dm):
    """magnitudes rounded to the nearest bin centre of width dm."""
    return np.round(np.round(np.asarray(m, float) / dm) * dm, 6)


def mc_maxc(mags, dm=0.1, correction=0.2):
    """maximum curvature (Wiemer and Wyss 2000) plus the usual 0.2 correction (Woessner and Wiemer 2005)."""
    m = _bins(mags, dm)
    if m.size < 20:
        return float("nan")
    vals, counts = np.unique(m, return_counts=True)
    return float(vals[np.argmax(counts)] + correction)


def b_value(mags, mc, dm=0.1):
    """Aki (1965) with Utsu's binning correction, and Shi and Bolt (1982) standard error."""
    m = _bins(mags, dm)
    m = m[m >= mc - 1e-9]
    if m.size < 20:
        return float("nan"), float("nan"), int(m.size)
    mean = m.mean()
    b = math.log10(math.e) / (mean - (mc - dm / 2))
    se = 2.3 * b * b * math.sqrt(np.sum((m - mean) ** 2) / (m.size * (m.size - 1)))
    return float(b), float(se), int(m.size)


def mc_gft(mags, dm=0.1, levels=(95.0, 90.0)):
    """goodness-of-fit Mc (Wiemer and Wyss 2000): the lowest Mc at which a Gutenberg-Richter law fitted above it
    explains the observed cumulative counts B with residual R = 100 - 100 * sum|B - S| / sum B at or above the level.
    Returns {level: Mc or nan}."""
    m = _bins(mags, dm)
    out = {lv: float("nan") for lv in levels}
    if m.size < 50:
        return out
    lo, hi = float(np.min(m)), float(np.max(m))
    for mc in np.round(np.arange(lo, hi - 1.0 + 1e-9, dm), 6):
        sel = m[m >= mc - 1e-9]
        if sel.size < 50:
            break
        b, _, n = b_value(sel, mc, dm)
        a = math.log10(n) + b * mc
        edges = np.round(np.arange(mc, float(sel.max()) + dm / 2, dm), 6)
        obs = np.array([np.sum(sel >= e - 1e-9) for e in edges], float)       # cumulative counts
        syn = 10 ** (a - b * edges)
        R = 100.0 - 100.0 * np.sum(np.abs(obs - syn)) / np.sum(obs)     # cumulative counts, as in ZMAP
        for lv in levels:
            if math.isnan(out[lv]) and R >= lv:
                out[lv] = float(mc)
    return out


def mc_moving(t, mags, n=250, step=50, dm=0.1):
    """Mc (maximum curvature + 0.2) in moving windows of n consecutive events: (start time, end time, Mc) rows.
    A large spread means the completeness changes with time (short-term incompleteness after large events)."""
    order = np.argsort(t)
    t, mags = np.asarray(t, float)[order], np.asarray(mags, float)[order]
    rows = []
    for a in range(0, max(1, t.size - n + 1), step):
        b = min(a + n, t.size)
        if b - a < n:
            break
        rows.append((float(t[a]), float(t[b - 1]), mc_maxc(mags[a:b], dm)))
    return rows


def close_pairs(t, x_km=None, y_km=None, dts=(0.5, 1.0, 2.0), r_km=None):
    """number of events whose next event (in origin time) follows within dt s, for each dt; with coordinates and
    r_km, only pairs closer than r_km count. Overlapping events are the target regime of the project."""
    order = np.argsort(t)
    t = np.asarray(t, float)[order]
    d = np.diff(t)
    if r_km is not None and x_km is not None:
        x, y = np.asarray(x_km, float)[order], np.asarray(y_km, float)[order]
        near = np.hypot(np.diff(x), np.diff(y)) <= r_km
    else:
        near = np.ones_like(d, bool)
    return {float(dt): int(np.sum((d < dt) & near)) for dt in dts}


def duplicates(t, lat, lon, dt=0.1, r_km=1.0):
    """likely duplicate entries: consecutive events within dt s and r_km km."""
    order = np.argsort(t)
    t, lat, lon = (np.asarray(v, float)[order] for v in (t, lat, lon))
    x = (lon - lon.mean()) * 111.19 * np.cos(np.radians(lat.mean()))
    y = (lat - lat.mean()) * 111.19
    return int(np.sum((np.diff(t) < dt) & (np.hypot(np.diff(x), np.diff(y)) < r_km)))


def _q(v, qs=(0.5, 0.9)):
    v = np.asarray([x for x in v if x is not None and not (isinstance(x, float) and math.isnan(x))], float)
    return {f"p{int(q * 100)}": (float(np.quantile(v, q)) if v.size else None) for q in qs} | {"n": int(v.size)}


def catalogue_qc(events, dm=0.1):
    """events: dicts with t (s), lat, lon, dep, M and, when known, erh, erz, rms, gap, nst, mag_type, status.
    Returns a JSON-able summary. It holds aggregates only, so it may be computed on a sealed catalogue."""
    if not events:
        return {"n": 0}
    t = np.array([e["t"] for e in events], float)
    M = np.array([e["M"] for e in events], float)
    lat = np.array([e["lat"] for e in events], float)
    lon = np.array([e["lon"] for e in events], float)
    x = (lon - lon.mean()) * 111.19 * np.cos(np.radians(lat.mean()))
    y = (lat - lat.mean()) * 111.19
    mc = mc_maxc(M, dm)
    gft = mc_gft(M, dm)
    b, bse, nb = b_value(M, mc, dm) if not math.isnan(mc) else (float("nan"),) * 2 + (0,)
    mov = mc_moving(t, M, dm=dm) if t.size >= 250 else []
    mcs = [r[2] for r in mov]
    return {
        "n": int(t.size),
        "span_s": float(t.max() - t.min()),
        "M_min": float(M.min()), "M_max": float(M.max()),
        "mag_types": dict(Counter(e.get("mag_type") for e in events)),
        "status": dict(Counter(e.get("status") for e in events)),
        "mc_maxc": mc, "mc_gft95": gft[95.0], "mc_gft90": gft[90.0],
        "b": b, "b_se": bse, "n_above_mc": nb,
        "mc_moving_min": (float(min(mcs)) if mcs else None), "mc_moving_max": (float(max(mcs)) if mcs else None),
        "erh_km": _q([e.get("erh") for e in events]), "erz_km": _q([e.get("erz") for e in events]),
        "rms_s": _q([e.get("rms") for e in events]), "gap_deg": _q([e.get("gap") for e in events]),
        "nst": _q([e.get("nst") for e in events], (0.1, 0.5)),
        "next_within_s": close_pairs(t, x, y),
        "next_within_s_10km": close_pairs(t, x, y, r_km=10.0),
        "duplicates": duplicates(t, lat, lon),
    }


# ================================================================================================= waveforms
@dataclass
class ChannelQC:
    nslc: str
    span: tuple                       # requested (start, end) as ISO strings
    coverage: float = 0.0
    n_gaps: int = 0
    n_overlaps: int = 0
    max_gap_s: float = 0.0
    sampling_rates: list = field(default_factory=list)
    timing_quality_mean: float = None
    clock_unlocked_frac: float = None
    time_tag_questionable: int = None
    flatline_frac: float = None
    issues: list = field(default_factory=list)

    @property
    def usable(self):
        return not self.issues


def _flag_stats(path):
    """miniSEED flags of one file (obspy.io.mseed.util.get_flags); None when the file is not miniSEED."""
    try:
        from obspy.io.mseed.util import get_flags
        return get_flags(path, io_flags=True, activity_flags=False, data_quality_flags=True, timing_quality=True)
    except Exception:
        return None


def _flat_fraction(tr, run_s):
    x = np.asarray(tr.data)
    if x.size < 2:
        return 1.0
    n = int(round(run_s * tr.stats.sampling_rate))
    same = np.concatenate([[False], np.diff(x) == 0])
    # lengths of runs of equal consecutive samples
    edges = np.flatnonzero(np.diff(np.concatenate([[0], same.astype(int), [0]])))
    starts, ends = edges[::2], edges[1::2]
    flat = sum(e - s + 1 for s, e in zip(starts, ends) if e - s + 1 >= n)
    return float(flat / x.size)


def waveform_qc(paths, t0, t1, samples=False, limits=LIMITS, metadata_rate=None, coords=None):
    """QC of the channels in the miniSEED files ``paths`` over [t0, t1) (obspy UTCDateTime). One ChannelQC per
    NSLC. Without ``samples`` only headers are read (allowed on sealed partitions). With ``samples`` every channel
    passes ``registry.check_read`` before its samples are decoded; ``coords`` {"NET.STA": (lat, lon)} places the
    stations (a station without coordinates counts as inside every sealed area)."""
    from obspy import Stream, read
    from .registry import check_read
    if samples:
        for p in paths:
            for tr in read(str(p), headonly=True, starttime=t0, endtime=t1, nearest_sample=False):
                lat, lon = (coords or {}).get(f"{tr.stats.network}.{tr.stats.station}", (None, None))
                check_read(tr.stats.network, lat, lon, float(t0.timestamp), float(t1.timestamp), purpose="qc-samples")
    st = Stream()
    flags = defaultdict(list)
    for p in paths:
        s = read(str(p), headonly=not samples, starttime=t0, endtime=t1, nearest_sample=False)
        st += s
        f = _flag_stats(str(p))
        if f is not None and len(s):
            flags[s[0].id].append(f)
    span = float(t1 - t0)
    out = {}
    for nslc in sorted({tr.id for tr in st}):
        trs = sorted(st.select(id=nslc), key=lambda tr: tr.stats.starttime)
        q = ChannelQC(nslc, (str(t0), str(t1)))
        q.sampling_rates = sorted({float(tr.stats.sampling_rate) for tr in trs})
        # coverage from the union of the trace intervals clipped to the span
        iv = sorted((max(float(tr.stats.starttime - t0), 0.0),
                     min(float(tr.stats.endtime - t0) + tr.stats.delta, span)) for tr in trs)
        covered, cur_a, cur_b, gaps, overlaps = 0.0, None, None, [], 0
        for a, b in iv:
            if b <= a:
                continue
            if cur_b is None:
                cur_a, cur_b = a, b
                if a > 0:
                    gaps.append(a)
            elif a > cur_b + 1e-6:
                covered += cur_b - cur_a
                gaps.append(a - cur_b)
                cur_a, cur_b = a, b
            else:
                if a < cur_b - 1.5 / max(q.sampling_rates):
                    overlaps += 1
                cur_b = max(cur_b, b)
        if cur_b is not None:
            covered += cur_b - cur_a
            if cur_b < span - 1e-6:
                gaps.append(span - cur_b)
        q.coverage = covered / span if span > 0 else 0.0
        q.n_gaps, q.n_overlaps = len(gaps), overlaps
        q.max_gap_s = float(max(gaps)) if gaps else 0.0
        fl = flags.get(nslc, [])
        tq = [f["timing_quality"]["mean"] for f in fl if f.get("timing_quality")]
        if tq:
            q.timing_quality_mean = float(np.mean(tq))
        dq = [f.get("data_quality_flags_counts", {}) for f in fl]
        io = [f.get("io_and_clock_flags_counts", {}) for f in fl]
        nrec = sum(f.get("record_count", 0) for f in fl) if fl else 0
        if io and nrec:
            locked = sum(d.get("clock_locked", 0) for d in io)
            # networks that never set the flag report zero locked records; count it only if any record sets it
            q.clock_unlocked_frac = (1.0 - locked / nrec) if locked else None
        if dq:
            q.time_tag_questionable = int(sum(d.get("time_tag_questionable", 0) for d in dq))
        if samples:
            n_all = sum(tr.stats.npts for tr in trs)
            q.flatline_frac = (float(sum(_flat_fraction(tr, limits["flatline_s"]) * tr.stats.npts for tr in trs)
                                     / n_all) if n_all else 1.0)
        # issues
        if q.coverage < limits["coverage_min"]:
            q.issues.append(f"coverage {q.coverage:.3f}")
        if q.max_gap_s > limits["gap_max_s"]:
            q.issues.append(f"gap {q.max_gap_s:.0f} s")
        if len(q.sampling_rates) > 1:
            q.issues.append(f"sampling rate changes {q.sampling_rates}")
        if metadata_rate and nslc in metadata_rate and q.sampling_rates and \
                abs(q.sampling_rates[0] - metadata_rate[nslc]) > 1e-6:
            q.issues.append(f"rate {q.sampling_rates[0]} != metadata {metadata_rate[nslc]}")
        if q.timing_quality_mean is not None and q.timing_quality_mean < limits["timing_quality_min"]:
            q.issues.append(f"timing quality {q.timing_quality_mean:.0f}")
        if q.clock_unlocked_frac is not None and q.clock_unlocked_frac > limits["clock_unlocked_max"]:
            q.issues.append(f"clock unlocked {q.clock_unlocked_frac:.2f}")
        if q.time_tag_questionable:
            q.issues.append(f"time tag questionable in {q.time_tag_questionable} records")
        if q.flatline_frac is not None and q.flatline_frac > limits["flatline_max"]:
            q.issues.append(f"flat {q.flatline_frac:.3f}")
        out[nslc] = q
    return out


# ================================================================================================= stations
@dataclass
class StationQC:
    station: str                      # NET.STA
    channels: dict = field(default_factory=dict)   # NSLC -> list of issues ([] = fine)
    three_component: list = field(default_factory=list)   # NET.STA.LOC.BAND+INSTR codes with a full 3-C set
    issues: list = field(default_factory=list)

    @property
    def usable(self):
        return bool(self.three_component) and not self.issues


def _check_channel(cha, t0, t1, limits):
    iss = []
    resp = cha.response
    if resp is None or resp.instrument_sensitivity is None:
        return ["no response"]
    s = resp.instrument_sensitivity
    if not s.value:
        iss.append("sensitivity zero")
    units = (s.input_units or "").upper()
    if units not in ("M/S", "M/S**2", "M/S2", "M", "COUNTS", "PA", "V"):
        iss.append(f"input units {s.input_units!r}")
    gains = [st.stage_gain for st in resp.response_stages if st.stage_gain]
    if gains and s.value:
        prod = float(np.prod(gains))
        if abs(prod / s.value - 1.0) > limits["sensitivity_rel_tol"]:
            iss.append(f"sensitivity {s.value:.4g} != product of stage gains {prod:.4g}")
    elif not resp.response_stages:
        iss.append("no response stages")
    if cha.azimuth is None or cha.dip is None:
        iss.append("orientation missing")
    if cha.sample_rate in (None, 0):
        iss.append("sample rate missing")
    return iss


def station_qc(inventory, t0, t1, limits=LIMITS):
    """QC of the StationXML ``inventory`` (obspy Inventory) for the span [t0, t1). One StationQC per station."""
    out = {}
    for net in inventory:
        for sta in net:
            q = StationQC(f"{net.code}.{sta.code}")
            groups = defaultdict(list)
            epochs = defaultdict(list)
            for cha in sta:
                end = cha.end_date or t1
                if cha.start_date > t1 or end < t0:
                    continue
                nslc = f"{net.code}.{sta.code}.{cha.location_code}.{cha.code}"
                epochs[nslc].append(cha)
            for nslc, chs in epochs.items():
                iss = []
                covering = [c for c in chs if c.start_date <= t0 and (c.end_date is None or c.end_date >= t1)]
                if not covering:
                    iss.append(f"{len(chs)} epochs inside the span (metadata change mid-span)")
                    cha = chs[-1]
                else:
                    cha = covering[0]
                iss += _check_channel(cha, t0, t1, limits)
                q.channels[nslc] = iss
                loc, code = nslc.split(".")[2], nslc.split(".")[3]
                groups[(loc, code[:2])].append((code[2], cha))
            for (loc, bi), comps in groups.items():
                codes = {c for c, _ in comps}
                vert = [ch for c, ch in comps if c == "Z"]
                hor = [ch for c, ch in comps if c in "NE12"]
                if not vert or len(hor) < 2:
                    continue
                ok = all(not q.channels[f"{net.code}.{sta.code}.{loc}.{bi}{c}"] for c in codes if c in "ZNE12")
                az = sorted(float(h.azimuth or 0.0) % 180 for h in hor[:2])
                orth = abs((az[1] - az[0]) - 90.0) <= limits["orthogonality_tol_deg"]
                vdip = abs(abs(float(vert[0].dip or 0.0)) - 90.0) <= limits["orthogonality_tol_deg"]
                tag = f"{net.code}.{sta.code}.{loc}.{bi}"
                if ok and orth and vdip:
                    q.three_component.append(tag)
                elif not (orth and vdip):
                    q.issues.append(f"{tag}: components not orthogonal (azimuths {az}, dip {vert[0].dip})")
            if not q.three_component and not q.issues:
                q.issues.append("no complete three-component set with valid metadata")
            out[q.station] = q
    return out


def report_dict(obj):
    """JSON-able form of a QC result (dict of ChannelQC / StationQC, or a catalogue summary)."""
    if isinstance(obj, dict):
        return {k: report_dict(v) for k, v in obj.items()}
    if isinstance(obj, (ChannelQC, StationQC)):
        return asdict(obj) | {"usable": obj.usable}
    return obj
