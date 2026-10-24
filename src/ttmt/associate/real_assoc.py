"""REAL (Zhang et al., 2019; Dal-mzhang/REAL, REAL.c of the Aug. 2021 version) on picks in the atlas frame.

REAL is a C program that reads picks per station from files and works in geographic coordinates. The wrapper maps the
local frame to degrees around (0, 0) with REAL's own scale (x = lon * 111.19 km, y = lat * 111.19 km; -D latitude
centre 0, so a degree is the same distance on both axes), writes station.dat, the pick files (seconds from the record
start, which REAL takes as seconds of its day; REAL accepts times up to one month), and the 1-D table ttdb.txt of the
frame's model (velocity.Velocity1D.real_ttdb), runs the binary in a fresh folder with OMP_NUM_THREADS = 1, and reads
phase_sel.txt (events after REAL's second selection, with their picks). Location is REAL's grid location of
phase_sel.txt, not the simulated-annealing relocation (hypolocSA.dat), because only the former comes with the picks.

Binary: TTMT_REAL_BIN or Frame.fixed["real"]["binary"] (built with gcc -O2 -fopenmp REAL.c -lm).

Tuned (space.REAL): as_np0 -> np0 (minimum P picks), as_nps0 -> nps0 (minimum P + S picks), as_tint -> tint (REAL's
minimum interval between two events, s; log-uniform in [0.1, 5] s by ATLAS-0 amendment 1, ADR-0009). Fixed: the -S defaults of REAL's synthetic demo
(ns0 2, npsboth0 2, std0 0.5 s, dtps 0, nrt 1.5), a -R search radius of 25 km around the station of the initiating
pick with 2 km grid steps, and the -G table grid (0.01 deg, 1 km). The run folder is deleted after reading.
"""
import os
import shutil
import subprocess
import tempfile

import numpy as np

from .base import detections, empty_result
from .velocity import KM_PER_DEG

FIXED = {"binary": None, "ns0": 2, "npsboth0": 2, "std0": 0.5, "dtps": 0.0, "nrt": 1.5, "drt": 0.5, "nxd": 1.0,
         "rsel": 5.0, "rx_km": 25.0, "tdx_km": 2.0, "tdh_km": 2.0, "V": None, "table_dx_deg": 0.01, "table_dz_km": 1.0, "one_d": True,
         "timeout": 3600}


def binary(fx):
    b = fx["binary"] or os.environ.get("TTMT_REAL_BIN")
    if not b or not os.path.exists(b):
        raise FileNotFoundError("REAL binary not found: set TTMT_REAL_BIN or Frame.fixed['real']['binary']")
    return b


def associate(frame, picks, cfg):
    fx = frame.fixed_for("real", FIXED)
    if len(picks) < int(cfg["as_nps0"]):
        return empty_result(picks)
    codes = frame.codes
    xyz = frame.stations.xyz
    run = tempfile.mkdtemp(prefix="real_", dir=frame.workdir)
    pkdir = os.path.join(run, "pk")
    os.makedirs(pkdir)
    with open(os.path.join(run, "station.dat"), "w") as fh:
        for c, (x, y, _z) in zip(codes, xyz):
            fh.write(f"{x / KM_PER_DEG:.6f} {y / KM_PER_DEG:.6f} XX {c} HHZ 0.0\n")
    for i, c in enumerate(codes):
        for ph, name in ((0, "P"), (1, "S")):
            k = np.flatnonzero((picks.sta == i) & (picks.phase == ph))
            k = k[np.argsort(picks.t[k], kind="stable")]
            with open(os.path.join(pkdir, f"XX.{c}.{name}.txt"), "w") as fh:
                for j in k:
                    fh.write(f"{picks.t[j]:.4f} {picks.p[j]:.4f} 0.0\n")
    # search box around the station of the initiating pick (REAL searches near the first-picked station)
    rx = fx["rx_km"] / KM_PER_DEG
    rh = float(frame.zlim[1])
    tdx, tdh = fx["tdx_km"] / KM_PER_DEG, fx["tdh_km"]
    R = f"-R{rx:.4f}/{rh:g}/{tdx:.5f}/{tdh:g}/{float(cfg['as_tint']):.4f}"
    vel = fx["V"] or (*frame.velocity.mean_velocities(rh), *frame.velocity.mean_velocities(1.0))
    V = "-V{:.4f}/{:.4f}/{:.4f}/{:.4f}/0".format(*vel)
    S = (f"-S{int(cfg['as_np0'])}/{fx['ns0']}/{int(cfg['as_nps0'])}/{fx['npsboth0']}/{fx['std0']:g}/{fx['dtps']:g}/"
         f"{fx['nrt']:g}/{fx['drt']:g}/{fx['nxd']:g}/{fx['rsel']:g}/0")
    args = [binary(fx), "-D2000/01/01/0.0", R, S, V]
    if fx["one_d"]:
        xmax = rx + frame.diagonal() / KM_PER_DEG + 0.05
        key = abs(hash((frame.velocity.layers, frame.velocity.vpvs, round(xmax, 3), rh + 2 * tdh,
                        fx["table_dx_deg"], fx["table_dz_km"])))
        tab = os.path.join(frame.workdir, f"real_ttdb_{key:x}.txt")
        if not os.path.exists(tab):
            frame.velocity.real_ttdb(tab + ".tmp", xmax, rh + 2 * tdh, fx["table_dx_deg"], fx["table_dz_km"])
            os.replace(tab + ".tmp", tab)
        shutil.copyfile(tab, os.path.join(run, "ttdb.txt"))
        args += [f"-G{xmax:.3f}/{rh + 2 * tdh:g}/{fx['table_dx_deg']:g}/{fx['table_dz_km']:g}"]
    args += ["station.dat", "pk"] + (["ttdb.txt"] if fx["one_d"] else [])
    env = dict(os.environ, OMP_NUM_THREADS="1")
    try:
        subprocess.run(args, cwd=run, env=env, check=True, capture_output=True, timeout=fx["timeout"])
        return read_phase_sel(os.path.join(run, "phase_sel.txt"), picks, codes)
    finally:
        shutil.rmtree(run, ignore_errors=True)


def read_phase_sel(path, picks, codes):
    """events and their picks from REAL's phase_sel.txt; a pick line is matched to the input pick of the same
    station and phase whose time is nearest (REAL prints 4 decimals)."""
    lut = {c: i for i, c in enumerate(codes)}
    rows, assign = [], []
    if not os.path.exists(path):
        return empty_result(picks)
    for line in open(path):
        f = line.split()
        if not f:
            continue
        if len(f) >= 17 and f[0].isdigit():         # num, otime (4 tokens), atime, std, lat, lon, dep, ...
            atime, lat, lon, dep = float(f[5]), float(f[7]), float(f[8]), float(f[9])
            rows.append((lon * KM_PER_DEG, lat * KM_PER_DEG, dep, atime))
        else:                                        # net sta phase abs_pk pk amp res weig baz
            i, ph, t = lut[f[1]], 0 if f[2].upper().startswith("P") else 1, float(f[3])
            k = np.flatnonzero((picks.sta == i) & (picks.phase == ph))
            if len(k):
                j = k[np.argmin(np.abs(picks.t[k] - t))]
                if abs(picks.t[j] - t) < 1e-3:
                    assign.append((j, len(rows) - 1))
    return detections(rows, picks, assign)
