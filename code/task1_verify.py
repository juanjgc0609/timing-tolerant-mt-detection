"""
task1_verify.py — verifica la matematica de polarizacion ANTES de tocar dato real.

Construye ondas sinteticas con back-azimut conocido y comprueba que el estimador de
`task1_analyze.polarization_at` lo recupera. Si esto falla, el numero de la Tarea 1 no vale y no
habria forma de saberlo mirando dato real, porque ahi no hay verdad conocida.

    python task1_verify.py
"""
import sys
import numpy as np
from obspy import Trace, Stream, UTCDateTime

sys.path.insert(0, ".")
from task1_analyze import polarization_at          # noqa: E402

FS = 100.0
T0 = UTCDateTime("2024-01-01T00:00:00")
results = []


def check(name, ok, detail=""):
    results.append((name, ok))
    print(f"   [{'PASS' if ok else 'FAIL'}] {name}" + (f"  — {detail}" if detail else ""))
    return ok


def synth_station(baz_deg, inc_deg=35.0, rect=1.0, sta="TST", noise=0.0, seed=0):
    """
    Onda P rectilinea llegando con back-azimut `baz_deg` medido desde el norte.
    Vector de movimiento: horizontal a lo largo del back-azimut, vertical segun la incidencia.
    """
    rng = np.random.default_rng(seed)
    n = int(8 * FS)
    t = np.arange(n) / FS
    fc = 2.0
    x = np.exp(-((t - 4.0) ** 2) / (2 * 0.25 ** 2)) * np.sin(2 * np.pi * fc * (t - 4.0))
    b = np.radians(baz_deg)
    i = np.radians(inc_deg)
    e = np.sin(b) * np.sin(i)
    nn = np.cos(b) * np.sin(i)
    z = np.cos(i)
    comps = {"E": e * x, "N": nn * x, "Z": z * x}
    if rect < 1.0:                       # añade una segunda direccion ortogonal
        amp = np.sqrt(max(0.0, 1.0 - rect))
        y = np.roll(x, 7) * amp
        comps["E"] = comps["E"] + np.cos(b) * y
        comps["N"] = comps["N"] - np.sin(b) * y
    st = Stream()
    for k, v in comps.items():
        v = v + (rng.normal(scale=noise, size=n) if noise > 0 else 0.0)
        tr = Trace(data=v.astype(float))
        tr.stats.station = sta
        tr.stats.channel = "HH" + k
        tr.stats.sampling_rate = FS
        tr.stats.starttime = T0
        st += tr
    return st


def fold(a, b):
    d = abs((a % 180.0) - (b % 180.0))
    return min(d, 180.0 - d)


def main():
    print("=" * 68)
    print("VERIFICACION DE LA POLARIZACION (sobre sintetico, antes del dato real)")
    print("=" * 68)

    print("\n1 — recuperacion del back-azimut, sin ruido")
    errs = []
    for baz in (0., 23., 47., 90., 130., 175., 210., 300.):
        st = synth_station(baz)
        out = polarization_at(st, "TST", T0 + 4.0, 0.6)
        if out is None:
            check(f"baz={baz:.0f}", False, "sin salida")
            continue
        rect, az, inc, _n = out
        e = fold(az, baz)
        errs.append(e)
        check(f"baz={baz:6.1f} -> medido {az:6.1f}", e <= 2.0, f"error {e:.2f} grados")
    if errs:
        print(f"   error maximo {max(errs):.2f} grados")

    print("\n2 — rectilinearidad")
    for r_true in (1.0, 0.7, 0.4):
        st = synth_station(60.0, rect=r_true)
        out = polarization_at(st, "TST", T0 + 4.0, 0.6)
        rect = out[0] if out else np.nan
        check(f"rect objetivo {r_true:.1f} -> medida {rect:.2f}",
              (not np.isnan(rect)) and rect > (0.3 if r_true < 1.0 else 0.9))

    print("\n3 — degradacion con ruido (debe ser gradual, no abrupta)")
    prev = 0.0
    for nz in (0.0, 0.02, 0.05, 0.10, 0.20):
        es = []
        for s in range(8):
            st = synth_station(75.0, noise=nz, seed=s)
            out = polarization_at(st, "TST", T0 + 4.0, 0.6)
            if out:
                es.append(fold(out[1], 75.0))
        m = float(np.median(es)) if es else np.nan
        print(f"   ruido {nz:4.2f} -> error mediano {m:5.2f} grados")
        prev = m
    check("error acotado con ruido 0.20", prev < 25.0, f"{prev:.1f} grados")

    print("\n4 — ambiguedad de 180 grados declarada, no un fallo")
    st = synth_station(40.0)
    out = polarization_at(st, "TST", T0 + 4.0, 0.6)
    az = out[1]
    check("el estimador devuelve azimut plegado a [0,180)", 0.0 <= az < 180.0,
          f"az={az:.1f}. El signo de un autovector es arbitrario, por eso todas las "
          f"comparaciones de la Tarea 1 se hacen plegadas a [0,90].")

    nf = sum(1 for _n, ok in results if not ok)
    print("\n" + "=" * 68)
    if nf == 0:
        print("TODO PASA. El estimador de polarizacion recupera el back-azimut conocido.")
        print("Puedes interpretar el numero que salga de task1_decide.py.")
    else:
        print(f"{nf} COMPROBACION(ES) FALLAN. NO ejecutes la Tarea 1 sobre dato real:")
        for n, ok in results:
            if not ok:
                print(f"   {n}")
    print("=" * 68)
    sys.exit(0 if nf == 0 else 1)


if __name__ == "__main__":
    main()
