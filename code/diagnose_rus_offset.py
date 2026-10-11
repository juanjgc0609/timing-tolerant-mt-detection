"""
diagnose_rus_offset.py — prueba si el desajuste de RUS es un DESFASE ANGULAR CONSTANTE (canales
N/E invertidos o rotados respecto al norte real -- error de orientacion de sensor, tipico cuando
se reinstala un equipo), en vez de dispersion aleatoria genuina.

Metodo: para cada estacion, barre un offset de 0 a 175 grados (paso 5) y ve con cual offset
'medido + offset' coincide mejor (menor desviacion circular mod 180) con el azimut verdadero,
usando SOLO los picks de control (asociados, rectilinearidad alta). Si existe un offset con una
fraccion coherente muchisimo mayor que el offset=0 actual, confirma un error de orientacion
sistematico y CORREGIBLE (rotar los datos ese angulo antes de medir). Si ningun offset mejora
mucho, es dispersion real, no una rotacion fija.

    python diagnose_rus_offset.py
"""
import pandas as pd
import numpy as np

df = pd.read_csv("t1_data/picks_polarization.csv")
is_assoc = df["associated"].astype(bool) if df["associated"].dtype != object else \
    df["associated"].astype(str).str.strip().str.upper().isin(["TRUE", "1"])
a = df[is_assoc].copy()
a = a[a["rectilinearity"] >= 0.60].copy()

for sta in sorted(a["sta"].unique()):
    sub = a[a["sta"] == sta]
    meas = sub["backazimuth_meas"].to_numpy(float)
    true = sub["ev_azimuth"].to_numpy(float) % 180.0
    print(f"\n{'='*60}\nestacion {sta}  (n={len(sub)})\n{'='*60}")
    best = (0.0, -1.0, np.nan)
    rows = []
    for off in range(0, 180, 5):
        m2 = (meas + off) % 180.0
        d = np.abs(m2 - true)
        d = np.minimum(d, 180.0 - d)
        coh = float((d <= 20.0).mean())
        med = float(np.median(d))
        rows.append((off, coh, med))
        if coh > best[1]:
            best = (off, coh, med)
    for off, coh, med in rows:
        marca = "  <-- mejor" if off == best[0] else ""
        if coh >= 0.4 or off == 0:
            print(f"   offset={off:3d}  fraccion_coherente={coh:.3f}  mediana_dev={med:5.1f}{marca}")
    print(f"\n   offset=0 (actual): fraccion={rows[0][1]:.3f}")
    print(f"   MEJOR offset={best[0]}: fraccion={best[1]:.3f}  mediana={best[2]:.1f}")
    if best[1] >= 0.60 and best[0] != 0:
        print(f"   -> ROTACION FIJA DETECTADA: aplicar +{best[0]} grados a la orientacion de "
              f"{sta} explicaria el desajuste. Revisa el inventario/instalacion de esta estacion.")
    elif best[1] < 0.60:
        print(f"   -> NINGUN offset fijo llega a 0.60: no es una rotacion simple, parece "
              f"dispersion real (o un problema mas complejo que una rotacion constante).")
