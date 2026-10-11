"""task2e_explore.py - EXPLORATORIO (escrito despues de ver los resultados de 2E; no decide ninguna regla).
Para TS-SFW+SC y BOMP-PK-SLIDE+SC en las 4800 ventanas de evaluacion: desplazamiento epicentral de cada atomo desde el
nodo de la rejilla en que entro hasta su posicion final, y cambio del estado de acierto (5 km) entre el nodo de entrada
y la posicion final; tiempos reales por ventana y metodo.  Salida: task2e_explore.txt."""
import glob
import json
import os

import numpy as np

import task2e_methods as E
from task2c_methods import hit_of
from task0_common import DEFAULT_TAUS, default_grid

HERE = os.path.dirname(os.path.abspath(__file__))
garr = np.asarray(default_grid())
L = ["# EXPLORATORIO (posterior a los resultados): deslizamiento de los atomos y tiempos reales de 2E", ""]
for reg in ("thin", "strong"):
    for m in ("TS-SFW+SC", "BOMP-PK-SLIDE+SC"):
        disp, flip_in, flip_out, n_at, n_first_disp = [], 0, 0, 0, []
        for p in sorted(glob.glob(os.path.join(HERE, "task2e_parts", f"ev_seed*_{reg}_*.json"))):
            e = json.load(open(p))
            for w in e["windows"]:
                d = w["diag"][m]
                for k, (it, pos) in enumerate(zip(d["iters"], d["pos"])):
                    node = garr[it[1]]
                    dd = float(np.linalg.norm(np.asarray(pos[:2]) - node[:2]))
                    disp.append(dd)
                    h0 = hit_of(node, w["truth"])
                    h1 = hit_of(np.asarray(pos[:3]), w["truth"])
                    flip_in += (h0 < 0 and h1 >= 0)
                    flip_out += (h0 >= 0 and h1 < 0)
                    n_at += 1
        disp = np.array(disp)
        L.append(f"- {reg} {m}: {n_at} atomos; desplazamiento epicentral desde el nodo de entrada mediana {np.median(disp):.2f} km, "
                 f"P90 {np.percentile(disp, 90):.2f} km, > 1.5 km en {np.mean(disp > 1.5):.3f}, > 5 km en {np.mean(disp > 5):.3f}; "
                 f"estado de acierto: falso en el nodo -> acierto al final {flip_in}, acierto en el nodo -> falso al final {flip_out}")
L.append("")
tim = {}
for p in glob.glob(os.path.join(HERE, "task2e_parts", "ev_seed*.json")):
    e = json.load(open(p))
    for m, v in e["timing"].items():
        tim.setdefault(m, []).extend(v)
L.append("Tiempo real por ventana (s, media / P90) en la evaluacion: " + "; ".join(
    f"{m} {np.mean(v):.2f} / {np.percentile(v, 90):.2f}" for m, v in sorted(tim.items())))
txt = "\n".join(L)
print(txt)
open(os.path.join(HERE, "task2e_explore.txt"), "w").write(txt + "\n")
