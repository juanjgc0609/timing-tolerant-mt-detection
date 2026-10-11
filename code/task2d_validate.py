"""
task2d_validate.py - validaciones de los metodos de la Tarea 2D (task2d_methods.py), solo con la semilla 999.

  V1  peak_table: la muestra de pico tabulada es identica al tpk de node_vectors_c en todos los nodos, origenes,
      estaciones y fases (MFn y un residuo deflactado por un paso de MT-FW+SC).
  V2  reducciones exactas: con la lectura en el pico desactivada (peak=False), BOMP-PK, MTSTACK-PK y GL-FW-PK dan
      EXACTAMENTE los registros de task2c_methods (BOMP+SC, MTSTACK+SC, GL-FW+SC); fw_record_v(peak, cap) da
      exactamente certified_fw_record (MT-FW+SC).  Ventanas reales 0-1 de la red 0 de 2A de la semilla 999.
  V3  identidades del diccionario con tiempos de lectura desplazados: atomos explicitos en el dominio de los datos
      a_k = G[:, k] (x) ricker(t - I_bi DT) / sqrt(E_R) con I = muestra predicha + desplazamiento al azar en [-4, 4]:
      A^T y = PKDict.corr_atoms, A^T A = PKDict.gram, K_atom = gram([a]), y la seudoinversa por subconjuntos de
      PKDict.energy igual a la de K_atom.
  V4  estadisticos por caminos independientes: primer puntaje de MT-FW-NOCAP = maximo del campo de MTSTACK-PK (proyeccion
      por pinv(G) frente a la base Q); primer puntaje de MT-FW-EX = maximo de ||Q^T vt||^2 / 2 con vt leido en la muestra
      exacta y con tope (vectorizado); primer puntaje de BOMP-PK = lambda_max^2 / 2 de GL-FW-PK; ortogonalidad del
      reajuste de BOMP-PK en sus tiempos fijos.
  V5  GL-FW-PK en ventanas reales: KKT del oraculo de pico (max c / lambda <= 1 + 1e-6), |c / lambda - 1| en el soporte,
      brecha de dualidad relativa en cada lambda resuelto, registro con parada temprana = con 2 lambda mas.
  V6  error de tiempo inyectado: generador ideal de 1x (task1x_mtglrt.generate_x, sin perturbaciones, clutter 0.35),
      fuentes movidas al nodo mas cercano y origen a la rejilla de 0.1 s (la E3 de task2c_explore.py: solo queda el
      redondeo a la muestra), 30 ventanas por celda.  Los detectores se corren con tcorr = 0 (sin error) y con
      tcorr = +-0.2 s de signo al azar por estacion y fase (error de tiempo de 0.2 s = 3.2 muestras, dentro de la
      ventana de +-0.25 s).  Se espera que los metodos de tiempo exacto pierdan recall con el error y que sus variantes
      PK lo recuperen.

    OPENBLAS_CORETYPE=Prescott python task2d_validate.py [--workers 2]
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import task2b_ablation as B  # noqa: E402
import task2c_methods as M  # noqa: E402
import task2d_methods as D2  # noqa: E402
from certified_fw import node_vectors_c  # noqa: E402
from froc_exact import exact_HF  # noqa: E402
from task0_common import DEFAULT_TAUS, DT, NT, T as TT, default_grid, make_network, ricker  # noqa: E402
from task1o_consistent import draw_sources, matched  # noqa: E402

T = B.T
X = T.X
SEED = 999
OUT = os.path.join(HERE, "task2d_validation.txt")
N_V6 = 30
ERR_S = 0.2


def setup(cell, regime="thin", j=0):
    ci, nst, nclut, STA, statics, cat, drawn, corr = B.setup_net(SEED, regime, cell, j)
    grid = T.default_grid()
    Gs = [X.radiation_operator(g, STA) for g in np.asarray(grid)]
    Qs = [X.basis(G) for G in Gs]
    BD = M.BlockDict(STA, grid, T.DEFAULT_TAUS, corr, Gs=Gs)
    return STA, statics, nclut, corr, grid, Gs, Qs, BD, M.Stacker(BD), D2.PKDict(BD)


def as_rec(r):
    return [[float(s), int(h)] for s, h in r]


def recall(recs, k):
    _t, H, F = exact_HF([[tuple(x) for x in r] for r in recs])
    H, F = np.asarray(H, float), np.asarray(F, float)
    out = []
    for cap in (0.1, 0.25):
        ok = F.sum(0) / len(recs) <= cap + 1e-12
        out.append(float(H.sum(0)[ok].max() / max(1, sum(k))) if ok.any() else 0.0)
    return out


# ------------------------------------------------------------------------------------------------ V1-V5
def v1_to_v5(args):
    cell, regime = args
    L = []
    ok = True
    STA, statics, nclut, corr, grid, Gs, Qs, BD, ST, PD = setup(cell, regime)
    garr = BD.garr
    rngs = T.window_rngs(SEED, cell, 0, regime)
    rng = np.random.default_rng([SEED, 2711, T.CELLS.index(cell), T.REGIMES.index(regime)])
    v1_bad = v1_n = 0
    v2 = {"BOMP": 0, "MTSTACK": 0, "GL-FW": 0, "MT-FW": 0}
    v3 = []
    v4 = {"nocap_vs_stackpk": 0.0, "ex_vs_vec": 0.0, "bomp_vs_lam": 0.0, "orth": 0.0}
    v5 = {"kkt_off": 0.0, "kkt_on": 0.0, "gap": 0.0, "same": 0, "fw": [], "capped": 0, "stall": 0}
    nwin = 2
    for w in range(nwin):
        srcs, taus, D, _p, tr, _a = T.generate_window(STA, statics, nclut, T.CLUT_AMP[regime], rngs)
        MF, sd = T.matched(D)
        MFn = MF / sd
        # V1: MFn y un residuo deflactado (un paso de MT-FW+SC)
        R1 = MFn.copy()
        r1 = T.certified_fw_record(STA, MFn, grid, T.DEFAULT_TAUS, Gs, Qs, tr, kmax=1, tcorr=corr)
        gi = None
        best = -np.inf
        for g, src in enumerate(garr):
            vv, tp, _s, _e = node_vectors_c(MFn, STA, src, T.DEFAULT_TAUS, corr)
            sc = np.sum((vv.reshape(len(T.DEFAULT_TAUS), -1) @ Qs[g]) ** 2, axis=1) / 2.0
            if sc.max() > best:
                best, gi, jb, vb, tb = float(sc.max()), g, int(sc.argmax()), vv, tp
        assert abs(best - r1[0][0]) == 0.0
        m, *_ = np.linalg.lstsq(Gs[gi], vb[jb].reshape(-1), rcond=None)
        pred = (Gs[gi] @ m).reshape(2, BD.ns, 3)
        for b in range(2):
            for i in range(BD.ns):
                t0 = int(tb[jb, b, i])
                lo, hi = max(0, t0 - X.AC_HALF), min(NT, t0 + X.AC_HALF + 1)
                R1[i, :, lo:hi] -= np.outer(pred[b, i], X.AC[X.AC_HALF - (t0 - lo):X.AC_HALF + (hi - t0)])
        for R in (MFn, R1):
            I = PD.read_idx(R, True)
            for g, src in enumerate(garr):
                _v, tp, _s, _e = node_vectors_c(R, STA, src, T.DEFAULT_TAUS, corr)
                v1_bad += int(np.sum(tp != I[g]))
                v1_n += tp.size
        # V2
        v2["BOMP"] += as_rec(D2.bomp_pk_record(PD, MFn, tr, peak=False)) == as_rec(M.bomp_record(BD, MFn, tr))
        v2["MTSTACK"] += as_rec(D2.mtstack_pk_record(PD, ST, MFn, tr, peak=False)) == as_rec(ST.record(MFn, tr))
        v2["GL-FW"] += as_rec(D2.glfw_pk_record(PD, MFn, D, sd, tr, peak=False)) == as_rec(M.glfw_record(BD, MFn, D, sd, tr))
        rc = as_rec(T.certified_fw_record(STA, MFn, grid, T.DEFAULT_TAUS, Gs, Qs, tr, tcorr=corr))
        v2["MT-FW"] += as_rec(D2.fw_record_v(STA, MFn, grid, T.DEFAULT_TAUS, Gs, Qs, tr, corr, True, True)) == rc
        # V3: atomos explicitos con tiempos desplazados
        y = D * np.sqrt(M.E_R) / sd
        atoms = []
        for _ in range(6):
            n, t = int(rng.integers(BD.N)), int(rng.integers(BD.T))
            I = np.clip(BD.idx[n, t] + rng.integers(-4, 5, size=(2, BD.ns)), 0, NT - 1)
            atoms.append((n, t, I))
        n, t = int(rng.integers(BD.N)), int(rng.integers(BD.T - 1))
        I = np.clip(BD.idx[n, t] + rng.integers(-4, 5, size=(2, BD.ns)), 0, NT - 1)
        atoms += [(n, t, I), (n, t + 1, np.clip(I + rng.integers(-2, 3, size=I.shape), 0, NT - 1))]
        A = []
        for n, t, I in atoms:
            cols = np.zeros((6, BD.ns, 3, NT))
            for b in range(2):
                for i in range(BD.ns):
                    pulse = ricker(TT - I[b, i] * DT) / np.sqrt(M.E_R)
                    cols[:, i, :, :] += np.einsum("ck,s->kcs", BD.Gb[n, b, i], pulse)
            A.append(cols.reshape(6, -1))
        A = np.concatenate(A, 0)
        h = PD.corr_atoms(MFn, atoms).ravel()
        e1 = float(np.max(np.abs(A @ y.reshape(-1) - h)) / np.max(np.abs(h)))
        K = PD.gram(atoms)
        e2 = float(np.max(np.abs(A @ A.T - K)) / np.max(np.abs(K)))
        e3 = max(float(np.max(np.abs(PD.K_atom(a[0], a[2]) - PD.gram([a])))) for a in atoms)
        e_all, Ipk, hpk = PD.energy(MFn, True)
        e4 = 0.0
        for _ in range(20):
            n, t = int(rng.integers(BD.N)), int(rng.integers(BD.T))
            Ka = PD.K_atom(n, Ipk[n, t])
            ee = float(hpk[n, t] @ np.linalg.pinv(Ka) @ hpk[n, t])
            e4 = max(e4, abs(ee - e_all[n, t]) / max(1e-300, abs(ee)))
        v3.append((e1, e2, e3, e4))
        # V4
        rn = D2.fw_record_v(STA, MFn, grid, T.DEFAULT_TAUS, Gs, Qs, tr, corr, True, False)
        Ipk0 = PD.read_idx(MFn, True)
        vv = PD.sample_idx(MFn, Ipk0).reshape(BD.N, BD.T, -1)
        mm = np.matmul(vv, np.swapaxes(ST.Pinv, 1, 2))
        Fpk = np.sum(np.matmul(mm, np.swapaxes(BD.Gf, 1, 2)) ** 2, axis=-1) / 2.0
        v4["nocap_vs_stackpk"] = max(v4["nocap_vs_stackpk"], abs(rn[0][0] - Fpk.max()) / Fpk.max())
        rx = D2.fw_record_v(STA, MFn, grid, T.DEFAULT_TAUS, Gs, Qs, tr, corr, False, True)
        ve = PD.sample_idx(MFn, BD.idx)                                          # (N, T, 2, ns, 3)
        nv = np.linalg.norm(ve, axis=-1, keepdims=True)
        ve = (ve * np.minimum(1.0, X.CAPV / np.maximum(nv, 1e-12))).reshape(BD.N, BD.T, -1)
        Qst = np.stack([np.pad(Q, ((0, 0), (0, 6 - Q.shape[1]))) for Q in Qs])
        Fex = np.sum(np.einsum("ntr,nrk->ntk", ve, Qst) ** 2, axis=-1) / 2.0
        v4["ex_vs_vec"] = max(v4["ex_vs_vec"], abs(rx[0][0] - Fex.max()) / Fex.max())
        dgb = {}
        rb = D2.bomp_pk_record(PD, MFn, tr, diag=dgb)
        v4["orth"] = max(v4["orth"], max(dgb["orth"]))
        # V5
        dgf, dge = {}, {}
        rf = D2.glfw_pk_record(PD, MFn, D, sd, tr, diag=dgf, extra=2)
        re_ = D2.glfw_pk_record(PD, MFn, D, sd, tr, diag=dge)
        v4["bomp_vs_lam"] = max(v4["bomp_vs_lam"], abs(dgf["lam_max"] ** 2 / 2 - rb[0][0]) / rb[0][0])
        pl = dgf["per_lambda"]
        v5["kkt_off"] = max(v5["kkt_off"], max(p["kkt_off_max"] for p in pl))
        v5["kkt_on"] = max(v5["kkt_on"], max(p["kkt_on_maxdev"] for p in pl))
        v5["gap"] = max(v5["gap"], max(p["rel_gap"] for p in pl if np.isfinite(p["rel_gap"])))
        v5["same"] += rf == re_
        v5["fw"].append(sum(p["n_fw"] for p in pl))
        v5["capped"] += sum(p["fw_capped"] for p in pl)
        v5["stall"] += sum(p["stall"] for p in pl)
    e = np.array(v3).max(axis=0)
    L.append(f"V1 {regime} {cell}: muestras de pico tabuladas distintas del tpk de node_vectors_c: {v1_bad} de {v1_n} "
             f"(MFn y residuo tras un paso de MT-FW+SC, {nwin} ventanas)")
    L.append(f"V2 {regime} {cell}: registros identicos con peak=False en {nwin} ventanas: BOMP {v2['BOMP']}, MTSTACK "
             f"{v2['MTSTACK']}, GL-FW {v2['GL-FW']}; fw_record_v(pico, tope) = certified_fw_record: {v2['MT-FW']}")
    L.append(f"V3 {regime} {cell}: 8 atomos con tiempos desplazados por ventana: max|A^T y - h|/max|h| {e[0]:.1e}; "
             f"max|A^T A - K|/max|K| {e[1]:.1e}; max|K_atom - gram([a])| {e[2]:.1e}; energia por subconjuntos contra "
             f"pinv(K_atom) rel {e[3]:.1e}")
    L.append(f"V4 {regime} {cell}: |NOCAP primer puntaje - max MTSTACK-PK| rel {v4['nocap_vs_stackpk']:.1e}; |EX primer "
             f"puntaje - max vectorizado| rel {v4['ex_vs_vec']:.1e}; |BOMP-PK primer puntaje - lambda_max^2/2 de "
             f"GL-FW-PK| rel {v4['bomp_vs_lam']:.1e}; ortogonalidad del reajuste de BOMP-PK max {v4['orth']:.1e}")
    L.append(f"V5 {regime} {cell}: GL-FW-PK con 2 lambda mas: max c/lambda del oraculo {v5['kkt_off']:.9f}; max "
             f"|c/lambda - 1| en el soporte {v5['kkt_on']:.1e}; max brecha relativa {v5['gap']:.1e}; registro con parada "
             f"= con 2 lambda mas {v5['same']}/{nwin}; pasos FW por ventana {v5['fw']}; lambdas con tope de pasos FW "
             f"{v5['capped']}; re-solves {v5['stall']}")
    ok &= v1_bad == 0 and all(v == nwin for v in v2.values()) and e[0] < 1e-6 and e[1] < 1e-6 and e[2] < 1e-12 \
        and e[3] < 1e-8 and v4["nocap_vs_stackpk"] < 1e-9 and v4["ex_vs_vec"] < 1e-9 and v4["bomp_vs_lam"] < 1e-9 \
        and v4["orth"] < 1e-9 and v5["kkt_off"] <= 1 + 1e-6 and v5["kkt_on"] <= 1e-6 and v5["gap"] <= 1e-6 \
        and v5["same"] == nwin and v5["capped"] == 0
    return ok, L


# ----------------------------------------------------------------------------------------------------- V6
V6_METHODS = ("MT-FW", "MT-FW-EX", "MT-FW-NOCAP", "BOMP", "BOMP-PK", "MTSTACK", "MTSTACK-PK", "GL-FW", "GL-FW-PK")


def v6_cell(args):
    cell, cond = args
    grid = default_grid()
    garr = np.asarray(grid)
    ci = T.CELLS.index(cell)
    kind, nst, nclut = T.SCEN[cell]
    rng = np.random.default_rng([SEED, 2710, ci])
    rng_m = np.random.default_rng([SEED, 2711, ci])
    rng_h = np.random.default_rng([SEED, 2712, ci])
    rng_e = np.random.default_rng([SEED, 2714, ci])
    recs = {m: [] for m in V6_METHODS}
    ks = []
    t0 = time.time()
    for w in range(N_V6):
        S = make_network(kind, nst, seed=[SEED, 2713, ci, w])
        srcs, taus = draw_sources(rng)
        srcs = [garr[int(np.argmin(np.linalg.norm(garr - s_, axis=1)))].copy() for s_ in srcs]
        taus = [float(DEFAULT_TAUS[int(np.argmin(np.abs(DEFAULT_TAUS - t_)))]) for t_ in taus]
        D, _p, tr, _rad = X.generate_x(S, srcs, taus, nclut, rng, rng_m, False, rng_h)
        err = ERR_S * rng_e.choice([-1.0, 1.0], size=(nst, 2))
        tc = err if cond == "error" else np.zeros((nst, 2))
        MF, sd = matched(D)
        MFn = MF / sd
        Gs = [X.radiation_operator(g, S) for g in garr]
        Qs = [X.basis(G) for G in Gs]
        BD = M.BlockDict(S, grid, DEFAULT_TAUS, tc, Gs=Gs)
        ST = M.Stacker(BD)
        PD = D2.PKDict(BD)
        recs["MT-FW"].append(T.certified_fw_record(S, MFn, grid, DEFAULT_TAUS, Gs, Qs, tr, tcorr=tc))
        recs["MT-FW-EX"].append(D2.fw_record_v(S, MFn, grid, DEFAULT_TAUS, Gs, Qs, tr, tc, False, True))
        recs["MT-FW-NOCAP"].append(D2.fw_record_v(S, MFn, grid, DEFAULT_TAUS, Gs, Qs, tr, tc, True, False))
        recs["BOMP"].append(M.bomp_record(BD, MFn, tr))
        recs["BOMP-PK"].append(D2.bomp_pk_record(PD, MFn, tr))
        recs["MTSTACK"].append(ST.record(MFn, tr))
        recs["MTSTACK-PK"].append(D2.mtstack_pk_record(PD, ST, MFn, tr))
        recs["GL-FW"].append(M.glfw_record(BD, MFn, D, sd, tr))
        recs["GL-FW-PK"].append(D2.glfw_pk_record(PD, MFn, D, sd, tr))
        ks.append(len(srcs))
    return cell, cond, {m: recall(r, ks) for m, r in recs.items()}, sum(ks), time.time() - t0


def main(workers=2):
    t0 = time.time()
    L = [f"validacion de task2d_methods.py (semilla {SEED}); nucleo {B.blas_arch()} "
         f"(OPENBLAS_CORETYPE={os.environ.get('OPENBLAS_CORETYPE')})"]
    ok = True
    jobs = [(cell, reg) for reg in T.REGIMES for cell in T.CELLS]
    jobs6 = [(cell, cond) for cell in T.CELLS for cond in ("sin error", "error")]
    with ProcessPoolExecutor(workers) as ex:
        res = list(ex.map(v1_to_v5, jobs))
        res6 = list(ex.map(v6_cell, jobs6))
    for o, l_ in res:
        ok &= o
        L += l_
    R6 = {(c, cd): (r, n, s) for c, cd, r, n, s in res6}
    L.append(f"V6 generador ideal, fuentes en nodos y origen en la rejilla, {N_V6} ventanas por celda; error de tiempo "
             f"inyectado +-{ERR_S} s por estacion y fase (en tcorr de los detectores). Recall a 0.1 / 0.25:")
    ok6 = True
    for cell in T.CELLS:
        for cond in ("sin error", "error"):
            r, n, s = R6[(cell, cond)]
            L.append(f"   {cell:15s} {cond:9s} ({n} eventos, {s:.0f} s): " + " | ".join(
                f"{m} {r[m][0]:.2f}/{r[m][1]:.2f}" for m in V6_METHODS))
        r0, re_ = R6[(cell, "sin error")][0], R6[(cell, "error")][0]
        lost = {m: r0[m][1] - re_[m][1] for m in V6_METHODS}
        rec_ = {pk: re_[pk][1] - re_[ex][1] for ex, pk in (("MT-FW-EX", "MT-FW"), ("BOMP", "BOMP-PK"),
                                                            ("MTSTACK", "MTSTACK-PK"), ("GL-FW", "GL-FW-PK"))}
        L.append(f"   {cell:15s} perdida por el error a 0.25 (sin error - con error): " + ", ".join(
            f"{m} {v:+.2f}" for m, v in lost.items()) + " ; ganancia de la lectura en el pico con error (PK - exacto): "
            + ", ".join(f"{m} {v:+.2f}" for m, v in rec_.items()))
        ok6 &= all(v >= 0.0 for v in rec_.values())
    L.append(f"V6: la lectura en el pico recupera (PK >= exacto con error a 0.25) en todas las celdas y metodos: "
             f"{'si' if ok6 else 'NO'}")
    L.append(f"validacion {'CORRECTA' if ok else 'CON FALLAS'} (V1-V5; V6 es descriptiva) ({time.time() - t0:.0f} s)")
    txt = "\n".join(L)
    print(txt)
    with open(OUT, "w") as fh:
        fh.write(txt + "\n")
    return ok


if __name__ == "__main__":
    w = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 2
    raise SystemExit(0 if main(w) else 1)
