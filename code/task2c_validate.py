"""
task2c_validate.py - validaciones de los metodos de la Tarea 2C (task2c_methods.py), solo con la semilla 999.

  V1  identidad del diccionario: atomos explicitos en el dominio de los datos, a_k = G[:, k] (x) ricker(t - t_bi) /
      sqrt(E_R), y = D sqrt(E_R) / sd; A^T y contra BlockDict.correlate(MF / sd) y A^T A contra BlockDict.cross_gram
      (bloques al azar, incluidos pares con solape temporal).
  V2  BOMP recupera nodo, origen y tensor (a escala) de una fuente unica sin ruido: (a) sintetica en el modelo (la
      imagen exacta de un bloque en el dominio del filtro adaptado); (b) forma de onda de la fisica de
      radiation_operator (Ricker en el tiempo continuo, tiempos de arribo no enteros), filtrada con matched: se exige
      nodo y origen, el coseno entre u y m se reporta (el muestreo en la muestra predicha entera, sin busqueda de pico,
      distorsiona las amplitudes); (c) la misma forma de onda con los arribos en la muestra predicha: tensor exacto.
  V3  GL-FW: KKT (certificado <= 1 + tol fuera del soporte, = 1 en el soporte) y brecha de dualidad en todos los
      lambda resueltos, en ventanas reales de la semilla 999, resolviendo 2 lambda mas alla de la parada temprana
      (extra=2); el registro con parada temprana es igual al de esa trayectoria mas larga; el primer bloque que entra es el primero de BOMP y
      lambda_max^2 / 2 es el primer puntaje de BOMP; en la sintetica sin ruido el primer bloque es el verdadero.
  V4  MTSTACK = energia de la proyeccion ||Q(x)^T v||^2 / 2 (el estadistico de WAVE-MT) con v = MF / sd en los tiempos
      predichos exactos (sin pico ni tope), en todos los nodos y origenes.
  (El control de reproduccion de MT-FW+SC y QM-3C+SC contra task2a_parts se hace en la prueba de humo y en cada red.)

    OPENBLAS_CORETYPE=Prescott python task2c_validate.py
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import task2b_ablation as B  # noqa: E402
import task2c_methods as M  # noqa: E402
from task0_common import AP, AS, DT, NT, SNR, T as TT, VP, VS, ricker  # noqa: E402
from task1w_mech import CP, CS, rand_dc  # noqa: E402

T = B.T
X = T.X
SEED = 999
OUT = os.path.join(HERE, "task2c_validation.txt")


def mvec(Mt):
    return np.array([Mt[0, 0], Mt[1, 1], Mt[2, 2], Mt[0, 1], Mt[0, 2], Mt[1, 2]])


def setup(cell, regime="thin", j=0):
    ci, nst, nclut, STA, statics, cat, drawn, corr = B.setup_net(SEED, regime, cell, j)
    grid = T.default_grid()
    Gs = [X.radiation_operator(g, STA) for g in np.asarray(grid)]
    BD = M.BlockDict(STA, grid, T.DEFAULT_TAUS, corr, Gs=Gs)
    return STA, statics, nclut, corr, grid, Gs, BD


def v1(L, cell, rng):
    STA, statics, nclut, corr, grid, Gs, BD = setup(cell)
    rngs = T.window_rngs(SEED, cell, 0, "thin")
    for _w in range(2):
        srcs, taus, D, _p, tr, _a = T.generate_window(STA, statics, nclut, T.CLUT_AMP["thin"], rngs)
    MF, sd = T.matched(D)
    MFn = MF / sd
    y = D * np.sqrt(M.E_R) / sd
    # bloques: 12 al azar + 4 pares con el mismo nodo y origenes vecinos (solape P/P y P/S en el tiempo)
    blocks = [(int(rng.integers(BD.N)), int(rng.integers(BD.T))) for _ in range(12)]
    for _ in range(4):
        n, t = int(rng.integers(BD.N)), int(rng.integers(BD.T - 1))
        blocks += [(n, t), (n, t + 1)]
    ns = BD.ns
    A = []
    for n, t in blocks:
        cols = np.zeros((6, ns, 3, NT))
        for b in range(2):
            for i in range(ns):
                pulse = ricker(TT - BD.idx[n, t, b, i] * DT) / np.sqrt(M.E_R)
                g = BD.Gb[n, b, i]                                     # (3, 6)
                cols[:, i, :, :] += np.einsum("ck,s->kcs", g, pulse)
        A.append(cols.reshape(6, -1))
    A = np.concatenate(A, 0)                                              # (6 nb, ns*3*NT)
    Aty = A @ y.reshape(-1)
    h = BD.corr_blocks(MFn, blocks).ravel()
    e1 = float(np.max(np.abs(Aty - h)) / np.max(np.abs(h)))
    AtA = A @ A.T
    K = BD.gram(blocks)
    e2 = float(np.max(np.abs(AtA - K)) / np.max(np.abs(K)))
    offd = float(np.max(np.abs(K - np.kron(np.eye(len(blocks)), np.ones((6, 6))) * K)))
    # diagonal de bloques de la Gram (BlockDict.K) contra cross_gram
    e3 = max(float(np.max(np.abs(BD.K[n, t] - BD.gram([(n, t)])))) for n, t in blocks)
    print(".", flush=True)
    L.append(f"V1 {cell}: {len(blocks)} bloques; max|A^T y - h| / max|h| = {e1:.1e}; max|A^T A - K| / max|K| = {e2:.1e} "
             f"(max |K| fuera de la diagonal de bloques {offd:.2f}); max|K_bb - gram([b])| = {e3:.1e}; "
             f"rango minimo de K_bb {int(BD.rank.min())}/6; min |t_S - t_P| {BD.min_ps_lag} muestras "
             f"(AC_HALF {M.AC_HALF}), max |AC(t_S - t_P)| {BD.max_ps_ac:.1e}")
    return e1 < 1e-6 and e2 < 1e-6 and e3 < 1e-12 and BD.rank.min() == 6


def synth_block(BD, n, t, u):
    R = np.zeros((BD.ns, 3, NT))
    BD.synth_sub(R, n, t, -np.asarray(u, float))
    return R


def cosine(a, b):
    return float(abs(a @ b) / (np.linalg.norm(a) * np.linalg.norm(b)))


def v2_v3_synth(L, cell, rng):
    """(a) sintetica en el modelo; (b) forma de onda fisica sin ruido (tcorr = 0, fuente en un nodo, origen en la
    rejilla)."""
    STA, statics, nclut, corr, grid, Gs, BD = setup(cell)
    BD0 = M.BlockDict(STA, grid, T.DEFAULT_TAUS, None, Gs=Gs)
    ok = True
    ra, rb, rc = [], [], []
    for rep in range(4):
        n0, t0 = int(rng.integers(BD.N)), int(rng.integers(3, BD.T - 3))
        m0 = mvec(rand_dc(rng)) * 40.0
        # (a)
        MFn = synth_block(BD, n0, t0, m0)
        dg = {}
        rec = M.bomp_record(BD, MFn, [BD.garr[n0]], diag=dg)
        (n1, t1) = dg["blocks"][0]
        u = dg["U"][0]
        err = float(np.linalg.norm(u - m0) / np.linalg.norm(m0))
        s2 = rec[1][0] / rec[0][0] if len(rec) > 1 else 0.0
        dgl = {}
        M.glfw_record(BD, MFn, MFn, 1.0, [BD.garr[n0]], diag=dgl)
        first_gl = tuple(dgl["first_blocks"][0][:2])
        ra.append(((n1, t1) == (n0, t0), err, s2, first_gl == (n0, t0)))
        ok &= (n1, t1) == (n0, t0) and err < 1e-9 and s2 < 1e-12 and first_gl == (n0, t0)
        # (b) forma de onda fisica: desplazamiento G m (x) Ricker en el tiempo de arribo continuo, SNR como el generador
        src, tau = BD.garr[n0], float(BD.taus[t0])
        D = np.zeros((BD.ns, 3, NT))
        d = np.linalg.norm(STA - src, axis=1)
        G = Gs[n0].reshape(2, BD.ns, 3, 6)
        mu = m0 / 40.0
        for b, V in enumerate((VP, VS)):
            for i in range(BD.ns):
                D[i] += SNR * np.outer(G[b, i] @ mu, ricker(TT - (d[i] / V + tau)))
        MF, _sd = T.matched(D)
        dg = {}
        M.bomp_record(BD0, MF, [src], diag=dg)
        (n2, t2) = dg["blocks"][0]
        c = cosine(dg["U"][0], mu)
        dgl = {}
        M.glfw_record(BD0, MF, D, 1.0, [src], diag=dgl)
        rb.append(((n2, t2) == (n0, t0), c, tuple(dgl["first_blocks"][0][:2]) == (n0, t0)))
        ok &= (n2, t2) == (n0, t0)
        # (c) la misma forma de onda con los arribos en la muestra predicha (tiempos enteros): recuperacion exacta
        D = np.zeros((BD.ns, 3, NT))
        for b in range(2):
            for i in range(BD.ns):
                D[i] += SNR * np.outer(G[b, i] @ mu, ricker(TT - BD0.idx[n0, t0, b, i] * DT))
        MF, _sd = T.matched(D)
        dg = {}
        M.bomp_record(BD0, MF, [src], diag=dg)
        c3 = cosine(dg["U"][0], mu)
        rc.append((tuple(dg["blocks"][0]) == (n0, t0), c3))
        ok &= tuple(dg["blocks"][0]) == (n0, t0) and c3 > 1 - 1e-10
    print(L[-1] if L else "", flush=True)
    L.append(f"V2 {cell} (a) en el modelo, 4 fuentes: bloque recuperado {[r[0] for r in ra]}; |u - m| / |m| max "
             f"{max(r[1] for r in ra):.1e}; segundo puntaje / primero max {max(r[2] for r in ra):.1e}; primer bloque de GL-FW "
             f"= verdadero {[r[3] for r in ra]}")
    L.append(f"V2 {cell} (b) forma de onda fisica sin ruido: (nodo, origen) recuperado {[r[0] for r in rb]}; coseno(u, m) "
             f"min {min(r[1] for r in rb):.5f}; primer bloque de GL-FW = verdadero {[r[2] for r in rb]}")
    L.append(f"V2 {cell} (c) la misma forma de onda con arribos en la muestra predicha: bloque recuperado "
             f"{[r[0] for r in rc]}; 1 - coseno(u, m) max {max(1 - r[1] for r in rc):.1e}")
    return ok


def v3_windows(L, cell, regime, nwin):
    STA, statics, nclut, corr, grid, Gs, BD = setup(cell, regime)
    rngs = T.window_rngs(SEED, cell, 0, regime)
    ok = True
    kk_off, kk_on, gaps, gap_lw, same, first_eq, lam_eq, nlam, stalls, sup = [], [], [], [], [], [], [], [], 0, []
    t_full = t_early = 0.0
    for _w in range(nwin):
        srcs, taus, D, _p, tr, _a = T.generate_window(STA, statics, nclut, T.CLUT_AMP[regime], rngs)
        MF, sd = T.matched(D)
        MFn = MF / sd
        t0 = time.time()
        dgf = {}
        rf = M.glfw_record(BD, MFn, D, sd, tr, diag=dgf, extra=2)
        t_full += time.time() - t0
        print(f"   V3 {regime} {cell} ventana {_w}: {time.time() - t0:.1f} s, lambdas {len(dgf['per_lambda'])}, soporte final {dgf['per_lambda'][-1]['n_support']}", flush=True)
        t0 = time.time()
        dge = {}
        re_ = M.glfw_record(BD, MFn, D, sd, tr, diag=dge)
        t_early += time.time() - t0
        same.append(rf == re_)
        pl = dgf["per_lambda"]
        kk_off.append(max(p["kkt_off_max"] for p in pl))
        kk_on.append(max(p["kkt_on_maxdev"] for p in pl))
        gaps.append(max(p["rel_gap"] for p in pl if np.isfinite(p["rel_gap"])))
        gap_lw.append(max(p["gap_over_lam_omega"] for p in pl))
        stalls += sum(p["stall"] for p in pl)
        sup.append(pl[-1]["n_support"])
        nlam.append(len(pl))
        dg = {}
        rb = M.bomp_record(BD, MFn, tr, diag=dg)
        lam1 = max(lam for lam, _nv in dgf["entry"].values())
        first_eq.append(dgf["entry"].get(tuple(dg["blocks"][0]), (None,))[0] == lam1)
        lam_eq.append(abs(dgf["lam_max"] ** 2 / 2 - rb[0][0]) / rb[0][0])
    L.append(f"V3 {regime} {cell}: {nwin} ventanas, trayectoria hasta 2 lambda despues de completar el registro: max KKT fuera del soporte "
             f"{max(kk_off):.9f} (tol 1 + {M.TOL_FW:.0e}); max |c/lambda - 1| en el soporte {max(kk_on):.1e}; max brecha "
             f"relativa (brecha / (F(0) - F(u))) {max(gaps):.1e}; max brecha / (lambda Omega) {max(gap_lw):.1e}; re-solves "
             f"{stalls}; soporte en lambda_max/30: {min(sup)}-{max(sup)}; registro con parada temprana = con 2 lambda mas "
             f"{sum(same)}/{nwin}; el primer bloque de BOMP entra en el primer lambda con soporte no vacio {sum(first_eq)}/{nwin}; max |lambda_max^2/2 - "
             f"puntaje BOMP 1| rel {max(lam_eq):.1e}; tiempo por ventana con 2 lambda mas {t_full / nwin:.2f} s, con parada "
             f"{t_early / nwin:.2f} s")
    ok &= max(kk_off) <= 1 + 1e-6 and max(kk_on) <= 1e-6 and max(gaps) <= 1e-6 and all(same) and all(first_eq) \
        and max(lam_eq) < 1e-10 and stalls == 0
    return ok


def v4(L, cell, regime):
    STA, statics, nclut, corr, grid, Gs, BD = setup(cell, regime)
    Qs = [X.basis(G) for G in Gs]
    ST = M.Stacker(BD)
    rngs = T.window_rngs(SEED, cell, 0, regime)
    err, ok = 0.0, True
    for _w in range(3):
        srcs, taus, D, _p, tr, _a = T.generate_window(STA, statics, nclut, T.CLUT_AMP[regime], rngs)
        MF, sd = T.matched(D)
        MFn = MF / sd
        F = ST.field(MFn)
        st = np.arange(BD.ns)
        for gi, src in enumerate(BD.garr):
            idx = M._idx_node(STA, src, BD.taus, corr)                    # (nt, 2, ns)
            v = MFn[st[None, None, :], :, idx].reshape(len(BD.taus), -1)  # (nt, 6ns), orden (b, i, c)
            e = np.sum((v @ Qs[gi]) ** 2, axis=1) / 2.0
            err = max(err, float(np.max(np.abs(F[gi] - e) / np.maximum(e, 1e-300))))
        rec = ST.record(MFn, tr)
        ok &= abs(rec[0][0] - F.max()) <= 1e-12 * F.max()
    L.append(f"V4 {regime} {cell}: 3 ventanas x {BD.N} nodos x {BD.T} origenes: max |MTSTACK - ||Q^T v||^2/2| / "
             f"(||Q^T v||^2/2) = {err:.1e}; primer puntaje del registro = maximo del campo: {ok}")
    return ok and err < 1e-10


def main():
    t0 = time.time()
    L = [f"validacion de task2c_methods.py (semilla {SEED}); nucleo {B.blas_arch()} "
         f"(OPENBLAS_CORETYPE={os.environ.get('OPENBLAS_CORETYPE')})"]
    ok = True
    rng = np.random.default_rng([SEED, 2601])
    for cell in T.CELLS:
        ok &= v1(L, cell, rng)
    for cell in T.CELLS:
        ok &= v2_v3_synth(L, cell, rng)
    for regime in T.REGIMES:
        for cell in T.CELLS:
            ok &= v3_windows(L, cell, regime, 4)
    for regime in T.REGIMES:
        for cell in T.CELLS:
            ok &= v4(L, cell, regime)
    L.append(f"validacion {'CORRECTA' if ok else 'CON FALLAS'} ({time.time() - t0:.0f} s)")
    txt = "\n".join(L)
    print(txt)
    with open(OUT, "w") as fh:
        fh.write(txt + "\n")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
