"""C5 analysis: location errors (NODE, TS, CR) of hits and coverage of the U1 (AA.6 curvature) and U3 (travel-time
residual) uncertainty regions; inflation factor fitted on networks of split A, coverage checked on split B."""
import json, sys
import numpy as np

Q68_2, Q95_2 = 2.2789, 5.9915
Q68_1, Q95_1 = 0.9945, 1.9600


def stats(e):
    e = np.asarray(e)
    return f"median {np.median(e):.2f} P90 {np.percentile(e, 90):.2f} (n {len(e)})"


def cov_ok(C):
    if C is None:
        return False
    C = np.asarray(C)
    return np.all(np.isfinite(C)) and np.all(np.linalg.eigvalsh(C[:2, :2]) > 0) and C[2, 2] > 0


def main(path):
    rows = [r for r in json.load(open(path)) if r.get("hit")]
    allr = json.load(open(path))
    out = [f"windows {len(allr)}, BOMP-PK hits {len(rows)}"]
    tr = np.array([r["truth"] for r in rows])
    node = np.array([r["node"][:3] for r in rows])
    cr = np.array([r["cr"][:3] for r in rows])
    ts = [r for r in rows if r["ts"] is not None]
    tsl = np.array([r["ts"][:3] for r in ts])
    trs = np.array([r["truth"] for r in ts])
    for nm, L, TT in (("NODE", node, tr), ("CR", cr, tr), ("TS", tsl, trs)):
        eh = np.linalg.norm(L[:, :2] - TT[:, :2], axis=1)
        ez = L[:, 2] - TT[:, 2]
        out.append(f"  {nm:5s} epicentral {stats(eh)}; depth |err| {stats(np.abs(ez))}, bias {np.median(ez):+.2f} km")
    # paired improvement CR vs NODE
    d = np.linalg.norm(cr[:, :2] - tr[:, :2], axis=1) - np.linalg.norm(node[:, :2] - tr[:, :2], axis=1)
    out.append(f"  CR - NODE epicentral error: median {np.median(d):+.2f} km, CR better in {np.mean(d < 0):.2f}")
    # coverage
    for U in ("U1", "U3"):
        m2, zz, nets = [], [], []
        for r in rows:
            if U == "U1":
                try:
                    C = np.linalg.inv(np.asarray(r["H"]))
                except np.linalg.LinAlgError:
                    continue
            else:
                C = r["C_tt"]
            if not cov_ok(C):
                continue
            C = np.asarray(C)
            e = np.asarray(r["cr"][:3]) - np.asarray(r["truth"])
            m2.append(float(e[:2] @ np.linalg.solve(C[:2, :2], e[:2])))
            zz.append(abs(e[2]) / np.sqrt(C[2, 2]))
            nets.append((r["cell"], r["net"]))
        m2, zz = np.array(m2), np.array(zz)
        out.append(f"  {U}: n {len(m2)}; nominal coverage epicentral 68/95: {np.mean(m2 <= Q68_2):.2f}/{np.mean(m2 <= Q95_2):.2f}; "
                   f"depth 68/95: {np.mean(zz <= Q68_1):.2f}/{np.mean(zz <= Q95_1):.2f}; median sqrt(m2) {np.median(np.sqrt(m2)):.2f}")
        # inflation fitted on even networks (split A), tested on odd (split B)
        A = np.array([n[1] % 2 == 0 for n in nets])
        if A.any() and (~A).any():
            kh = np.sqrt(np.quantile(m2[A], 0.95) / Q95_2)
            kz = np.quantile(zz[A], 0.95) / Q95_1
            kh68 = np.sqrt(np.quantile(m2[A], 0.68) / Q68_2)
            out.append(f"     inflation fitted on split A: horizontal x{kh:.2f} (x{kh68:.2f} from the 68% quantile), depth x{kz:.2f}; "
                       f"split B coverage 68/95 epicentral {np.mean(m2[~A] <= Q68_2 * kh ** 2):.2f}/{np.mean(m2[~A] <= Q95_2 * kh ** 2):.2f}, "
                       f"depth {np.mean(zz[~A] <= Q68_1 * kz):.2f}/{np.mean(zz[~A] <= Q95_1 * kz):.2f} (n_B {int((~A).sum())})")
            # sharpness: median semi-axis of inflated 95% epicentral ellipse
    sd = []
    for r in rows:
        try:
            C = np.linalg.inv(np.asarray(r["H"]))
            if cov_ok(C):
                sd.append(np.sqrt(np.linalg.eigvalsh(C[:2, :2])).max())
        except np.linalg.LinAlgError:
            pass
    out.append(f"  U1 major semi-axis sd (km): median {np.median(sd):.3f}")
    sd = [np.sqrt(np.linalg.eigvalsh(np.asarray(r["C_tt"])[:2, :2])).max() for r in rows if cov_ok(r["C_tt"])]
    out.append(f"  U3 major semi-axis sd (km): median {np.median(sd):.3f}")
    return "\n".join(out)


if __name__ == "__main__":
    print(main(sys.argv[1]))
