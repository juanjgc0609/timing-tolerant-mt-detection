"""C11 POST HOC (written after the evaluation analysis): paired contrasts between the candidates themselves (the plan
only fixed contrasts to S0).  -> out/c11_posthoc.txt"""
import c10_analyze as A
import c11_analyze as CA

lines = []


def P(s=""):
    print(s, flush=True)
    lines.append(s)


names = ["MT@S0", "MT@S0+T3", "MT@S1+T3", "MT@S2+T3", "MT@S3+T3", "MT@S2@40+T3", "MT@S2@80+T3", "MT@ORA+T3"]
pairs = [("MT@S3+T3", "MT@S2+T3"), ("MT@S2@40+T3", "MT@S2+T3"), ("MT@S2@80+T3", "MT@S2+T3"), ("MT@ORA+T3", "MT@S2+T3"),
         ("MT@S2+T3", "MT@S1+T3"), ("MT@ORA+T3", "MT@S2@80+T3")]
P("C11 POST HOC: paired differences between candidates, @0.1 / @0.25 [95%], B = 2000, regimes pooled")
for design in ("STD", "DENSE"):
    rows, nets = CA.load(design, "c11_eval")
    for rule, tol in (("LOC+T0.5", 0.5), ("LOC", None)):
        recs = CA.build(rows, names, tol)
        bt = A.Boot(rows, recs, lambda r: True, 2000)
        P(f"== {design}, {rule}")
        for a, b in pairs:
            p, c, sd = bt.diff(a, b)
            P(f"   {a:12s} - {b:10s} {CA.fmt(p[0], c[0])} / {CA.fmt(p[1], c[1])} (sd {sd[1]:.4f})")
open("out/c11_posthoc.txt", "w").write("\n".join(lines) + "\n")
