"""C7 power calculation from the pilot bootstrap SDs (out/c7_bootsd_<design>.json, written by c7_analyze.py).
SE(n) = SD_pilot sqrt(n_pilot / n), n = networks per cell and regime (5 windows per network, as in the pilot);
one-sided alpha 0.025; power(n, delta) = Phi(delta / SE(n) - z_0.975); n(delta, power) = n_pilot (SD (z_0.975 + z_pow)
/ delta)^2.
    python c7_power.py > out/c7_power.txt"""
import json
import os

import numpy as np
from scipy.stats import norm

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
NP = 20
CONTRASTS = [("STD", "BOMP-PK-PH|BOMP-PK"), ("STD", "BOMP-PK-PHF|BOMP-PK"), ("STD", "BOMP-PK-PH-CAP|BOMP-PK-CAP"),
             ("DENSE", "BOMP-PK-PH|BOMP-PK"), ("DENSE", "BOMP-PK-PH-CAP|BOMP-PK-CAP"), ("DENSE", "MT-FW|BOMP-PK"),
             ("DENSE", "BOMP-PK-CAP|BOMP-PK"), ("GABOR90", "BOMP-PK-PH|BOMP-PK"), ("GABOR90", "BOMP-PK-PHF|BOMP-PK"),
             ("BERLAGE90", "BOMP-PK-PH|BOMP-PK"), ("BERLAGE90", "BOMP-PK-PHF|BOMP-PK")]
za = norm.ppf(0.975)
print(f"pilot: {NP} networks per cell and regime, 5 windows per network; one-sided alpha 0.025")
print("design     contrast                         pool    SD@0.1  SD@0.25   n(d=.05,pw=.95)@0.25  n(d=.03,pw=.90)@0.25  power@0.25 n=100 d=.05 / d=.03")
for des, key in CONTRASTS:
    p = os.path.join(OUT, f"c7_bootsd_{des}.json")
    if not os.path.exists(p):
        continue
    d = json.load(open(p))
    for pool in ("thin", "strong", "both"):
        if pool not in d or key not in d[pool]:
            continue
        sd = d[pool][key]
        npil = NP * (2 if pool == "both" else 1)
        def n_req(delta, pw):
            return int(np.ceil(npil * (sd[1] * (za + norm.ppf(pw)) / delta) ** 2))
        def pwr(n, delta):
            return norm.cdf(delta / (sd[1] * np.sqrt(npil / n)) - za)
        nn = 100 * (2 if pool == "both" else 1)
        print(f"{des:10s} {key:32s} {pool:6s}  {sd[0]:.4f}  {sd[1]:.4f}   {n_req(0.05, 0.95):6d}"
              f"                {n_req(0.03, 0.90):6d}              {pwr(nn, 0.05):.3f} / {pwr(nn, 0.03):.3f}")
print("(for pool 'both' n counts networks per cell summed over the two regimes)")
