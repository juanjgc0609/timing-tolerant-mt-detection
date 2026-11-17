"""Machine check (descriptive, no ledger): frozen configurations without terms, re-scored in this cloud container, so
the paired a6 differences are not a Windows-vs-Linux effect."""
import json, sys, time
sys.path.insert(0, "experiments/exploratory/COMPARATORS")
from ttmt.associate import pipeline as PL, tune as TU
out = {}
for m, f in (("phasenet+pyocto", "phasenet+pyocto.json"), ("phasenet+gamma", "phasenet+gamma_a4.json")):
    fr = json.load(open("experiments/exploratory/COMPARATORS/frozen/" + f))
    recs = [r for p in ("/mnt/user-data/uploads/tune/ss_ridgecrest_2019_egf_v1/tune_phasenet.json",
                        "/mnt/user-data/uploads/tune/ss_hinet_2012_egf_v1/tune_phasenet.json")
            for r in TU.load_manifest(p, m)]
    run = PL.runner(m); t0 = time.time()
    o = TU.objective(recs, [run(fr["params"], r) for r in recs])
    out[m] = dict(frozen_ledger=[fr["objective"], fr["recall_T"], fr["recall_C"]], cloud=[round(x, 6) for x in o[:3]],
                  wall_s=round(time.time() - t0, 1))
    print(m, out[m], flush=True)
json.dump(out, open("results/st_retune/machine_check.json", "w"), indent=1)
