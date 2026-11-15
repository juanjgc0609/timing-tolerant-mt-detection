import sys, numpy as np, multiprocessing as mp
sys.argv=["x"]
import importlib.util
spec=importlib.util.spec_from_file_location("cal","" + __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "run.py") + "")
C=importlib.util.module_from_spec(spec); spec.loader.exec_module(C)
def f(c):
    srcs,Y,q=C._scene((C.SEED,2,2,1,c),"C3",True)
    a=C.m6.adcg_loop(C.E,Y,q,4,C.LAG,thr=48.5,free="joint")
    out=[]
    for x in C.MB.spurious(a,srcs):
        d=[ (round(float(np.linalg.norm(x["x"][:2]-s[0][:2])),1), round(x["tau"]-s[1],2)) for s in srcs]
        dn=min(float(np.linalg.norm(x["x"][:2]-b["x"][:2])) for b in a if b is not x)
        out.append((c, x["pair_gain"], d, round(dn,1)))
    return out
with mp.get_context("fork").Pool(4) as p:
    for o in p.map(f, range(200)):
        for r in o: print(r)
