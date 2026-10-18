import json, itertools, sys
sys.path.insert(0,"/home/claude/work/r1_verif_indep")
from v3_endpoint import op, truth
pk=json.load(open("/home/claude/work/r1/run/picks_all.json"))
tgc={i for i,e in enumerate(truth) if 0<=e["t"]<840 and e["M"]>=0.5}
GG=[dict(min_picks=m, sigma=s) for m, s in itertools.product([3, 4, 5, 6, 8], [1.0, 2.0])]
PG=[dict(n_picks=n, tol=t) for n, t in itertools.product([3, 4, 6, 8, 10], [0.5, 1.0, 2.0])]
rows=[]
for fam,(key,grid) in enumerate((("pyocto_1d",PG),("gamma_1d",GG))):
    assert list(pk[key].keys())==[json.dumps(c) for c in grid], key
    for gi,c in enumerate(grid):
        v=[e for e in pk[key][json.dumps(c)] if 0<=e[3]<840]
        g,thr,nf,n=op(v,0,840,0.5,30,tgc)
        rows.append((-len(g),n,fam,gi,key,json.dumps(c),nf))
rows.sort()
for r in rows: print("%-9s %-32s CAL hits %2d/30  n %2d  false %d"%(r[4],r[5],-r[0],r[1],r[6]))
print("SELECTED",rows[0][4],rows[0][5])
sel=json.load(open("/home/claude/work/r1/run/hyb_eval.json"))["selection"]
print("matches run table:", [(a,b,round(c*30),n) for a,b,c,n in sel]==[(r[4],r[5],-r[0],r[1]) for r in rows])
# BASE list in hyb_eval equals picks_all gamma config restricted to EVAL?
d=json.load(open("/home/claude/work/r1/run/hyb_eval.json"))
b=[e[:5] for e in pk[rows[0][4]][rows[0][5]] if 840<=e[3]<1960]
print("BASE list identical:", b==d["BASE"])
ok=all(d["BESTPICK"][f"{k}|{c}"]==[e[:5] for e in v if 840<=e[3]<1960] for k in ("pyocto_1d","gamma_1d") for c,v in pk[k].items())
print("BESTPICK lists identical:", ok, len(d["BESTPICK"]))
# edge truncation: any associator events near/after 1960 or before 840?
for k in ("pyocto_1d","gamma_1d"):
    for c,v in pk[k].items():
        ts=sorted(e[3] for e in v)
        pass
allt=[e[3] for k in ("pyocto_1d","gamma_1d") for v in pk[k].values() for e in v]
print("associator origin range %.1f %.1f; n in [1960,2000) %d; n in [835,840) %d"%(min(allt),max(allt),sum(1960<=t for t in allt),sum(835<=t<840 for t in allt)))
