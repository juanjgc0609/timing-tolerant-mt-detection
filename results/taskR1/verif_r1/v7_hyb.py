"""Own pick terms and own PICKSCORE-UNION from picks_all.json, p1_<seg>.json and the BASE list (own travel times)."""
import json, sys, math
import numpy as np
sys.path.insert(0,"/home/claude/work/r1_verif_indep")
from v6_tt import tt
from v3_endpoint import truth, lat0, lon0, op
SEG=sys.argv[1] if len(sys.argv)>1 else "eval"
RUN="/home/claude/work/r1/run"
pk=json.load(open(RUN+"/picks_all.json"))
from obspy import read
D="/home/claude/realdata/MatchLocate2/example/Trace/20120902/"
names=pk["names"]
KX=111.19*math.cos(math.radians(lat0))
STA=[]
for n in names:
    h=read(D+f"N.{n}.U",headonly=True)[0].stats.sac; STA.append(((h.stlo-lon0)*KX,(h.stla-lat0)*111.19))
STA=np.array(STA)
PA=[{b:np.array([p[0] for p in pl if p[1]==b]) for b in "PS"} for pl in pk["picks"]]
def T(src):
    X=np.hypot(STA[:,0]-src[0],STA[:,1]-src[1]); return np.stack([tt(src[2],X,1.0),tt(src[2],X,1.73)])
# pick terms over CAL truth events with M >= 1.0
res=[[[],[]] for _ in names]
for e in truth:
    if not (0<=e["t"]<840 and e["M"]>=1.0): continue
    t2=T((e["x"],e["y"],e["dep"]))
    for i in range(len(names)):
        for b,ph in enumerate("PS"):
            d=PA[i][ph]-(e["t"]+t2[b,i]); d=d[np.abs(d)<=2.0]
            if len(d): res[i][b].append(float(d[np.argmin(np.abs(d))]))
terms=np.array([[np.median(x) if len(x)>=2 else 0.0 for x in r] for r in res])
nread=[[len(x) for x in r] for r in res]
print("own pick terms P",np.round(terms[:,0],3).tolist(),"S",np.round(terms[:,1],3).tolist())
print("run pick terms P",np.round(np.array(pk["pick_terms"])[:,0],3).tolist(),"S",np.round(np.array(pk["pick_terms"])[:,1],3).tolist())
print("max |diff| %.4f s; readings per station-phase"%np.abs(terms-np.array(pk["pick_terms"])).max(),nread)
if len(sys.argv)>2 and sys.argv[2]=="terms": sys.exit()
lo,hi=(0.0,840.0) if SEG=="cal" else (840.0,1960.0)
h=json.load(open(RUN+f"/hyb_{SEG}.json"))
base=h["BASE"]
atoms=[a[:5] for a in json.load(open(RUN+f"/p1_{SEG}.json"))["atoms"]["NOCAP"] if lo<=a[3]<hi]
# own dedupe 5 km / 1.0 s
keep=[]
for a in sorted(atoms,key=lambda r:-r[4]):
    if all(math.hypot(a[0]-k[0],a[1]-k[1])>5.0 or abs(a[3]-k[3])>1.0 for k in keep): keep.append(a)
print("P1 after own dedupe",len(keep),"run",len(h["P1"]),"identical",keep==h["P1"])
p1=keep
def pct(x):
    x=np.asarray(x,float); r=np.argsort(np.argsort(x,kind="stable"),kind="stable"); return (r+1)/len(x)
wp=pct([a[4] for a in p1])
# waveform percentile of BASE events: not recomputed (needs the engine); recovered from the run's HYB output where
# the BASE event survives, else unknown. First pass: my own pick-support counts.
def support(cands,terms_):
    used=[{ph:np.zeros(len(PA[i][ph]),bool) for ph in "PS"} for i in range(len(names))]
    out=[]
    for c in sorted(cands,key=lambda r:-r[4]):
        t2=T(c[:3]); n=0
        for i in range(len(names)):
            for b,ph in enumerate("PS"):
                p=PA[i][ph]
                if not len(p): continue
                d=np.abs(p-(c[3]+t2[b,i]+terms_[i,b])); d[used[i][ph]]=np.inf; k=int(np.argmin(d))
                if d[k]<=0.5: used[i][ph][k]=True; n+=1
        out.append(c[:4]+[n,c[5]])
    return out
WB=None
if len(sys.argv)>2: WB=json.load(open(sys.argv[2]))
wb=WB if WB is not None else [0.0]*len(base)
cands=[e[:4]+[e[4]+100.0,wb[i]] for i,e in enumerate(base)]+[a[:4]+[wp[i],wp[i]] for i,a in enumerate(p1)]
sc=support(cands,terms)
dets=sorted([c[:4]+[c[4]+0.5*c[5]] for c in sc],key=lambda r:-r[4])
hyb=[]
for r in dets:
    if all(math.hypot(r[0]-k[0],r[1]-k[1])>5.0 or abs(r[3]-k[3])>1.5 for k in hyb): hyb.append(r)
# compare
run=h["HYB"]
print("own HYB n",len(hyb),"run HYB n",len(run))
def key(r): return (round(r[0],4),round(r[1],4),round(r[2],4),round(r[3],4))
ro={key(r):r for r in run}; mo={key(r):r for r in hyb}
print("common",len(set(ro)&set(mo)),"only run",len(set(ro)-set(mo)),"only own",len(set(mo)-set(ro)))
dn=[(math.floor(ro[k][4]),math.floor(mo[k][4])) for k in set(ro)&set(mo)]
print("integer support equal on common:",sum(a==b for a,b in dn),"of",len(dn))
for k in sorted(set(ro)^set(mo),key=lambda k:k[3])[:20]: print("   diff",k,ro.get(k,[None]*5)[4],mo.get(k,[None]*5)[4])
tg={i for i,e in enumerate(truth) if lo<=e["t"]<hi and e["M"]>=0.5}
for nm,L in (("own",hyb),("run",run)):
    print(nm,"LOC+T 10/30/60:",[len(op(L,lo,hi,0.5,r,tg)[0]) for r in (10,30,60)],"LOC:",[len(op(L,lo,hi,3.0,r,tg)[0]) for r in (10,30,60)])
# recover the BASE waveform percentiles from the run where a BASE event survived
bfr=[]
for e in base:
    k=key(e); 
    if k in ro: bfr.append(round(2*(ro[k][4]-math.floor(ro[k][4]))*len(base)))
print("recovered BASE waveform ranks (x n_base):",sorted(bfr))
json.dump(hyb,open(f"/home/claude/work/r1_verif_indep/own_hyb_{SEG}.json","w"))
