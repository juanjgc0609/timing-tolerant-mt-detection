import json, math, numpy as np, sys
sys.path.insert(0,"/home/claude/work/r1_verif_indep")
from v3_endpoint import labels, truth, TM
RUN="/home/claude/work/r1/run"
d=json.load(open(RUN+"/hyb_eval.json"))
lo,hi=840.0,1960.0; BL=60.0; nb=int(math.ceil((hi-lo)/BL))
tg=np.array([i for i,e in enumerate(truth) if lo<=e["t"]<hi and e["M"]>=0.5]); tblk=((np.array([truth[i]["t"] for i in tg])-lo)//BL).astype(int)
def prep(v):
    v=[x for x in v if lo<=x[3]<hi]; lab=labels(v,0.5)
    s=np.array([x[4] for x in v]); blk=((np.array([x[3] for x in v])-lo)//BL).astype(int) if v else np.zeros(0,int)
    kind=np.array([l[0] for l in lab]); idx=np.array([l[1] for l in lab])
    return s,blk,kind,idx
P={"BASE":[prep(d["BASE"])],"HYB":[prep(d["HYB"])],"P1":[prep(d["P1"])],"BESTPICK":[prep(v) for v in d["BESTPICK"].values()]}
pos={int(t):k for k,t in enumerate(tg)}
def rec(sets,w,seconds):
    cap=int(math.floor(30*seconds/3600+1e-9)); tw=w[tblk]; best=0.0
    for s,blk,kind,idx in sets:
        if not len(s): continue
        f=kind=="false"; fs=s[f]; fw=w[blk[f]]; o=np.argsort(-fs,kind="stable")
        # threshold: largest set of top false scores with weighted count <= cap; ties at threshold excluded
        thr=-np.inf; cum=0.0
        for j in o:
            if cum+fw[j]>cap: thr=fs[j]; break
            cum+=fw[j]
        hit=(kind=="hit")&(s>thr)
        got=np.zeros(len(tg),bool)
        for i in idx[hit]:
            if i in pos: got[pos[i]]=True
        best=max(best,float((got*tw).sum()/tw.sum()))
    return best
w1=np.ones(nb); print("point (1120 s):",{m:round(rec(P[m],w1,1120.0),3) for m in P})
rng=np.random.default_rng(777); B=2000
for mode in ("seconds=sum(w)*60","seconds=sum of actual block lengths"):
    vals={("HYB","BASE"):[],("HYB","BESTPICK"):[],("P1","BASE"):[]}
    blen=np.array([min(BL,hi-lo-b*BL) for b in range(nb)])
    rng=np.random.default_rng(777)
    for _ in range(B):
        w=np.bincount(rng.integers(0,nb,nb),minlength=nb).astype(float)
        sec=w.sum()*BL if mode.startswith("seconds=sum(w)") else float((w*blen).sum())
        r={m:rec(P[m],w,sec) for m in P}
        for a,b in vals: vals[(a,b)].append(r[a]-r[b])
    print(mode)
    for (a,b),v in vals.items():
        v=np.array(v); print("  %s-%s median %+.3f 95%% [%+.3f, %+.3f]  P(diff>0) %.3f"%(a,b,np.median(v),np.percentile(v,2.5),np.percentile(v,97.5),(v>0).mean()))
