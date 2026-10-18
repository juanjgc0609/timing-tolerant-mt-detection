"""Own recomputation of the R1 endpoint from raw detection lists (no project imports)."""
import json, math, sys
import numpy as np
from math import comb
W="/home/claude/work/r1_verif_indep/"
RUN=sys.argv[1] if len(sys.argv)>1 else "/home/claude/work/r1/run"
truth=json.load(open(W+"truth_own.json")); jma=json.load(open(W+"jma_own.json"))
cal=[j for j in jma if 0<=j["t"]<840]
lat0=sum(j["lat"] for j in cal)/len(cal); lon0=sum(j["lon"] for j in cal)/len(cal)
KX=111.19*math.cos(math.radians(lat0)); KY=111.19
for e in truth: e["x"]=(e["lon"]-lon0)*KX; e["y"]=(e["lat"]-lat0)*KY
TX=np.array([e["x"] for e in truth]); TY=np.array([e["y"] for e in truth]); TT=np.array([e["t"] for e in truth]); TM=np.array([e["M"] for e in truth])
MT=0.5
def labels(dets, tol):
    """per detection: ('hit', idx) if within 5 km/tol of an event with M>=0.5 (nearest such), ('neutral',idx) if only small, ('false',-1)"""
    out=[]
    for d in dets:
        dd=np.hypot(TX-d[0],TY-d[1]); ok=(dd<=5.0)&(np.abs(TT-d[3])<=tol+1e-9)
        big=ok&(TM>=MT)
        if big.any(): out.append(("hit",int(np.where(big,dd,np.inf).argmin())))
        elif ok.any(): out.append(("neutral",int(np.where(ok,dd,np.inf).argmin())))
        else: out.append(("false",-1))
    return out
def op(dets, lo, hi, tol, rate, targets):
    dets=[d for d in dets if lo<=d[3]<hi]
    hours=(hi-lo)/3600.0; cap=int(math.floor(rate*hours+1e-9))
    lab=labels(dets,tol)
    fs=sorted([d[4] for d,l in zip(dets,lab) if l[0]=="false"],reverse=True)
    thr=fs[cap] if cap<len(fs) else -np.inf
    got={l[1] for d,l in zip(dets,lab) if l[0]=="hit" and d[4]>thr} & targets
    nfalse=sum(1 for d,l in zip(dets,lab) if l[0]=="false" and d[4]>thr)
    return got,thr,nfalse,len(dets)
def mcn(b,c):
    n=b+c
    return 1.0 if n==0 else sum(comb(n,k) for k in range(b,n+1))/2.0**n
if __name__=="__main__":
    d=json.load(open(RUN+"/hyb_eval.json"))
    lo,hi=840.0,1960.0
    tg={i for i,e in enumerate(truth) if lo<=e["t"]<hi and e["M"]>=MT}
    tj={i for i,e in enumerate(truth) if lo<=e["t"]<hi and e["jma"]}
    print("targets",len(tg),"jma",len(tj),"centre",lat0,lon0)
    meth={"BASE":[("BASE",d["BASE"])],"HYB":[("HYB",d["HYB"])],"P1":[("P1",d["P1"])],"BESTPICK":list(d["BESTPICK"].items())}
    R={}
    for m,L in meth.items():
        row=[]
        for rule,tol in (("LOC+T",0.5),("LOC",3.0)):
            for rate in (10,30,60):
                best=None
                for k,v in L:
                    g,thr,nf,n=op(v,lo,hi,tol,rate,tg)
                    if best is None or len(g)>len(best[1]): best=(k,g,thr,nf)
                R[(m,rule,rate)]=best; row.append("%.3f (%2d)"%(len(best[1])/len(tg),len(best[1])))
        print("%-9s"%m,"  ".join(row), "n det",len([x for x in L[0][1] if lo<=x[3]<hi]) if m!="BESTPICK" else "")
    # ties for BESTPICK at LOC+T 30
    print("BESTPICK configs at LOC+T 30/h:")
    for k,v in d["BESTPICK"].items():
        g,thr,nf,n=op(v,lo,hi,0.5,30,tg); print("   %-45s %2d  nfalse %d thr %s n %d"%(k,len(g),nf,thr,n))
    bk=R[("BESTPICK","LOC+T",30)][0]; print("chosen BESTPICK",bk)
    tim=lambda s: ", ".join("t %.1f M %.2f"%(truth[i]["t"],truth[i]["M"]) for i in sorted(s,key=lambda i: truth[i]["t"]))
    for tag,a,b in (("R1-a","HYB","BASE"),("R1-b","HYB","BESTPICK"),("R1-c","P1","BASE")):
        ga,gb=R[(a,"LOC+T",30)][1],R[(b,"LOC+T",30)][1]
        p=mcn(len(ga-gb),len(gb-ga))
        print("%s %s-%s: %d vs %d, diff %+.3f, disc %d/%d, p=%.4f %s"%(tag,a,b,len(ga),len(gb),(len(ga)-len(gb))/len(tg),len(ga-gb),len(gb-ga),p,"MET" if p<=0.025 and len(ga)>len(gb) else "NOT MET"))
        print("    only",a,":",tim(ga-gb)); print("    only",b,":",tim(gb-ga))
    # all tied BESTPICK configs: McNemar vs HYB
    gh=R[("HYB","LOC+T",30)][1]
    for k,v in d["BESTPICK"].items():
        g=op(v,lo,hi,0.5,30,tg)[0]
        if len(g)==len(R[("BESTPICK","LOC+T",30)][1]): print("   tied",k,"disc",len(gh-g),len(g-gh),"p %.4f"%mcn(len(gh-g),len(g-gh)))
    # JMA-only: targets = jma-flag truth events; detections hitting other truth events not false
    print("JMA-only")
    for m,L in meth.items():
        row=[]
        for rate in (10,30,60):
            best=0
            for k,v in L:
                dets=[x for x in v if lo<=x[3]<hi]; lab=labels(dets,0.5)
                # registered: hit any truth event -> not false; hit counted if label event is jma
                # use nearest within tol among all truth (any M) for this reading
                lab2=[]
                for x in dets:
                    dd=np.hypot(TX-x[0],TY-x[1]); ok=(dd<=5.0)&(np.abs(TT-x[3])<=0.5+1e-9)
                    lab2.append(int(np.where(ok,dd,np.inf).argmin()) if ok.any() else -1)
                cap=int(math.floor(rate*(hi-lo)/3600+1e-9))
                fs=sorted([x[4] for x,l in zip(dets,lab2) if l<0],reverse=True); thr=fs[cap] if cap<len(fs) else -np.inf
                got={l for x,l in zip(dets,lab2) if l>=0 and x[4]>thr}&tj
                # also any-jma-within-tol variant
                got2=set()
                for x in dets:
                    if x[4]<=thr: continue
                    dd=np.hypot(TX-x[0],TY-x[1]); ok=(dd<=5.0)&(np.abs(TT-x[3])<=0.5+1e-9)
                    got2|=set(np.where(ok)[0].tolist())&tj
                best=max(best,len(got))
                if len(got2)!=len(got): print("   note",m,k,rate,"nearest-label",len(got),"any-within",len(got2))
            row.append("%d/%d=%.3f"%(best,len(tj),best/len(tj)))
        print("  ",m,row)
    # deployable threshold
    dc=json.load(open(RUN+"/hyb_cal.json"))
    tgc={i for i,e in enumerate(truth) if 0<=e["t"]<840 and e["M"]>=MT}
    print("deployable (CAL targets %d)"%len(tgc))
    for m in ("HYB","BASE","P1"):
        gc,thr,nfc,nc=op(dc[m],0.0,840.0,0.5,30,tgc)
        dets=[x for x in d[m] if lo<=x[3]<hi]; lab=labels(dets,0.5)
        got={l[1] for x,l in zip(dets,lab) if l[0]=="hit" and x[4]>thr}&tg
        nf=sum(1 for x,l in zip(dets,lab) if l[0]=="false" and x[4]>thr)
        print("   %s CAL recall %d/%d thr %.3f -> EVAL %d (%.3f), %d false (%.1f/h)"%(m,len(gc),len(tgc),thr,len(got),len(got)/len(tg),nf,nf/((hi-lo)/3600)))
    json.dump({"%s|%s|%d"%k:[v[0],sorted(v[1]),float(v[2]) if np.isfinite(v[2]) else None,v[3]] for k,v in R.items()},open(W+"own_endpoint.json","w"))
