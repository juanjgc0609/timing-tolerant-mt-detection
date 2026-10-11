import json, time, numpy as np, warnings, sys
warnings.filterwarnings("ignore")
sys.path.insert(0,"/home/claude/work/gamma_src")
from task0_common import generate, make_network
from task1o_consistent import draw_sources, score_events
from task1r_fair import SCEN, degrade2
from task1r_baselines import run_gamma, GAMMA_GRID
import os; os.chdir("/home/claude/work"); seed=2024; NW=int(sys.argv[2]) if len(sys.argv)>2 else 40
for scen in sys.argv[1].split(","):
    kind,nst,nclut=SCEN[scen]
    rng=np.random.default_rng(seed+41); rng_d=np.random.default_rng(seed+43)
    def gen(base,k):
        S=make_network(kind,nst,seed=seed+base+k); srcs,taus=draw_sources(rng)
        D,p,tr=generate(S,srcs,taus,n_clutter=nclut,rng=rng); p2,lab=degrade2(S,srcs,p,"consistent",rng_d)
        return S,[[(t,b) for (t,b) in pl if t>-100] for pl in p2],tr,lab
    for k in range(60): gen(1000,k)
    ev=[gen(5000,k) for k in range(NW)]
    ntrue=sum(len(e[2]) for e in ev)
    # number of true sources with >=3 detected picks  (labels)
    npk=[]
    for S,p,tr,lab in ev:
        for j in range(len(tr)):
            npk.append(sum(1 for ll in lab for l in ll if l==j+1))
    npk=np.array(npk)
    print(scen,"true sources",ntrue,"detected picks per source: mean %.1f; frac with >=3: %.2f, >=4: %.2f, >=5: %.2f"%(npk.mean(),(npk>=3).mean(),(npk>=4).mean(),(npk>=5).mean()))
    t0=time.time()
    for ci,cfg in enumerate(GAMMA_GRID):
        h5=f5=h10=f10=0; dmiss=[]
        for S,p,tr,lab in ev:
            evs=run_gamma(S,p,**cfg)
            a,b=score_events(evs,tr); h5+=a; f5+=b
            # tolerance 10 km
            hit=set(); nf=0
            for e in evs:
                if not tr: nf+=1; continue
                dd=[np.linalg.norm(e[:2]-s[:2]) for s in tr]; j=int(np.argmin(dd))
                if dd[j]<=10: hit.add(j)
                else: nf+=1
            h10+=len(hit); f10+=nf
            for j,s in enumerate(tr):
                d=min([np.linalg.norm(e[:2]-s[:2]) for e in evs] or [np.inf])
                if d>5: dmiss.append(d)
        dmiss=np.array(dmiss)
        print(f"  cfg {cfg}: tol5 rec {h5/ntrue:.2f} fr {f5/NW:.2f} | tol10 rec {h10/ntrue:.2f} fr {f10/NW:.2f} | missed@5: {len(dmiss)}, of which nearest GaMMA event in (5,10] km: {np.sum(dmiss<=10)}, (10,20]: {np.sum((dmiss>10)&(dmiss<=20))}, none: {np.sum(~np.isfinite(dmiss))}")
    print("  time",time.time()-t0)
