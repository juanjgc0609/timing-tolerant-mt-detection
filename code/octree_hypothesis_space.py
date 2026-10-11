"""
Resolution-matched octree over the source hypothesis space.

Three search strategies compared at equal task:
  A  uniform fine grid                       (baseline, expensive)
  B  score-refined octree                    (Qseek / NonLinLoc style: refine where score is high)
  C  observability-matched octree            (refine until cell size <= c * Delta_res(x),
                                              with Delta_res from the LOCAL coercivity constant)
C never consults the data to decide where to refine: the tree is built from geometry alone,
before any waveform is read.  That is the difference from B.
Metrics: number of score evaluations, epicentral error, false-event rate.
"""
import numpy as np, statistics
rng=np.random.default_rng(2024)
FC=2.0; VS=3.0; VP=1.73*VS; DT=1/(8*FC); T=np.arange(0,40,DT); NT=len(T)
SNR=7.0; TOL=0.6; AP,AS=1.0,1.5; VMIN=2.0; TAU=6.0
# ONE-SIDED network: coverage, and therefore gamma, varies strongly across the region
a=np.array([0.15,0.45,0.75,1.05,1.35,1.65,1.95,2.25])
r=np.array([10.,14.,19.,25.,32.,39.,47.,55.])
STA=np.stack([r*np.cos(a),r*np.sin(a),np.zeros(8)],axis=1)
DIJ=np.linalg.norm(STA[:,None,:]-STA[None,:,:],axis=2)
def ricker(t):
    x=(np.pi*FC*t)**2; return (1-2*x)*np.exp(-x)
def pol(src,rr,b):
    v=rr-src; n=v/np.linalg.norm(v)
    if b=='P': return n
    ez=np.array([0.,0.,1.]); sv=ez-np.dot(ez,n)*n; nn=np.linalg.norm(sv)
    return sv/nn if nn>1e-9 else np.array([1.,0.,0.])
def gamma_lead(src):
    M=np.zeros((4,4))
    for rr in STA:
        d=np.linalg.norm(rr-src); n=(rr-src)/d
        for V,A in ((VP,AP),(VS,AS)):
            v=np.concatenate([-n/V,[1.]]); M+=(A/d)**2*np.outer(v,v)
    return max(np.linalg.eigvalsh(M)[0],1e-300)
def gen(src,nc=30):
    D=np.zeros((8,3,NT)); pk=[[] for _ in range(8)]
    d=np.linalg.norm(STA-src,axis=1)
    for i in range(8):
        for b,A,V in (('P',AP,VP),('S',AS,VS)):
            ta=d[i]/V+TAU; pk[i].append(ta)
            D[i]+=SNR*A*np.outer(pol(src,STA[i],b),ricker(T-ta))/d[i]
    for _ in range(nc):
        i=rng.integers(8); tc=rng.uniform(2,36); q=rng.normal(size=3); q/=np.linalg.norm(q)
        D[i]+=SNR*.35*np.outer(q,ricker(T-tc))/28.; pk[i].append(tc)
    D+=rng.normal(scale=SNR*.02,size=D.shape); return D,pk
def pre_(D,pk):
    """vectorized: sliding covariance via cumulative sums, batched 3x3 eigendecomposition"""
    W=int(.6/DT); h=W//2
    PD=np.zeros((8,NT,3)); HIT=np.zeros((8,NT),bool); NEAR=np.zeros((8,NT))
    for i in range(8):
        x=D[i]                                    # (3,NT)
        prod=np.einsum('at,bt->abt',x,x)          # (3,3,NT)
        cs=np.concatenate([np.zeros((3,3,1)),np.cumsum(prod,axis=2)],axis=2)
        lo=np.clip(np.arange(NT)-h,0,NT); hi=np.clip(np.arange(NT)+h,0,NT)
        C=(cs[:,:,hi]-cs[:,:,lo]).transpose(2,0,1)   # (NT,3,3)
        w,V=np.linalg.eigh(C)
        PD[i]=V[:,:,-1]
        p=np.sort(np.array(pk[i])); idx=np.clip(np.searchsorted(p,T),0,len(p)-1)
        cand=np.stack([p[np.clip(idx-1,0,len(p)-1)],p[idx]])
        best=cand[np.argmin(np.abs(cand-T),axis=0),np.arange(NT)]
        NEAR[i]=best; HIT[i]=np.abs(best-T)<TOL
    return PD,HIT,NEAR
TAUS=np.arange(5.0,7.01,0.25)
POLC={}
def polcache(src):
    k=tuple(np.round(src,3))
    if k not in POLC:
        POLC[k]=np.stack([[pol(src,STA[i],b) for i in range(8)] for b in ('P','S')])
    return POLC[k]
def score(src,pre):
    PD,HIT,NEAR=pre; st=np.arange(8); best=0.
    d=np.linalg.norm(STA-src,axis=1); PC=polcache(src)
    for tau in TAUS:
        tt=np.stack([d/VP,d/VS])+tau; idx=np.clip((tt/DT).astype(int),0,NT-1)
        hit=HIT[st[None,:],idx]
        if hit.sum()==0: continue
        cl=NEAR[st[None,:],idx]; ok=True
        for b in range(2):
            ii=np.where(hit[b])[0]
            if len(ii)>1 and np.any(np.abs(cl[b][ii][:,None]-cl[b][ii][None,:])>DIJ[np.ix_(ii,ii)]/VMIN+1e-9): ok=False;break
        if not ok: continue
        s=hit.sum()/16.
        dp=np.abs(np.einsum('bij,bij->bi',PC,PD[st[None,:],idx]))
        cs=dp[hit]
        s*= cs.mean()**2 if cs.size else 0.
        best=max(best,s)
    return best


ROOT=(np.array([-24.,-24.,3.]), np.array([24.,24.,15.]))
GREF=gamma_lead(np.array([0.,0.,9.]))
KAPPA=0.12
def dres(src):
    """local resolution limit from the coercivity constant"""
    return KAPPA*VS/FC*np.sqrt(GREF/gamma_lead(src))

def uniform(pre,h):
    lo,hi=ROOT; n=0; best=(-1,None)
    for x in np.arange(lo[0]+h/2,hi[0],h):
        for y in np.arange(lo[1]+h/2,hi[1],h):
            for z in np.arange(lo[2]+h/2,hi[2],h):
                p=np.array([x,y,z]); sc=score(p,pre); n+=1
                if sc>best[0]: best=(sc,p)
    return best,n

def octree(pre,stop,keep_frac=0.25,maxlev=7,minlev=2):
    """branching is ALWAYS score-guided (keep the best fraction).
       stop='fixed'   -> descend a fixed number of levels
       stop='resmatch'-> stop a branch once its cell size reaches the LOCAL resolution limit"""
    lo,hi=ROOT; cells=[(lo,hi)]; n=0; best=(-1,None); done=[]
    for lev in range(maxlev):
        scored=[]
        for (a_,b_) in cells:
            c=(a_+b_)/2; sc=score(c,pre); n+=1
            if sc>best[0]: best=(sc,c)
            scored.append((sc,a_,b_,c))
        if stop=='resmatch':
            alive=[t for t in scored if np.max(t[2]-t[1])>dres(t[3])]
            done+=[t for t in scored if np.max(t[2]-t[1])<=dres(t[3])]
        else:
            alive=scored if lev<4 else []
        if not alive: break
        alive.sort(key=lambda t:-t[0])
        k=max(1,int(len(alive)*keep_frac)) if lev>=minlev else len(alive)
        new=[]
        for (_,a_,b_,c) in alive[:k]:
            m=(a_+b_)/2
            for dx in (0,1):
                for dy in (0,1):
                    for dz in (0,1):
                        na=np.array([a_[0] if dx==0 else m[0], a_[1] if dy==0 else m[1], a_[2] if dz==0 else m[2]])
                        nb=na+(b_-a_)/2; new.append((na,nb))
        cells=new
    return best,n

import statistics
gs=[gamma_lead(np.array([x,y,9.])) for x in np.arange(-20,21,5.) for y in np.arange(-20,21,5.)]
print(f"one-sided network: gamma varies {max(gs)/min(gs):.0f}x across the region\n")
print(f"{'strategy':>28} {'evals':>7} {'median err km':>14} {'false-event':>12}")
print("-"*66)
REPS=6
SRCS=[np.array([6.,-9.,9.]), np.array([-14.,11.,7.]), np.array([2.,3.,12.]), np.array([-5.,-16.,6.])]
for name,fn in [("uniform grid 4 km",lambda pr:uniform(pr,4.0)),
                ("uniform grid 2 km",lambda pr:uniform(pr,2.0)),
                ("octree, fixed depth",lambda pr:octree(pr,'fixed')),
                ("octree, observability stop",lambda pr:octree(pr,'resmatch'))]:
    errs=[];ns=[];fe=0
    for rep in range(REPS):
        S0=SRCS[rep%len(SRCS)]
        D,pk=gen(S0); pre=pre_(D,pk); POLC.clear()
        (sv,bp),n=fn(pre)
        e=float(np.hypot(bp[0]-S0[0],bp[1]-S0[1])); errs.append(e); ns.append(n); fe+= int(e>5.)
    print(f"{name:>28} {int(np.median(ns)):7d} {statistics.median(errs):14.2f} {fe/REPS:12.2f}")
