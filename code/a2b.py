"""
A2: is a proper subset of the physical channels statistically indistinguishable from the full
stack?   A1: is the ordering of channels by marginal contribution stable across cells?

Replicated over 4 geometry families x 2 clutter levels x 6 realizations.
Pick matching and nearest-pick lookup are precomputed per station as arrays indexed by time
bin, so scoring a candidate is O(1) per arrival and the replication is affordable.

Endpoints per realization:
  false_event = 1 if some hypothesis >5 km from the truth scores >= the truth
  margin      = (score at truth - best far score) / sd of the score field
"""
import numpy as np, statistics
from scipy.stats import spearmanr
rng=np.random.default_rng(101)
FC=2.0; VS=3.0; VP=1.73*VS; DT=1/(8*FC); T=np.arange(0,40,DT); NT=len(T)
SNR=8.0; TOL=0.6; AP,AS=1.0,1.5; VMIN=2.0
SRC=np.array([0.,0.,8.]); TAU=6.0
def ricker(t):
    a=(np.pi*FC*t)**2; return (1-2*a)*np.exp(-a)
def pol(src,r,b):
    v=r-src; n=v/np.linalg.norm(v)
    if b=='P': return n
    ez=np.array([0.,0.,1.]); sv=ez-np.dot(ez,n)*n; nn=np.linalg.norm(sv)
    return sv/nn if nn>1e-9 else np.array([1.,0.,0.])
def geom(kind):
    a=np.linspace(0,2*np.pi,8,endpoint=False)+rng.uniform(0,1)
    if kind=='spread':   r=np.array([6.,11.,18.,27.,36.,45.,54.,62.])
    elif kind=='ring':   r=np.full(8,28.)
    elif kind=='onesided':a=rng.uniform(0,np.pi*.7,8); r=45*np.sqrt(rng.uniform(.05,1,8))
    else:                a=rng.uniform(0,2*np.pi,8); r=45*np.sqrt(rng.uniform(.05,1,8))
    return np.stack([r*np.cos(a),r*np.sin(a),np.zeros(8)],axis=1)

def gen(STA,nc,scatter_deg=90.):
    """scatter_deg = angular scatter of clutter polarization around the parent-event direction.
       0 deg  -> clutter perfectly mimics the true event's particle motion (coda-like, hardest)
       90 deg -> clutter polarization is random (easiest for C-IV)"""
    D=np.zeros((8,3,NT)); pk=[[] for _ in range(8)]
    d=np.linalg.norm(STA-SRC,axis=1)
    for i in range(8):
        for b,A,V in (('P',AP,VP),('S',AS,VS)):
            ta=d[i]/V+TAU; pk[i].append(ta)
            D[i]+=SNR*A*np.outer(pol(SRC,STA[i],b),ricker(T-ta))/d[i]
    sc=np.radians(scatter_deg)
    for _ in range(nc):
        i=rng.integers(8); tc=rng.uniform(2,36)
        base=pol(SRC,STA[i],'P' if rng.random()<.5 else 'S')
        pert=rng.normal(size=3); pert-=np.dot(pert,base)*base
        nn=np.linalg.norm(pert); pert=pert/nn if nn>1e-9 else np.array([1.,0,0])
        ang=abs(rng.normal(0,sc/2)) if sc>0 else 0.
        q=np.cos(ang)*base+np.sin(ang)*pert; q/=np.linalg.norm(q)
        D[i]+=SNR*.35*np.outer(q,ricker(T-tc))/np.median(np.linalg.norm(STA,axis=1)); pk[i].append(tc)
    D+=rng.normal(scale=SNR*.02,size=D.shape)
    return D,pk

def precompute(D,pk):
    W=int(.6/DT); AMP=np.zeros((8,NT)); PD=np.zeros((8,NT,3))
    HIT=np.zeros((8,NT),bool); NEAR=np.zeros((8,NT))
    for i in range(8):
        for j in range(NT):
            seg=D[i][:,max(0,j-W//2):min(NT,j+W//2)]; C=seg@seg.T
            AMP[i,j]=np.sqrt(np.trace(C)); PD[i,j]=np.linalg.eigh(C)[1][:,-1]
        p=np.sort(np.array(pk[i]))
        idx=np.searchsorted(p,T); idx=np.clip(idx,0,len(p)-1)
        cand=np.stack([p[np.clip(idx-1,0,len(p)-1)],p[idx]])
        best=cand[np.argmin(np.abs(cand-T),axis=0),np.arange(NT)]
        NEAR[i]=best; HIT[i]=np.abs(best-T)<TOL
    return AMP,PD,HIT,NEAR

def score(src,tau,STA,DIJ,pre,use):
    AMP,PD,HIT,NEAR=pre
    d=np.linalg.norm(STA-src,axis=1)
    tt=np.stack([d/VP,d/VS])+tau
    idx=np.clip((tt/DT).astype(int),0,NT-1)
    st=np.arange(8)
    hit=HIT[st[None,:],idx]
    n=hit.sum()
    if n==0: return 0.
    if 'CIII' in use:
        cl=NEAR[st[None,:],idx]
        for b in range(2):
            ii=np.where(hit[b])[0]
            if len(ii)>1:
                dt=np.abs(cl[b][ii][:,None]-cl[b][ii][None,:])
                if np.any(dt>DIJ[np.ix_(ii,ii)]/VMIN+1e-9): return 0.
        both=hit[0]&hit[1]
        if np.any(cl[1][both]<=cl[0][both]): return 0.
    s=n/16.
    if 'CII' in use:
        m=[AMP[i,idx[b,i]]*d[i]/A for b,(bb,A) in enumerate((('P',AP),('S',AS))) for i in range(8) if hit[b,i]]
        M0=np.median(m) if m else 0.
        nf=np.percentile(AMP,10)
        for b,(bb,A) in enumerate((('P',AP),('S',AS))):
            for i in range(8):
                if hit[b,i]: continue
                if M0*A/d[i]>3*nf and AMP[i,idx[b,i]]<1.5*nf: return 0.
    if 'CIV' in use:
        cs=[abs(np.dot(pol(src,STA[i],bb),PD[i,idx[b,i]])) for b,bb in enumerate(('P','S')) for i in range(8) if hit[b,i]]
        s*= np.mean(cs)**2 if cs else 0.
    if 'CV' in use:
        m=np.array([AMP[i,idx[b,i]]*d[i]/A for b,(bb,A) in enumerate((('P',AP),('S',AS))) for i in range(8) if hit[b,i]])
        s*= 1./(1.+m.std()/max(m.mean(),1e-9)) if len(m)>=2 else .5
    return s

gx=np.arange(-21,21.1,3.5); gz=np.array([4.,8.,12.]); TAUS=np.arange(4.5,7.51,0.2)
GP=[np.array([x,y,z]) for x in gx for y in gx for z in gz]
FAR=np.array([np.linalg.norm(g-SRC)>5. for g in GP])
def evaluate(STA,DIJ,pre,use):
    v=np.array([max(score(g,tt,STA,DIJ,pre,use) for tt in TAUS) for g in GP])
    st=max(score(SRC,tt,STA,DIJ,pre,use) for tt in TAUS)
    bf=v[FAR].max()
    return int(bf>=st), (st-bf)/max(v.std(),1e-9)

SETS=[('TT',set()),('C-IV',{'CIV'}),('MIN',{'CIII','CIV'}),('ALL',{'CII','CIII','CIV','CV'})]
REPS=8
print("false-event rate vs POLARIZATION COHERENCE of the clutter")
print("(0 deg = clutter mimics the true event's particle motion; 90 deg = random)\n")
print(f"{'geometry':>9} {'scatter':>8} | " + " ".join(f"{n:>6}" for n,_ in SETS))
print("-"*46)
for kind in ['spread','onesided']:
    STA=geom(kind); DIJ=np.linalg.norm(STA[:,None,:]-STA[None,:,:],axis=2)
    for sd in [0.,15.,30.,60.,90.]:
        fe={n:0 for n,_ in SETS}
        for _ in range(REPS):
            D,pk=gen(STA,30,scatter_deg=sd); pre=precompute(D,pk)
            for n,u in SETS:
                f,m=evaluate(STA,DIJ,pre,u); fe[n]+=f
        print(f"{kind:>9} {sd:7.0f}d | " + " ".join(f"{fe[n]/REPS:6.2f}" for n,_ in SETS))
    print()
