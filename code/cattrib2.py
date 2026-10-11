"""H3 re-run with C-I and C-II implemented so they can actually bind.

C-II is now the real test: from the MATCHED arrivals estimate the implied scalar moment;
then for every predicted arrival that is NOT matched, ask whether that moment would have
produced a detectable signal there.  If yes and the record is silent, the hypothesis is
infeasible.  This is what a human analyst does with non-detecting stations.
C-I now uses a non-trivial floor gamma_min = 20th percentile over the candidate grid.
"""
import numpy as np
rng=np.random.default_rng(20)
FC=2.0; VS=3.0; VP=1.73*VS; DT=1/(8*FC); T=np.arange(0,30,DT); NT=len(T)
SRC=np.array([0.,0.,8.]); TAU=5.0; SNR=8.0; TOL=0.5; AP,AS=1.0,1.5
th=np.linspace(0,2*np.pi,6,endpoint=False)+.3
STA=np.stack([25*np.cos(th),25*np.sin(th),np.zeros(6)],axis=1)
def ricker(t):
    a=(np.pi*FC*t)**2; return (1-2*a)*np.exp(-a)
def polar(src,r,b):
    v=r-src; n=v/np.linalg.norm(v)
    if b=='P': return n
    ez=np.array([0.,0.,1.]); sv=ez-np.dot(ez,n)*n; nn=np.linalg.norm(sv)
    return sv/nn if nn>1e-9 else np.array([1.,0.,0.])
def make_data(nc):
    D=np.zeros((6,3,NT)); picks=[[] for _ in range(6)]
    d=np.linalg.norm(STA-SRC,axis=1)
    for i in range(6):
        for b,A,V in (('P',AP,VP),('S',AS,VS)):
            ta=d[i]/V+TAU; picks[i].append(ta)
            D[i]+=SNR*A*np.outer(polar(SRC,STA[i],b),ricker(T-ta))/d[i]
    for _ in range(nc):
        i=rng.integers(6); tc=rng.uniform(2,26); v=rng.normal(size=3); v/=np.linalg.norm(v)
        D[i]+=SNR*0.35*np.outer(v,ricker(T-tc))/np.median(d); picks[i].append(tc)
    D+=rng.normal(scale=SNR*0.02,size=D.shape)
    return D,picks
def feats(D):
    W=int(0.6/DT); E=np.zeros((6,NT)); PD=np.zeros((6,NT,3))
    for i in range(6):
        for j in range(NT):
            seg=D[i][:,max(0,j-W//2):min(NT,j+W//2)]; C=seg@seg.T
            E[i,j]=np.trace(C); PD[i,j]=np.linalg.eigh(C)[1][:,-1]
    return np.sqrt(E),PD                                  # amplitude, principal direction
def glead(src):
    M=np.zeros((4,4))
    for r in STA:
        d=np.linalg.norm(r-src); n=(r-src)/d
        for V,A in ((VP,AP),(VS,AS)):
            v=np.concatenate([-n/V,[1.]]); M+=(A/d)**2*np.outer(v,v)
    return np.linalg.eigvalsh(M)[0]

gx=np.arange(-24,24.1,3.0); gz=np.array([4.,8.,12.,16.]); taus=np.arange(3.,7.01,0.1)
GRID=[np.array([x,y,z]) for x in gx for y in gx for z in gz]
GL=np.array([glead(s) for s in GRID]); GMIN=np.percentile(GL,20)

def score(src,tau,AMP,PD,picks,use,noise_floor):
    d=np.linalg.norm(STA-src,axis=1); tt=np.stack([d/VP,d/VS])+tau
    idx=np.clip((tt/DT).astype(int),0,NT-1)
    hit=np.array([[any(abs(tt[bi,i]-pk)<TOL for pk in picks[i]) for i in range(6)] for bi in range(2)])
    if hit.sum()==0: return -1
    s=hit.sum()/12.0
    m=[AMP[i,idx[bi,i]]*d[i]/A for bi,(b,A) in enumerate((('P',AP),('S',AS))) for i in range(6) if hit[bi,i]]
    M0=np.median(m) if m else 0.0
    if 'CII' in use:                      # silence as exclusion
        for bi,(b,A) in enumerate((('P',AP),('S',AS))):
            for i in range(6):
                if hit[bi,i]: continue
                pred=M0*A/d[i]
                if pred>3*noise_floor and AMP[i,idx[bi,i]]<1.5*noise_floor:
                    return -1
    if 'CIV' in use:
        cs=[abs(np.dot(polar(src,STA[i],b),PD[i,idx[bi,i]])) for bi,b in enumerate(('P','S')) for i in range(6) if hit[bi,i]]
        s*= np.mean(cs)**2 if cs else 0.
    if 'CV' in use:
        if len(m)>=2:
            m=np.array(m); s*=1./(1.+m.std()/max(m.mean(),1e-9))
        else: s*=0.5
    return s

def run(D,picks,use):
    AMP,PD=feats(D); nf=np.percentile(AMP,10)
    vals=[]; 
    for gi,src in enumerate(GRID):
        if 'CI' in use and GL[gi]<GMIN: continue
        b=max((score(src,tau,AMP,PD,picks,use,nf) for tau in taus),default=-1)
        vals.append((np.linalg.norm(src-SRC),b))
    V=np.array(vals); ok=V[:,1]>-0.5
    st=max(score(SRC,tau,AMP,PD,picks,use,nf) for tau in taus)
    far=(V[:,0]>5.)&ok
    bf=V[far,1].max() if far.sum() else 0.
    sd=V[ok,1].std() if ok.sum()>1 else 1.
    return st,bf,(st-bf)/max(sd,1e-9),int(ok.sum())

sets=[('TT only',set()),('ALL',{'CI','CII','CIV','CV'}),('no C-I',{'CII','CIV','CV'}),
      ('no C-II',{'CI','CIV','CV'}),('no C-IV',{'CI','CII','CV'}),('no C-V',{'CI','CII','CIV'}),
      ('C-IV only',{'CIV'}),('C-II only',{'CII'})]
print(f"gamma_min (20th pct) = {GMIN:.4f}\n")
print(f"{'clutter':>7} {'set':>10} | {'true':>6} {'bestfalse':>9} {'margin':>8} {'feasible':>9} {'d margin':>9}")
print("-"*70)
for nc in [15,40]:
    D,picks=make_data(nc); ref=None
    for name,use in sets:
        st,bf,mg,nf_=run(D,picks,use)
        if name=='ALL': ref=mg
        dd = "" if ref is None or name in('TT only','ALL') else f"{mg-ref:+9.2f}"
        print(f"{nc:7d} {name:>10} | {st:6.3f} {bf:9.3f} {mg:8.2f} {nf_:9d} {dd:>9}")
    print()
