"""
Overlapping events: does the constraint stack recover TWO sources where travel-time
scoring merges them or ranks false hypotheses above the true ones?

Two true sources with controlled time separation, plus clutter.  Candidates scored on a grid
(single-source scoring, i.e. each hypothesis explains what it can).  Metrics:
  rank2   = rank of the WORSE-ranked of the two true sources in the score ordering
  nfalse  = number of hypotheses >5 km from BOTH true sources that outrank that source
Lower is better for both.
"""
import numpy as np
rng=np.random.default_rng(7)
FC=2.0; VS=3.0; VP=1.73*VS; DT=1/(8*FC); T=np.arange(0,40,DT); NT=len(T)
SNR=8.0; TOL=0.5; AP,AS=1.0,1.5; VMIN=2.0
th=np.linspace(0,2*np.pi,8,endpoint=False)+.3
RR=np.array([6.,11.,18.,27.,36.,45.,54.,62.])
STA=np.stack([RR*np.cos(th),RR*np.sin(th),np.zeros(8)],axis=1)
DIJ=np.linalg.norm(STA[:,None,:]-STA[None,:,:],axis=2)
def ricker(t):
    a=(np.pi*FC*t)**2; return (1-2*a)*np.exp(-a)
def polar(src,r,b):
    v=r-src; n=v/np.linalg.norm(v)
    if b=='P': return n
    ez=np.array([0.,0.,1.]); sv=ez-np.dot(ez,n)*n; nn=np.linalg.norm(sv)
    return sv/nn if nn>1e-9 else np.array([1.,0.,0.])
def gen(S1,T1,S2,T2,nclutter):
    D=np.zeros((8,3,NT)); picks=[[] for _ in range(8)]
    for src,tau in ((S1,T1),(S2,T2)):
        d=np.linalg.norm(STA-src,axis=1)
        for i in range(8):
            for b,A,V in (('P',AP,VP),('S',AS,VS)):
                ta=d[i]/V+tau; picks[i].append(ta)
                D[i]+=SNR*A*np.outer(polar(src,STA[i],b),ricker(T-ta))/d[i]
    for _ in range(nclutter):
        i=rng.integers(8); tc=rng.uniform(2,36); v=rng.normal(size=3); v/=np.linalg.norm(v)
        D[i]+=SNR*0.35*np.outer(v,ricker(T-tc))/np.median(RR); picks[i].append(tc)
    D+=rng.normal(scale=SNR*0.02,size=D.shape)
    return D,picks
def feats(D):
    W=int(0.6/DT); E=np.zeros((8,NT)); PD=np.zeros((8,NT,3))
    for i in range(8):
        for j in range(NT):
            seg=D[i][:,max(0,j-W//2):min(NT,j+W//2)]; C=seg@seg.T
            E[i,j]=np.trace(C); PD[i,j]=np.linalg.eigh(C)[1][:,-1]
    return np.sqrt(E),PD
gx=np.arange(-18,18.1,3.0); gz=np.array([4.,8.,12.]); taus=np.arange(3.,9.01,0.2)
GRID=[np.array([x,y,z]) for x in gx for y in gx for z in gz]

def score(src,AMP,PD,picks,use):
    d=np.linalg.norm(STA-src,axis=1); best=-1
    for tau in taus:
        tt=np.stack([d/VP,d/VS])+tau; idx=np.clip((tt/DT).astype(int),0,NT-1)
        hit=np.array([[any(abs(tt[b,i]-pk)<TOL for pk in picks[i]) for i in range(8)] for b in range(2)])
        if hit.sum()==0: continue
        # C-III medium-independent admissibility on the CLAIMED picks
        if 'CIII' in use:
            ok=True
            for b in range(2):
                ii=[i for i in range(8) if hit[b,i]]
                for a_ in range(len(ii)):
                    for c_ in range(a_+1,len(ii)):
                        if abs(tt[b,ii[a_]]-tt[b,ii[c_]])>DIJ[ii[a_],ii[c_]]/VMIN+1e-9: ok=False;break
                    if not ok: break
                if not ok: break
            if not ok: continue
        s=hit.sum()/16.
        if 'CIV' in use:
            cs=[abs(np.dot(polar(src,STA[i],b),PD[i,idx[bi,i]])) for bi,b in enumerate(('P','S')) for i in range(8) if hit[bi,i]]
            s*= np.mean(cs)**2 if cs else 0.
        if 'CV' in use:
            m=[AMP[i,idx[bi,i]]*d[i]/A for bi,(b,A) in enumerate((('P',AP),('S',AS))) for i in range(8) if hit[bi,i]]
            if len(m)>=2:
                m=np.array(m); s*=1./(1.+m.std()/max(m.mean(),1e-9))
            else: s*=0.5
        best=max(best,s)
    return best

def run(S1,T1,S2,T2,nc,use):
    D,picks=gen(S1,T1,S2,T2,nc); AMP,PD=feats(D)
    v=np.array([score(g,AMP,PD,picks,use) for g in GRID])
    s1=score(S1,AMP,PD,picks,use); s2=score(S2,AMP,PD,picks,use)
    far=np.array([min(np.linalg.norm(g-S1),np.linalg.norm(g-S2))>4. for g in GRID])
    worse=min(s1,s2)
    amb=int(((v>=worse*0.98)&far).sum())          # ambiguity set: far hypotheses that TIE or beat
    top=v.max()
    ambtop=int(((v>=top*0.98)&far).sum())
    return s1,s2,amb,ambtop,len(GRID)

S1=np.array([0.,0.,8.])
print(f"{'sep km':>5} {'dt s':>6} {'clut':>6} {'set':>11} | {'s_1':>7} {'s_2':>7} {'ambig':>8} {'of grid':>8} {'at top':>7}")
print("ambig = far-from-both hypotheses scoring within 2% of the worse true source")
print("-"*80)
for sep,dt in [(3.0,0.3),(6.0,0.5),(12.0,1.0)]:
    S2=S1+np.array([sep,0.,0.])
    for nc in [30]:
        for name,use in [('TT only',set()),('CIII only',{'CIII'}),('CIV+CV',{'CIV','CV'}),
                         ('ALL',{'CIII','CIV','CV'})]:
            s1,s2,amb,ambtop,ng=run(S1,5.0,S2,5.0+dt,nc,use)
            print(f"{sep:5.1f} {dt:6.1f} {nc:6d} {name:>11} | {s1:7.3f} {s2:7.3f} {amb:8d} {100*amb/ng:7.1f}% {ambtop:7d}")
        print()
