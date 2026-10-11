"""
B3: does the observability constant predict WHERE discrimination is hard?

Appendix J showed the constant is not a discriminating constraint when used as a feasibility
mask.  Its proper role is different: gamma should be small exactly where the network sees the
source poorly, and that is where clutter should be able to win.  One fixed network; the true
source is moved over a grid of locations; for each location gamma is computed in closed form and
the false-hypothesis margin is measured over replicates.
"""
import numpy as np
from scipy.stats import spearmanr
rng=np.random.default_rng(606)
FC=2.0; VS=3.0; VP=1.73*VS; DT=1/(8*FC); T=np.arange(0,45,DT); NT=len(T)
SNR=6.0; TOL=0.6; AP,AS=1.0,1.5; VMIN=2.0
# deliberately asymmetric network so that gamma varies strongly across the region
ang=np.array([0.2,0.8,1.3,1.9,2.4,3.0]); rad=np.array([10.,16.,22.,30.,38.,46.])
STA=np.stack([rad*np.cos(ang),rad*np.sin(ang),np.zeros(6)],axis=1)
DIJ=np.linalg.norm(STA[:,None,:]-STA[None,:,:],axis=2)
def ricker(t):
    a=(np.pi*FC*t)**2; return (1-2*a)*np.exp(-a)
WAV=ricker(np.arange(-3,3,DT))
def pol(src,r,b):
    v=r-src; n=v/np.linalg.norm(v)
    if b=='P': return n
    ez=np.array([0.,0.,1.]); sv=ez-np.dot(ez,n)*n; nn=np.linalg.norm(sv)
    return sv/nn if nn>1e-9 else np.array([1.,0.,0.])
def gamma_lead(src):
    M=np.zeros((4,4))
    for r in STA:
        d=np.linalg.norm(r-src); n=(r-src)/d
        for V,A in ((VP,AP),(VS,AS)):
            v=np.concatenate([-n/V,[1.]]); M+=(A/d)**2*np.outer(v,v)
    return np.linalg.eigvalsh(M)[0]
def gen(src,tau,nc):
    D=np.zeros((6,3,NT)); pk=[[] for _ in range(6)]
    d=np.linalg.norm(STA-src,axis=1)
    for i in range(6):
        for b,A,V in (('P',AP,VP),('S',AS,VS)):
            ta=d[i]/V+tau; pk[i].append(ta)
            D[i]+=SNR*A*np.outer(pol(src,STA[i],b),ricker(T-ta))/d[i]
    for _ in range(nc):
        i=rng.integers(6); tc=rng.uniform(2,40); q=rng.normal(size=3); q/=np.linalg.norm(q)
        D[i]+=SNR*.35*np.outer(q,ricker(T-tc))/25.; pk[i].append(tc)
    D+=rng.normal(scale=SNR*.02,size=D.shape); return D,pk
def precompute(D,pk):
    W=int(.6/DT); AMP=np.zeros((6,NT)); PD=np.zeros((6,NT,3))
    HIT=np.zeros((6,NT),bool); NEAR=np.zeros((6,NT))
    for i in range(6):
        for j in range(NT):
            seg=D[i][:,max(0,j-W//2):min(NT,j+W//2)]; C=seg@seg.T
            AMP[i,j]=np.sqrt(np.trace(C)); PD[i,j]=np.linalg.eigh(C)[1][:,-1]
        p=np.sort(np.array(pk[i])); idx=np.clip(np.searchsorted(p,T),0,len(p)-1)
        cand=np.stack([p[np.clip(idx-1,0,len(p)-1)],p[idx]])
        best=cand[np.argmin(np.abs(cand-T),axis=0),np.arange(NT)]
        NEAR[i]=best; HIT[i]=np.abs(best-T)<TOL
    return AMP,PD,HIT,NEAR
def score(src,tau,pre):
    AMP,PD,HIT,NEAR=pre; st=np.arange(6)
    d=np.linalg.norm(STA-src,axis=1); tt=np.stack([d/VP,d/VS])+tau
    idx=np.clip((tt/DT).astype(int),0,NT-1); hit=HIT[st[None,:],idx]
    if hit.sum()==0: return 0.
    cl=NEAR[st[None,:],idx]
    for b in range(2):
        ii=np.where(hit[b])[0]
        if len(ii)>1 and np.any(np.abs(cl[b][ii][:,None]-cl[b][ii][None,:])>DIJ[np.ix_(ii,ii)]/VMIN+1e-9): return 0.
    both=hit[0]&hit[1]
    if np.any(cl[1][both]<=cl[0][both]): return 0.
    s=hit.sum()/12.
    cs=[abs(np.dot(pol(src,STA[i],bb),PD[i,idx[b,i]])) for b,bb in enumerate(('P','S')) for i in range(6) if hit[b,i]]
    return s*(np.mean(cs)**2 if cs else 0.)
gx=np.arange(-30,30.1,5.0); gz=np.array([5.,9.,13.])
GP=[np.array([x,y,z]) for x in gx for y in gx for z in gz]
TAUS=np.arange(4.5,7.51,0.1)
locs=[np.array([x,y,9.]) for x in np.arange(-24,24.1,8.) for y in np.arange(-24,24.1,8.)]
print(f"one fixed asymmetric network; true source moved over {len(locs)} locations, 5 realizations each\n")
G=[];M=[];FE=[]
for L0 in locs:
    g=gamma_lead(L0); mg=[]; fe=0
    for _ in range(5):
        D,pk=gen(L0,6.0,30); pre=precompute(D,pk)
        v=np.array([max(score(gp,tt,pre) for tt in TAUS) for gp in GP])
        stv=max(score(L0,tt,pre) for tt in TAUS)
        far=np.array([np.linalg.norm(gp-L0)>5. for gp in GP])
        bf=v[far].max(); mg.append((stv-bf)/max(v.std(),1e-9)); fe+= int(bf>=stv)
    G.append(g); M.append(float(np.median(mg))); FE.append(fe/5)
G=np.array(G);M=np.array(M);FE=np.array(FE)
print(f"gamma range: {G.min():.3e} to {G.max():.3e}   ({G.max()/max(G.min(),1e-30):.0f}x)")
print(f"margin range: {M.min():.2f} to {M.max():.2f}")
r1=spearmanr(G,M).statistic; r2=spearmanr(G,-FE).statistic
print(f"\nSpearman(gamma, local discrimination margin)  = {r1:+.3f}    <- B3 requires >= 0.5")
print(f"Spearman(gamma, local success = 1-false rate) = {r2:+.3f}")
q=np.argsort(G)
print(f"\nlowest-gamma quartile : median margin {np.median(M[q[:len(q)//4]]):.2f}, false-event rate {np.mean(FE[q[:len(q)//4]]):.2f}")
print(f"highest-gamma quartile: median margin {np.median(M[q[-len(q)//4:]]):.2f}, false-event rate {np.mean(FE[q[-len(q)//4:]]):.2f}")
