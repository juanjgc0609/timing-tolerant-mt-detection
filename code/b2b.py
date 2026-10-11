"""B2: does a network designed by maximizing the observability constant beat one designed by
minimizing azimuthal gap, at equal station count?  Generalized to arbitrary N."""
import numpy as np, statistics
rng=np.random.default_rng(77)
FC=2.0; VS=3.0; VP=1.73*VS; DT=1/(8*FC); T=np.arange(0,40,DT); NT=len(T)
SNR=4.0; TOL=0.6; AP,AS=1.0,1.5; VMIN=2.0; SRC=np.array([0.,0.,8.]); TAU=6.0
def ricker(t):
    a=(np.pi*FC*t)**2; return (1-2*a)*np.exp(-a)
def pol(src,r,b):
    v=r-src; n=v/np.linalg.norm(v)
    if b=='P': return n
    ez=np.array([0.,0.,1.]); sv=ez-np.dot(ez,n)*n; nn=np.linalg.norm(sv)
    return sv/nn if nn>1e-9 else np.array([1.,0.,0.])
def gen(S,nc):
    N=len(S); D=np.zeros((N,3,NT)); pk=[[] for _ in range(N)]
    d=np.linalg.norm(S-SRC,axis=1)
    for i in range(N):
        for b,A,V in (('P',AP,VP),('S',AS,VS)):
            ta=d[i]/V+TAU; pk[i].append(ta)
            D[i]+=SNR*A*np.outer(pol(SRC,S[i],b),ricker(T-ta))/d[i]
    for _ in range(nc):
        i=rng.integers(N); tc=rng.uniform(2,36); q=rng.normal(size=3); q/=np.linalg.norm(q)
        D[i]+=SNR*.35*np.outer(q,ricker(T-tc))/np.median(np.linalg.norm(S,axis=1)); pk[i].append(tc)
    D+=rng.normal(scale=SNR*.02,size=D.shape); return D,pk
def precompute(D,pk):
    N=len(D); W=int(.6/DT); AMP=np.zeros((N,NT)); PD=np.zeros((N,NT,3))
    HIT=np.zeros((N,NT),bool); NEAR=np.zeros((N,NT))
    for i in range(N):
        for j in range(NT):
            seg=D[i][:,max(0,j-W//2):min(NT,j+W//2)]; C=seg@seg.T
            AMP[i,j]=np.sqrt(np.trace(C)); PD[i,j]=np.linalg.eigh(C)[1][:,-1]
        p=np.sort(np.array(pk[i])); idx=np.clip(np.searchsorted(p,T),0,len(p)-1)
        cand=np.stack([p[np.clip(idx-1,0,len(p)-1)],p[idx]])
        best=cand[np.argmin(np.abs(cand-T),axis=0),np.arange(NT)]
        NEAR[i]=best; HIT[i]=np.abs(best-T)<TOL
    return AMP,PD,HIT,NEAR
def score(src,tau,S,DIJ,pre):
    AMP,PD,HIT,NEAR=pre; N=len(S); st=np.arange(N)
    d=np.linalg.norm(S-src,axis=1); tt=np.stack([d/VP,d/VS])+tau
    idx=np.clip((tt/DT).astype(int),0,NT-1); hit=HIT[st[None,:],idx]
    n=hit.sum()
    if n==0: return 0.
    cl=NEAR[st[None,:],idx]
    for b in range(2):
        ii=np.where(hit[b])[0]
        if len(ii)>1:
            if np.any(np.abs(cl[b][ii][:,None]-cl[b][ii][None,:])>DIJ[np.ix_(ii,ii)]/VMIN+1e-9): return 0.
    both=hit[0]&hit[1]
    if np.any(cl[1][both]<=cl[0][both]): return 0.
    s=n/(2*N)
    cs=[abs(np.dot(pol(src,S[i],bb),PD[i,idx[b,i]])) for b,bb in enumerate(('P','S')) for i in range(N) if hit[b,i]]
    return s*(np.mean(cs)**2 if cs else 0.)
def gamma_lead(S):
    M=np.zeros((4,4))
    for r in S:
        d=np.linalg.norm(r-SRC); n=(r-SRC)/d
        for V,A in ((VP,AP),(VS,AS)):
            v=np.concatenate([-n/V,[1.]]); M+=(A/d)**2*np.outer(v,v)
    return np.linalg.eigvalsh(M)[0]
def azgap(S):
    az=np.sort(np.degrees(np.arctan2(S[:,1]-SRC[1],S[:,0]-SRC[0]))%360)
    return np.diff(np.concatenate([az,[az[0]+360]])).max()
def rnet(N,rmax=45.):
    a=rng.uniform(0,2*np.pi,N); r=rmax*np.sqrt(rng.uniform(.02,1,N))
    return np.stack([r*np.cos(a),r*np.sin(a),np.zeros(N)],axis=1)

gx=np.arange(-18,18.1,6.0); gz=np.array([4.,8.,12.]); TAUS=np.arange(4.5,7.51,0.25)
GP=[np.array([x,y,z]) for x in gx for y in gx for z in gz]
FAR=np.array([np.linalg.norm(g-SRC)>5. for g in GP])
from scipy.optimize import minimize
def cont_err(S,DJ,pre):
    v=[max(score(gp,tt,S,DJ,pre) for tt in TAUS) for gp in GP]
    g0=GP[int(np.argmax(v))]
    t0=max(TAUS,key=lambda tt: score(g0,tt,S,DJ,pre))
    f=lambda p: -score(np.array([p[0],p[1],min(20.,max(1.,p[2]))]),p[3],S,DJ,pre)
    r=minimize(f,np.array([g0[0],g0[1],g0[2],t0]),method='Nelder-Mead',
               options=dict(maxiter=300,xatol=.05,fatol=1e-4))
    return float(np.hypot(r.x[0],r.x[1]))
print(f"{'N':>3} | {'design criterion':>22} {'gamma':>9} {'azgap':>7} | {'median err km':>14} {'false-event':>12}")
print("-"*76)
for N in [4,5,6]:
    pool=[rnet(N) for _ in range(3000)]
    g=np.array([gamma_lead(s) for s in pool]); a=np.array([azgap(s) for s in pool])
    picks={'max observability':pool[int(np.argmax(g))],
           'min azimuthal gap':pool[int(np.argmin(a))],
           'median-gamma random':pool[int(np.argsort(g)[len(g)//2])]}
    out={}
    for name,S in picks.items():
        DJ=np.linalg.norm(S[:,None,:]-S[None,:,:],axis=2); errs=[]; fe=0; R=10
        for _ in range(R):
            D,pk=gen(S,45); pre=precompute(D,pk)
            errs.append(cont_err(S,DJ,pre))
            v=np.array([max(score(gp,tt,S,DJ,pre) for tt in TAUS) for gp in GP])
            stv=max(score(SRC,tt,S,DJ,pre) for tt in TAUS)
            fe+= int(v[FAR].max()>=stv)
        m=statistics.median(errs); out[name]=m
        print(f"{N:3d} | {name:>22} {gamma_lead(S):9.2e} {azgap(S):7.1f} | {m:14.2f} {fe/R:12.2f}")
    imp=100*(out['min azimuthal gap']-out['max observability'])/max(out['min azimuthal gap'],1e-9)
    print(f"    -> observability design vs azimuthal-gap design: {imp:+.0f}% error reduction "
          f"(B2 requires >=20%)\n")
