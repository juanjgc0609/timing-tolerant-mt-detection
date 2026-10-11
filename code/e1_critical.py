"""
E-i: at what polarization coherence, and at what separation, does the polarization channel break?

The dangerous clutter is not random and not aligned with the true source (which reinforces it).
It is clutter coherent with a NEARBY decoy - an aftershock, or the coda of a neighbour.  Sweep the
decoy separation and its support level, and find where the true source stops being recovered.

Deliverable: the critical decoy separation for the polarization channel and for travel time, so
that a measurement of real clutter coherence can be compared against a known threshold.
"""
import numpy as np
rng=np.random.default_rng(909)
FC=2.0; VS=3.0; VP=1.73*VS; DT=1/(8*FC); T=np.arange(0,40,DT); NT=len(T)
SNR=7.0; TOL=0.6; AP,AS=1.0,1.5; VMIN=2.0; SRC=np.array([0.,0.,8.]); TAU=6.0
a=np.linspace(0,2*np.pi,8,endpoint=False)+.4; r=np.linspace(8,50,8)
STA=np.stack([r*np.cos(a),r*np.sin(a),np.zeros(8)],axis=1)
DIJ=np.linalg.norm(STA[:,None,:]-STA[None,:,:],axis=2)
def ricker(t):
    x=(np.pi*FC*t)**2; return (1-2*x)*np.exp(-x)
def pol(src,rr,b):
    v=rr-src; n=v/np.linalg.norm(v)
    if b=='P': return n
    ez=np.array([0.,0.,1.]); sv=ez-np.dot(ez,n)*n; nn=np.linalg.norm(sv)
    return sv/nn if nn>1e-9 else np.array([1.,0.,0.])
def gen(dec,ndec,dtau=1.4,nrand=8):
    D=np.zeros((8,3,NT)); pk=[[] for _ in range(8)]
    d=np.linalg.norm(STA-SRC,axis=1)
    for i in range(8):
        for b,A,V in (('P',AP,VP),('S',AS,VS)):
            ta=d[i]/V+TAU; pk[i].append(ta)
            D[i]+=SNR*A*np.outer(pol(SRC,STA[i],b),ricker(T-ta))/d[i]
    dd=np.linalg.norm(STA-dec,axis=1)
    slots=[(i,b) for i in range(8) for b in ('P','S')]; rng.shuffle(slots)
    for i,b in slots[:ndec]:
        A=AP if b=='P' else AS; V=VP if b=='P' else VS
        ta=dd[i]/V+TAU+dtau; pk[i].append(ta)
        D[i]+=SNR*A*np.outer(pol(dec,STA[i],b),ricker(T-ta))/dd[i]
    for _ in range(nrand):
        i=rng.integers(8); tc=rng.uniform(2,36); q=rng.normal(size=3); q/=np.linalg.norm(q)
        D[i]+=SNR*.3*np.outer(q,ricker(T-tc))/25.; pk[i].append(tc)
    D+=rng.normal(scale=SNR*.02,size=D.shape); return D,pk
def pre_(D,pk):
    W=int(.6/DT); AMP=np.zeros((8,NT)); PD=np.zeros((8,NT,3))
    HIT=np.zeros((8,NT),bool); NEAR=np.zeros((8,NT))
    for i in range(8):
        for j in range(NT):
            seg=D[i][:,max(0,j-W//2):min(NT,j+W//2)]; C=seg@seg.T
            AMP[i,j]=np.sqrt(np.trace(C)); PD[i,j]=np.linalg.eigh(C)[1][:,-1]
        p=np.sort(np.array(pk[i])); idx=np.clip(np.searchsorted(p,T),0,len(p)-1)
        cand=np.stack([p[np.clip(idx-1,0,len(p)-1)],p[idx]])
        best=cand[np.argmin(np.abs(cand-T),axis=0),np.arange(NT)]
        NEAR[i]=best; HIT[i]=np.abs(best-T)<TOL
    return AMP,PD,HIT,NEAR
def score(src,tau,pre,use):
    AMP,PD,HIT,NEAR=pre; st=np.arange(8)
    d=np.linalg.norm(STA-src,axis=1); tt=np.stack([d/VP,d/VS])+tau
    idx=np.clip((tt/DT).astype(int),0,NT-1); hit=HIT[st[None,:],idx]
    if hit.sum()==0: return 0.
    if 'CIII' in use:
        cl=NEAR[st[None,:],idx]
        for b in range(2):
            ii=np.where(hit[b])[0]
            if len(ii)>1 and np.any(np.abs(cl[b][ii][:,None]-cl[b][ii][None,:])>DIJ[np.ix_(ii,ii)]/VMIN+1e-9): return 0.
    s=hit.sum()/16.
    if 'CIV' in use:
        cs=[abs(np.dot(pol(src,STA[i],bb),PD[i,idx[b,i]])) for b,bb in enumerate(('P','S')) for i in range(8) if hit[b,i]]
        s*= np.mean(cs)**2 if cs else 0.
    return s
gx=np.arange(-20,20.1,2.5); gz=np.array([5.,8.,11.]); TAUS=np.arange(5.,8.51,0.15)
GP=[np.array([x,y,z]) for x in gx for y in gx for z in gz]
def outcome(dec,ndec,use,R=6):
    ok=0; merged=0
    for _ in range(R):
        D,pk=gen(dec,ndec); pre=pre_(D,pk)
        v=np.array([max(score(g,tt,pre,use) for tt in TAUS) for g in GP])
        j=int(np.argmax(v)); best=GP[j]
        if np.linalg.norm(best-SRC)<3.: ok+=1
        elif np.linalg.norm(best-dec)>3.: merged+=1     # neither: a merged/spurious solution
    return ok/R, merged/R
print("true source recovery rate as a function of decoy separation and decoy support\n")
print(f"{'sep km':>7} {'support':>8} | {'TT ok':>7} {'TT spur':>8} | {'C-IV ok':>8} {'C-IV spur':>10}")
print("-"*60)
for sep in [2.,4.,7.,12.]:
    dec=SRC+np.array([sep,0.,0.])
    for nd in [6,12]:
        o1,m1=outcome(dec,nd,{'CIII'}); o2,m2=outcome(dec,nd,{'CIII','CIV'})
        print(f"{sep:7.1f} {nd:6d}/16 | {o1:7.2f} {m1:8.2f} | {o2:8.2f} {m2:10.2f}")
    print()
