"""
F1 (gate): does a FULL WAVEFORM condition add anything beyond the minimal sufficient set?
C1: at what station count does the advantage of the constraint stack over travel time vanish?

The waveform channel is a genuine matched-filter score: normalized correlation between the
observed three-component record and the fully predicted waveform (both branches, correct
polarization, 1/d amplitude, correct wavelet).  It is precomputed by convolving each component
with the wavelet, so evaluating a candidate is O(1) per arrival.
"""
import numpy as np, statistics
rng=np.random.default_rng(202)
FC=2.0; VS=3.0; VP=1.73*VS; DT=1/(8*FC); T=np.arange(0,40,DT); NT=len(T)
SNR=6.0; TOL=0.6; AP,AS=1.0,1.5; VMIN=2.0; SRC=np.array([0.,0.,8.]); TAU=6.0
def ricker(t):
    a=(np.pi*FC*t)**2; return (1-2*a)*np.exp(-a)
WAV=ricker(np.arange(-3,3,DT))
def pol(src,r,b):
    v=r-src; n=v/np.linalg.norm(v)
    if b=='P': return n
    ez=np.array([0.,0.,1.]); sv=ez-np.dot(ez,n)*n; nn=np.linalg.norm(sv)
    return sv/nn if nn>1e-9 else np.array([1.,0.,0.])
def net(N,kind='random',rmax=45.):
    if kind=='spread': r=np.linspace(6,62,N); a=np.linspace(0,2*np.pi,N,endpoint=False)+.3
    else: a=rng.uniform(0,2*np.pi,N); r=rmax*np.sqrt(rng.uniform(.03,1,N))
    return np.stack([r*np.cos(a),r*np.sin(a),np.zeros(N)],axis=1)
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
    HIT=np.zeros((N,NT),bool); NEAR=np.zeros((N,NT)); R=np.zeros((N,NT,3))
    for i in range(N):
        for c in range(3):
            R[i,:,c]=np.convolve(D[i,c],WAV,mode='same')      # matched-filter response
        for j in range(NT):
            seg=D[i][:,max(0,j-W//2):min(NT,j+W//2)]; C=seg@seg.T
            AMP[i,j]=np.sqrt(np.trace(C)); PD[i,j]=np.linalg.eigh(C)[1][:,-1]
        p=np.sort(np.array(pk[i])); idx=np.clip(np.searchsorted(p,T),0,len(p)-1)
        cand=np.stack([p[np.clip(idx-1,0,len(p)-1)],p[idx]])
        best=cand[np.argmin(np.abs(cand-T),axis=0),np.arange(NT)]
        NEAR[i]=best; HIT[i]=np.abs(best-T)<TOL
    return AMP,PD,HIT,NEAR,R,np.linalg.norm(D)
def score(src,tau,S,DIJ,pre,use):
    AMP,PD,HIT,NEAR,R,DN=pre; N=len(S); st=np.arange(N)
    d=np.linalg.norm(S-src,axis=1); tt=np.stack([d/VP,d/VS])+tau
    idx=np.clip((tt/DT).astype(int),0,NT-1)
    if 'WF' in use:                                   # matched filter, no picks at all
        num=0.; nrm=0.
        for b,(bb,A) in enumerate((('P',AP),('S',AS))):
            for i in range(N):
                w=A/d[i]; num+= w*np.dot(pol(src,S[i],bb),R[i,idx[b,i]]); nrm+= w*w
        return max(num,0.)/ (np.sqrt(nrm)*DN + 1e-12)
    hit=HIT[st[None,:],idx]; n=hit.sum()
    if n==0: return 0.
    if 'CIII' in use:
        cl=NEAR[st[None,:],idx]
        for b in range(2):
            ii=np.where(hit[b])[0]
            if len(ii)>1 and np.any(np.abs(cl[b][ii][:,None]-cl[b][ii][None,:])>DIJ[np.ix_(ii,ii)]/VMIN+1e-9): return 0.
        both=hit[0]&hit[1]
        if np.any(cl[1][both]<=cl[0][both]): return 0.
    s=n/(2*N)
    if 'CIV' in use:
        cs=[abs(np.dot(pol(src,S[i],bb),PD[i,idx[b,i]])) for b,bb in enumerate(('P','S')) for i in range(N) if hit[b,i]]
        s*= np.mean(cs)**2 if cs else 0.
    return s
gx=np.arange(-18,18.1,4.5); gz=np.array([4.,8.,12.]); TAUS=np.arange(5.0,7.01,0.05)   # 1/10 of a period, so the matched-filter peak is resolved
GP=[np.array([x,y,z]) for x in gx for y in gx for z in gz]
FAR=np.array([np.linalg.norm(g-SRC)>5. for g in GP])
def fe_rate(S,pre,DIJ,use,):
    v=np.array([max(score(g,tt,S,DIJ,pre,use) for tt in TAUS) for g in GP])
    stv=max(score(SRC,tt,S,DIJ,pre,use) for tt in TAUS)
    return int(v[FAR].max()>=stv)
SETS=[('TT',set()),('C-IV',{'CIV'}),('MIN',{'CIII','CIV'}),('WF',{'WF'}),('MIN+WF',{'CIII','CIV','WF'})]
REPS=8
print(f"false-event rate, clutter=35, SNR={SNR}\n")
print(f"{'N sta':>6} | " + " ".join(f"{n:>7}" for n,_ in SETS) + "   | MIN-TT   WF-MIN")
print("-"*66)
res={}
for N in [3,4,6,12]:
    S=net(N); DIJ=np.linalg.norm(S[:,None,:]-S[None,:,:],axis=2)
    fe={n:0 for n,_ in SETS}
    for _ in range(REPS):
        D,pk=gen(S,35); pre=precompute(D,pk)
        for n,u in SETS: fe[n]+=fe_rate(S,pre,DIJ,u)
    r={n:fe[n]/REPS for n,_ in SETS}; res[N]=r
    print(f"{N:6d} | " + " ".join(f"{r[n]:7.2f}" for n,_ in SETS) +
          f"   | {r['MIN']-r['TT']:+7.2f} {r['WF']-r['MIN']:+8.2f}")
print("\nC1: advantage of MIN over travel time, by station count")
for N in res: print(f"   N={N:2d}: {res[N]['TT']-res[N]['MIN']:+.2f}")
print("\nF1 gate: does the waveform condition beat the minimal set by >=0.03?")
for N in res:
    d=res[N]['MIN']-res[N]['WF']
    print(f"   N={N:2d}: MIN {res[N]['MIN']:.2f}  WF {res[N]['WF']:.2f}  MIN+WF {res[N]['MIN+WF']:.2f}  -> WF gain {d:+.2f}")
