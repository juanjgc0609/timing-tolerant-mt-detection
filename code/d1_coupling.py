"""
D1: coupling between detection and association.

Instead of simulating a detector badly, the picker MISS RATE is the independent variable: each
true arrival is deleted from the pick list with probability m, while clutter picks remain (a
picker that misses buried arrivals still triggers on impulsive noise).  The waveform channel
uses no picks at all.

Endpoint: recovery rate, the fraction of realizations whose best hypothesis lies within 5 km of
the truth.  This captures both failure modes at once - unable to form the event, and forming a
false one.
Deliverable: the miss rate at which the pick-free channel overtakes the best pick-based channel.
"""
import numpy as np
rng=np.random.default_rng(505)
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
def net(N):
    a=np.linspace(0,2*np.pi,N,endpoint=False)+.4; r=np.linspace(8,50,N)
    return np.stack([r*np.cos(a),r*np.sin(a),np.zeros(N)],axis=1)
def gen(S,nc,miss):
    N=len(S); D=np.zeros((N,3,NT)); pk=[[] for _ in range(N)]
    d=np.linalg.norm(S-SRC,axis=1)
    for i in range(N):
        for b,A,V in (('P',AP,VP),('S',AS,VS)):
            ta=d[i]/V+TAU
            D[i]+=SNR*A*np.outer(pol(SRC,S[i],b),ricker(T-ta))/d[i]   # waveform always present
            if rng.random()>miss: pk[i].append(ta)                     # pick may be missed
    for _ in range(nc):
        i=rng.integers(N); tc=rng.uniform(2,36); q=rng.normal(size=3); q/=np.linalg.norm(q)
        D[i]+=SNR*.35*np.outer(q,ricker(T-tc))/np.median(np.linalg.norm(S,axis=1)); pk[i].append(tc)
    for i in range(N):
        if not pk[i]: pk[i].append(-99.)
    D+=rng.normal(scale=SNR*.02,size=D.shape); return D,pk
def precompute(D,pk):
    N=len(D); W=int(.6/DT); AMP=np.zeros((N,NT)); PD=np.zeros((N,NT,3))
    HIT=np.zeros((N,NT),bool); NEAR=np.zeros((N,NT)); R=np.zeros((N,NT,3))
    for i in range(N):
        for c in range(3): R[i,:,c]=np.convolve(D[i,c],WAV,mode='same')
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
    if 'WF' in use:
        num=0.;nrm=0.
        for b,(bb,A) in enumerate((('P',AP),('S',AS))):
            for i in range(N):
                w=A/d[i]; num+=w*np.dot(pol(src,S[i],bb),R[i,idx[b,i]]); nrm+=w*w
        return max(num,0.)/(np.sqrt(nrm)*DN+1e-12)
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
gx=np.arange(-18,18.1,4.5); gz=np.array([4.,8.,12.]); TAUS=np.arange(5.,7.01,0.05)
GP=[np.array([x,y,z]) for x in gx for y in gx for z in gz]
def recovered(S,DIJ,pre,use):
    v=np.array([max(score(g,tt,S,DIJ,pre,use) for tt in TAUS) for g in GP])
    j=int(np.argmax(v))
    return int(np.linalg.norm(GP[j]-SRC)<5.)
SETS=[('TT',set()),('C-IV',{'CIV'}),('MIN',{'CIII','CIV'}),('WF pick-free',{'WF'})]
S=net(6); DIJ=np.linalg.norm(S[:,None,:]-S[None,:,:],axis=2); REPS=12
print(f"recovery rate (best hypothesis within 5 km of truth), {REPS} realizations, 6 stations\n")
print(f"{'miss rate':>10} | " + " ".join(f"{n:>13}" for n,_ in SETS) + "  | WF - best pick-based")
print("-"*84)
for m in [0.0,0.1,0.2,0.3,0.4,0.6,0.8]:
    rec={n:0 for n,_ in SETS}
    for _ in range(REPS):
        D,pk=gen(S,30,m); pre=precompute(D,pk)
        for n,u in SETS: rec[n]+=recovered(S,DIJ,pre,u)
    r={n:rec[n]/REPS for n,_ in SETS}
    bestpick=max(r['TT'],r['C-IV'],r['MIN'])
    print(f"{m:10.1f} | " + " ".join(f"{r[n]:13.2f}" for n,_ in SETS) +
          f"  | {r['WF pick-free']-bestpick:+.2f}")
