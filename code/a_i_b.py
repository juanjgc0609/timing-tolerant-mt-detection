"""
A-i: an associator that uses NO velocity model.

Instead of predicting polarization from a candidate, use the MEASURED back-azimuths and intersect
them.  For a true source the rays from the stations meet at a point; for clutter they do not.
Feasibility comes from the eikonal Lipschitz bound of C-III, which needs only V_min.

  score(x) = mean over detected arrivals of cos^2(angle between measured back-azimuth
             and the station-to-candidate direction),  subject to the C-III admissibility test.

No V_P, no V_S, no traveltime prediction: only a lower bound V_min and observed directions.

Decisive test: behaviour under STRUCTURAL velocity-model error (data generated in a gradient
medium, inference under a homogeneous model), which is exactly where A3 was refuted for the
prediction-based channels.
"""
import numpy as np, statistics
rng=np.random.default_rng(1111)
FC=2.0; DT=1/(8*FC); T=np.arange(0,44,DT); NT=len(T)
VS0,VPR=2.6,1.73; SNR=7.0; TOL=0.7; AP,AS=1.0,1.5; VMIN=2.0
SRC=np.array([0.,0.,8.]); TAU=6.0
a=np.linspace(0,2*np.pi,8,endpoint=False)+.4; rr=np.linspace(8,50,8)
STA=np.stack([rr*np.cos(a),rr*np.sin(a),np.zeros(8)],axis=1)
DIJ=np.linalg.norm(STA[:,None,:]-STA[None,:,:],axis=2)
def ricker(t):
    x=(np.pi*FC*t)**2; return (1-2*x)*np.exp(-x)
def vel(z,k,b): return (VS0+k*z)*(VPR if b=='P' else 1.)
def ttrue(x1,x2,k,b):
    v1=vel(x1[2],k,b); v2=vel(x2[2],k,b); D=np.linalg.norm(x2-x1)
    return D/v1 if k<1e-9 else np.arccosh(1+k**2*D**2/(2*v1*v2))/k
def dir_rec(x1,x2,k,b,h=1e-4):
    g=np.zeros(3)
    for j in range(3):
        e=np.zeros(3); e[j]=h
        g[j]=(ttrue(x1,x2+e,k,b)-ttrue(x1,x2-e,k,b))/(2*h)
    return g/np.linalg.norm(g)
def polv(n,b):
    if b=='P': return n
    ez=np.array([0.,0.,1.]); sv=ez-np.dot(ez,n)*n; nn=np.linalg.norm(sv)
    return sv/nn if nn>1e-9 else np.array([1.,0.,0.])
def gen(k,nc=25):
    D=np.zeros((8,3,NT)); pk=[[] for _ in range(8)]
    for i in range(8):
        d=np.linalg.norm(STA[i]-SRC)
        for b,A in (('P',AP),('S',AS)):
            ta=ttrue(SRC,STA[i],k,b)+TAU; pk[i].append(ta)
            D[i]+=SNR*A*np.outer(polv(dir_rec(SRC,STA[i],k,b),b),ricker(T-ta))/d
    for _ in range(nc):
        i=rng.integers(8); tc=rng.uniform(2,40); q=rng.normal(size=3); q/=np.linalg.norm(q)
        D[i]+=SNR*.32*np.outer(q,ricker(T-tc))/25.; pk[i].append(tc)
    D+=rng.normal(scale=SNR*.02,size=D.shape); return D,pk
def pre_(D,pk):
    W=int(.6/DT); PD=np.zeros((8,NT,3)); HIT=np.zeros((8,NT),bool); NEAR=np.zeros((8,NT))
    for i in range(8):
        for j in range(NT):
            seg=D[i][:,max(0,j-W//2):min(NT,j+W//2)]; C=seg@seg.T
            PD[i,j]=np.linalg.eigh(C)[1][:,-1]
        p=np.sort(np.array(pk[i])); idx=np.clip(np.searchsorted(p,T),0,len(p)-1)
        cand=np.stack([p[np.clip(idx-1,0,len(p)-1)],p[idx]])
        best=cand[np.argmin(np.abs(cand-T),axis=0),np.arange(NT)]
        NEAR[i]=best; HIT[i]=np.abs(best-T)<TOL
    return PD,HIT,NEAR
VSH=VS0; VPH=VS0*VPR
def score_pred(src,tau,pre,use):
    """prediction-based channels: need the velocity model"""
    PD,HIT,NEAR=pre; st=np.arange(8)
    d=np.linalg.norm(STA-src,axis=1); tt=np.stack([d/VPH,d/VSH])+tau
    idx=np.clip((tt/DT).astype(int),0,NT-1); hit=HIT[st[None,:],idx]
    if hit.sum()==0: return 0.
    s=hit.sum()/16.
    if 'CIV' in use:
        cs=[abs(np.dot(polv((STA[i]-src)/d[i],bb),PD[i,idx[b,i]])) for b,bb in enumerate(('P','S')) for i in range(8) if hit[b,i]]
        s*= np.mean(cs)**2 if cs else 0.
    return s
def picks_and_dirs(pk,pre):
    """observed arrivals with their measured principal direction - no model used"""
    PD,HIT,NEAR=pre; out=[]
    for i in range(8):
        for tt in pk[i]:
            j=int(np.clip(tt/DT,0,NT-1)); out.append((i,tt,PD[i,j]))
    return out
COS_TOL=np.cos(np.radians(20.))
def score_azim(srcxy,obs,subset=None):
    """A-i: robust intersection of measured back-azimuths, with C-III admissibility.
       Only arrivals whose measured azimuth points at the candidate within a tolerance are
       counted (consensus, clutter-robust); the consistent set must also satisfy the
       Lipschitz test.  No velocity model at all."""
    keep=[]
    for q,(i,tt,pdv) in enumerate(obs):
        v=np.array([STA[i][0]-srcxy[0], STA[i][1]-srcxy[1]]); n=np.linalg.norm(v)
        if n<1e-6: continue
        w=pdv[:2]; nw=np.linalg.norm(w)
        if nw<1e-6: continue
        c=abs(np.dot(v/n,w/nw))
        if c>=COS_TOL: keep.append((q,c,i))
    if len(keep)<3: return 0.
    if not lipschitz_ok(obs,[q for q,_,_ in keep]): return 0.
    nsta=len({i for _,_,i in keep})
    return (nsta/8.)*float(np.mean([c for _,c,_ in keep])**2)
def lipschitz_ok(obs,subset):
    ii=[obs[q][0] for q in subset]; tt=[obs[q][1] for q in subset]
    for a_ in range(len(ii)):
        for b_ in range(a_+1,len(ii)):
            if ii[a_]==ii[b_]: continue
            if abs(tt[a_]-tt[b_])>DIJ[ii[a_],ii[b_]]/VMIN+1e-9: return False
    return True
gx=np.arange(-20,20.1,2.5)
GXY=[(x,y) for x in gx for y in gx]
GP=[np.array([x,y,z]) for x in gx for y in gx for z in [5.,8.,11.]]
TAUS=np.arange(5.,8.01,0.15)
def run(k,R=8):
    e_tt=[];e_civ=[];e_az=[]
    for _ in range(R):
        D,pk=gen(k); pre=pre_(D,pk); obs=picks_and_dirs(pk,pre)
        v=np.array([max(score_pred(g,tt,pre,set()) for tt in TAUS) for g in GP])
        j=int(np.argmax(v)); e_tt.append(np.hypot(GP[j][0],GP[j][1]))
        v=np.array([max(score_pred(g,tt,pre,{'CIV'}) for tt in TAUS) for g in GP])
        j=int(np.argmax(v)); e_civ.append(np.hypot(GP[j][0],GP[j][1]))
        # A-i: greedy admissible subset then azimuth intersection over the x-y grid only
        best=(-1,None)
        for xy in GXY:
            sc=score_azim(xy,obs)
            if sc>best[0]: best=(sc,xy)
        e_az.append(np.hypot(best[1][0],best[1][1]))
    return statistics.median(e_tt),statistics.median(e_civ),statistics.median(e_az)
print("epicentral error [km], data generated in a gradient medium, inference homogeneous\n")
print(f"{'k':>5} {'model err':>10} | {'TT':>8} {'C-IV':>8} {'A-i azimuth (no model)':>24}")
print("-"*58)
for k in [0.0,0.06,0.12,0.20]:
    e=[]
    for i in range(8):
        for b,Vh in (('P',VPH),('S',VSH)):
            tt_=ttrue(SRC,STA[i],k,b); e.append((np.linalg.norm(STA[i]-SRC)/Vh-tt_)/tt_)
    err=100*np.sqrt(np.mean(np.array(e)**2))
    a1,a2,a3=run(k)
    print(f"{k:5.2f} {err:9.1f}% | {a1:8.2f} {a2:8.2f} {a3:24.2f}")
