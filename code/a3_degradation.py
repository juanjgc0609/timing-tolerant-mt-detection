"""
A3: graceful degradation under velocity-model error.

Honest design: data are generated in a medium with a velocity GRADIENT; inference uses a
HOMOGENEOUS model.  The error is therefore structural - arrival times AND ray directions are
both wrong - not a trivial rescaling that would leave polarization untouched.

Channels compared as the gradient (hence the model error) increases:
  TT      travel-time matching only            (needs the velocity model)
  C-III   Lipschitz + P-S ordering             (needs only a lower bound Vmin)
  C-IV    polarization agreement               (needs ray directions)
  C-V     implied-moment consistency           (needs amplitudes)
  MIN     C-III + C-IV   (candidate minimal sufficient set of A2)
  ALL     everything
Metric: discrimination margin, true source vs best hypothesis >5 km away, in sd units.
"""
import numpy as np
rng=np.random.default_rng(17)
FC=2.0; DT=1/(8*FC); T=np.arange(0,44,DT); NT=len(T)
VS0,VPR=2.6,1.73; SNR=8.0; TOL=0.6; AP,AS=1.0,1.5; VMIN=2.0
SRC=np.array([0.,0.,8.]); TAU=6.0
th=np.linspace(0,2*np.pi,8,endpoint=False)+.3
RR=np.array([6.,11.,18.,27.,36.,45.,54.,62.])
STA=np.stack([RR*np.cos(th),RR*np.sin(th),np.zeros(8)],axis=1)
DIJ=np.linalg.norm(STA[:,None,:]-STA[None,:,:],axis=2)
def ricker(t):
    a=(np.pi*FC*t)**2; return (1-2*a)*np.exp(-a)
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
def pol_from_dir(n,b):
    if b=='P': return n
    ez=np.array([0.,0.,1.]); sv=ez-np.dot(ez,n)*n; nn=np.linalg.norm(sv)
    return sv/nn if nn>1e-9 else np.array([1.,0.,0.])

def gen(k,nc=25):
    D=np.zeros((8,3,NT)); picks=[[] for _ in range(8)]
    for i in range(8):
        d=np.linalg.norm(STA[i]-SRC)
        for b,A in (('P',AP),('S',AS)):
            ta=ttrue(SRC,STA[i],k,b)+TAU; picks[i].append(ta)
            p=pol_from_dir(dir_rec(SRC,STA[i],k,b),b)
            D[i]+=SNR*A*np.outer(p,ricker(T-ta))/d
    for _ in range(nc):
        i=rng.integers(8); tc=rng.uniform(2,40); q=rng.normal(size=3); q/=np.linalg.norm(q)
        D[i]+=SNR*0.35*np.outer(q,ricker(T-tc))/np.median(RR); picks[i].append(tc)
    D+=rng.normal(scale=SNR*0.02,size=D.shape)
    return D,picks
def feats(D):
    W=int(0.6/DT); E=np.zeros((8,NT)); PD=np.zeros((8,NT,3))
    for i in range(8):
        for j in range(NT):
            seg=D[i][:,max(0,j-W//2):min(NT,j+W//2)]; C=seg@seg.T
            E[i,j]=np.trace(C); PD[i,j]=np.linalg.eigh(C)[1][:,-1]
    return np.sqrt(E),PD

# --- inference model is HOMOGENEOUS with the surface velocity ---
VSH=VS0; VPH=VS0*VPR
gx=np.arange(-21,21.1,3.5); gz=np.array([4.,8.,12.]); taus=np.arange(4.,9.01,0.15)
GRID=[np.array([x,y,z]) for x in gx for y in gx for z in gz]

def score(src,AMP,PD,picks,use):
    d=np.linalg.norm(STA-src,axis=1); best=-1
    for tau in taus:
        tt=np.stack([d/VPH,d/VSH])+tau; idx=np.clip((tt/DT).astype(int),0,NT-1)
        hit=np.array([[any(abs(tt[b,i]-pk)<TOL for pk in picks[i]) for i in range(8)] for b in range(2)])
        if hit.sum()==0: continue
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
            cs=[abs(np.dot(pol_from_dir((STA[i]-src)/d[i],b),PD[i,idx[bi,i]]))
                for bi,b in enumerate(('P','S')) for i in range(8) if hit[bi,i]]
            s*= np.mean(cs)**2 if cs else 0.
        if 'CV' in use:
            m=[AMP[i,idx[bi,i]]*d[i]/A for bi,(b,A) in enumerate((('P',AP),('S',AS))) for i in range(8) if hit[bi,i]]
            if len(m)>=2:
                m=np.array(m); s*=1./(1.+m.std()/max(m.mean(),1e-9))
            else: s*=0.5
        best=max(best,s)
    return best
def margin(D,picks,use):
    AMP,PD=feats(D)
    v=np.array([score(g,AMP,PD,picks,use) for g in GRID])
    st=score(SRC,AMP,PD,picks,use)
    far=np.array([np.linalg.norm(g-SRC)>5. for g in GRID]); ok=v>-.5
    bf=v[far&ok].max() if (far&ok).sum() else 0.
    return (st-bf)/max(v[ok].std(),1e-9)

sets=[('TT',set()),('C-III',{'CIII'}),('C-IV',{'CIV'}),('C-V',{'CV'}),
      ('MIN=III+IV',{'CIII','CIV'}),('ALL',{'CIII','CIV','CV'})]
print(f"{'k':>5} {'model err':>10} | " + " ".join(f"{n:>10}" for n,_ in sets))
print("-"*82)
base={}
for k in [0.0,0.04,0.08,0.14,0.20]:
    # RMS relative traveltime error of the homogeneous model vs the true gradient medium
    e=[]
    for i in range(8):
        for b,Vh in (('P',VPH),('S',VSH)):
            tt_=ttrue(SRC,STA[i],k,b); e.append((np.linalg.norm(STA[i]-SRC)/Vh-tt_)/tt_)
    err=100*np.sqrt(np.mean(np.array(e)**2))
    D,picks=gen(k)
    row=[]
    for n,u in sets:
        m=margin(D,picks,u); row.append(m)
        if k==0.0: base[n]=m
    print(f"{k:5.2f} {err:9.1f}% | " + " ".join(f"{x:10.2f}" for x in row))
print("\nretention (margin / margin at zero model error)")
print(f"{'k':>5} {'model err':>10} | " + " ".join(f"{n:>10}" for n,_ in sets))
for k in [0.04,0.08,0.14,0.20]:
    e=[]
    for i in range(8):
        for b,Vh in (('P',VPH),('S',VSH)):
            tt_=ttrue(SRC,STA[i],k,b); e.append((np.linalg.norm(STA[i]-SRC)/Vh-tt_)/tt_)
    err=100*np.sqrt(np.mean(np.array(e)**2))
    D,picks=gen(k)
    print(f"{k:5.2f} {err:9.1f}% | " + " ".join(f"{margin(D,picks,u)/max(base[n],1e-9):10.2f}" for n,u in sets))
