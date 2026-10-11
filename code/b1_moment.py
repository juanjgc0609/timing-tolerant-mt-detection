"""
B1: physicality of the recovered moment tensor as a discriminant.

For a candidate location the representation is LINEAR in M, so M comes from closed-form least
squares - no search over origin time, no picking.  A tectonic earthquake is deviatoric: the
isotropic part of M is near zero, and the residual of the linear fit is small.  A false
hypothesis assembling clutter has no reason to produce either.

Two scalar tests per candidate, both closed form:
   iso   = |trace(M)| / ||M||          (should be ~0 for a real earthquake)
   resid = ||d - G M|| / ||d||         (should be small)
Compared against the polarization channel that carries the discrimination in the document.
"""
import numpy as np, statistics
rng=np.random.default_rng(808)
FC=2.0; VS=3.0; VP=1.73*VS; DT=1/(8*FC); T=np.arange(0,40,DT); NT=len(T)
SNR=7.0; TAU=6.0
a=np.linspace(0,2*np.pi,10,endpoint=False)+.3
r=np.array([9.,13.,18.,23.,29.,35.,41.,47.,52.,58.])
STA=np.stack([r*np.cos(a),r*np.sin(a),np.zeros(10)],axis=1); NS=len(STA)
def ricker(t):
    x=(np.pi*FC*t)**2; return (1-2*x)*np.exp(-x)
BASIS=[np.array([[1.,0,0],[0,-1.,0],[0,0,0]]), np.array([[0,1.,0],[1.,0,0],[0,0,0]]),
       np.array([[0,0,1.],[0,0,0],[1.,0,0]]), np.array([[0,0,0],[0,0,1.],[0,1.,0]]),
       np.array([[0,0,0],[0,1.,0],[0,0,-1.]]), np.eye(3)/np.sqrt(3)]   # 5 deviatoric + isotropic

def green_cols(src,tau):
    """columns of the linear operator: response of each moment basis element, 3C x stations"""
    cols=[]
    d=np.linalg.norm(STA-src,axis=1)
    for B in BASIS:
        out=np.zeros((NS,3,NT))
        for i in range(NS):
            n=(STA[i]-src)/d[i]
            ez=np.array([0.,0.,1.]); sv=ez-np.dot(ez,n)*n; nn=np.linalg.norm(sv)
            sv=sv/nn if nn>1e-9 else np.array([1.,0.,0.])
            # far-field radiation: P amplitude ~ n.B.n along n ; S ~ transverse projection
            aP=float(n@B@n); aS=float(sv@B@n)
            out[i]=np.outer(n,aP*ricker(T-(d[i]/VP+tau)))/d[i] \
                  +np.outer(sv,1.5*aS*ricker(T-(d[i]/VS+tau)))/d[i]
        cols.append(out.ravel())
    return np.stack(cols,axis=1)                       # (M, 6)

def gen(src,M,nclut=30):
    G=green_cols(src,TAU); d=(G@M)*SNR
    D=d.reshape(NS,3,NT)
    for _ in range(nclut):
        i=rng.integers(NS); tc=rng.uniform(2,36); q=rng.normal(size=3); q/=np.linalg.norm(q)
        D[i]+=SNR*0.30*np.outer(q,ricker(T-tc))/28.
    D+=rng.normal(scale=SNR*0.03,size=D.shape)
    return D.ravel()

def tests(src,tau,dvec):
    G=green_cols(src,tau)
    m,res,_,_=np.linalg.lstsq(G,dvec,rcond=None)
    fit=G@m
    resid=np.linalg.norm(dvec-fit)/np.linalg.norm(dvec)
    M=sum(mi*B for mi,B in zip(m,BASIS))
    iso=abs(np.trace(M))/max(np.linalg.norm(M),1e-12)
    expl=np.linalg.norm(fit)/np.linalg.norm(dvec)
    return iso,resid,expl

# true source: pure double couple (deviatoric, zero trace)
SRC=np.array([6.,-8.,9.])
Mtrue=np.array([1.0,0.6,-0.3,0.4,-0.5,0.0])       # last component (isotropic) = 0
d=gen(SRC,Mtrue)

TAUS=np.arange(5.,7.01,0.05)   # 1/10 of a period: the LSQ fit is matched-filter-like and is sampling-sensitive
gx=np.arange(-21,21.1,3.5); gz=np.array([5.,9.,13.])
GP=[np.array([x,y,z]) for x in gx for y in gx for z in gz]
def best_over_tau(p,dvec):
    out=[]
    for tau in TAUS: out.append(tests(p,tau,dvec))
    return min(out,key=lambda t:t[1])              # best-fitting tau

print("isotropic fraction and fit residual, true source vs far false hypotheses\n")
iso_t,res_t,ex_t=best_over_tau(SRC,d)
print(f"  TRUE source     : iso={iso_t:.3f}  resid={res_t:.3f}  explained={ex_t:.3f}")
vals=[]
for p in GP:
    if np.linalg.norm(p-SRC)<6.: continue
    vals.append((best_over_tau(p,d),p))
isos=np.array([v[0][0] for v in vals]); ress=np.array([v[0][1] for v in vals]); exps=np.array([v[0][2] for v in vals])
print(f"  false hypotheses: iso median={np.median(isos):.3f}  10th pct={np.percentile(isos,10):.3f}")
print(f"                    resid median={np.median(ress):.3f}  min={ress.min():.3f}")
print(f"                    explained median={np.median(exps):.3f}  max={exps.max():.3f}")
print(f"\n  separation on isotropic fraction : true {iso_t:.3f} vs false 10th pct {np.percentile(isos,10):.3f}")
print(f"  separation on residual           : true {res_t:.3f} vs false min {ress.min():.3f}")
nbetter=int((ress<res_t).sum())
print(f"  false hypotheses with SMALLER residual than the truth: {nbetter}/{len(vals)}")
niso=int((isos<iso_t).sum())
print(f"  false hypotheses with SMALLER isotropic fraction     : {niso}/{len(vals)}")
