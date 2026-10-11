"""
C-I verification, level 1: does a purely GEOMETRIC ray-coverage quantity predict the
coercivity constant of the reduced problem, and does the coercivity constant predict
actual estimator error?

Chain to test:   geometry  -->  gamma (coercivity)  -->  empirical epicenter error

Non-circularity: the geometric quantities are computed from station directions alone,
without ever forming the Jacobian or touching a waveform.  gamma is the smallest
eigenvalue of the whitened Fisher matrix of the FULL physical operator (P+S,
polarization, 1/d spreading).  The empirical error comes from Monte Carlo nonlinear
least squares on noisy data.  Three different objects.

theta = (x, y, z, tau).  Homogeneous half-space, straight rays, isotropic radiation.
"""
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import spearmanr
rng = np.random.default_rng(11)

FC=2.0; VS=3.0; VP=1.73*VS; DT=1/(8*FC)
SRC=np.array([0.,0.,8.,5.0])          # x,y,z,tau
SIGMA=0.02

def ricker(t):
    a=(np.pi*FC*t)**2; return (1-2*a)*np.exp(-a)

def atom(theta, STA, T, As=1.5):
    x,y,z,tau=theta; src=np.array([x,y,z])
    v=STA-src; d=np.linalg.norm(v,axis=1); n=v/d[:,None]
    ez=np.array([0.,0.,1.])
    out=np.zeros((len(STA),3,len(T)))
    for i in range(len(STA)):
        sv=ez-np.dot(ez,n[i])*n[i]; nn=np.linalg.norm(sv)
        sv=sv/nn if nn>1e-9 else np.array([1.,0.,0.])
        out[i]=np.outer(n[i],ricker(T-(d[i]/VP+tau)))/d[i] \
             + As*np.outer(sv,ricker(T-(d[i]/VS+tau)))/d[i]
    return out.ravel()

def jacobian(theta, STA, T, h=np.array([0.02,0.02,0.02,0.005])):
    cols=[]
    for k in range(4):
        e=np.zeros(4); e[k]=h[k]
        cols.append((atom(theta+e,STA,T)-atom(theta-e,STA,T))/(2*h[k]))
    return np.stack(cols,axis=1)

def coercivity(STA,T):
    J=jacobian(SRC,STA,T); F=J.T@J/SIGMA**2
    ev=np.linalg.eigvalsh(F)
    C=np.linalg.pinv(F,rcond=1e-15)
    crb_epi=np.sqrt(max(C[0,0]+C[1,1],0))
    return ev[0], crb_epi                      # gamma, CRB epicentral

# ---------- purely geometric quantities (no Jacobian, no waveform) ----------
def geometry_stats(STA):
    v=STA-SRC[:3]; d=np.linalg.norm(v,axis=1); n=v/d[:,None]
    Q=(n.T@n)/len(STA)                          # direction covariance, 3x3
    g_dir=np.linalg.eigvalsh(Q)[0]              # BLR-flavored: worst-covered direction
    az=np.sort(np.degrees(np.arctan2(STA[:,1]-SRC[1],STA[:,0]-SRC[0]))%360)
    gaps=np.diff(np.concatenate([az,[az[0]+360]]))
    return g_dir, gaps.max(), len(STA)          # coverage, azimuthal gap, N

# ---------- empirical estimator error ----------
def empirical_error(STA,T,ntrial=25):
    d0=atom(SRC,STA,T); errs=[]
    for _ in range(ntrial):
        d=d0+rng.normal(scale=SIGMA,size=d0.size)
        x0=SRC+np.array([rng.normal(0,.5),rng.normal(0,.5),rng.normal(0,.5),rng.normal(0,.05)])
        r=least_squares(lambda th: atom(th,STA,T)-d, x0, method='lm', xtol=1e-10, max_nfev=400)
        errs.append(np.hypot(r.x[0]-SRC[0], r.x[1]-SRC[1]))
    return float(np.median(errs))

# ---------- geometry generator: mixed families ----------
def make_geometry(kind,N):
    if kind=='random':
        a=rng.uniform(0,2*np.pi,N); r=40*np.sqrt(rng.uniform(0.02,1,N))
    elif kind=='onesided':
        a=rng.uniform(0,np.pi*0.7,N); r=40*np.sqrt(rng.uniform(0.05,1,N))
    elif kind=='ring':
        a=np.linspace(0,2*np.pi,N,endpoint=False)+rng.uniform(0,1); r=np.full(N,28.)
    elif kind=='clustered':
        a=rng.normal(1.0,0.25,N); r=rng.normal(30,4,N)
    elif kind=='linear':
        a=np.where(rng.random(N)<.5,0.35,0.35+np.pi); r=rng.uniform(5,40,N)
    elif kind=='near':
        a=rng.uniform(0,2*np.pi,N); r=rng.uniform(3,12,N)
    return np.stack([r*np.cos(a),r*np.sin(a),np.zeros(N)],axis=1)

T=np.arange(0,32,DT)
rows=[]
kinds=['random','onesided','ring','clustered','linear','near']
for kind in kinds:
    for N in [3,4,5,6,8,10,12,16]:
        for rep_ in range(2):
            STA=make_geometry(kind,N)
            g_dir,azgap,n_st=geometry_stats(STA)
            gam,crb=coercivity(STA,T)
            err=empirical_error(STA,T)
            rows.append((kind,N,g_dir,azgap,n_st,gam,crb,err))

import numpy as np
G=np.array([r[2] for r in rows]); AZ=np.array([r[3] for r in rows])
NS=np.array([r[4] for r in rows],float); GAM=np.array([r[5] for r in rows])
CRB=np.array([r[6] for r in rows]); ERR=np.array([r[7] for r in rows])
ok=np.isfinite(ERR)&np.isfinite(GAM)&(GAM>0)

print(f"{len(rows)} geometries, {ok.sum()} usable\n")
print("LEVEL 1  geometry -> coercivity      (does ray coverage predict gamma?)")
for name,q in [('directional coverage (BLR-like)',G),('azimuthal gap (classical)',-AZ),('station count',NS)]:
    rho,p=spearmanr(q[ok],GAM[ok]); print(f"   {name:34s} rho={rho:+.3f}  p={p:.2e}")
print("\nLEVEL 2  coercivity -> empirical error")
rho,p=spearmanr(GAM[ok],ERR[ok]); print(f"   gamma vs median epicentral error      rho={rho:+.3f}  p={p:.2e}")
rho,p=spearmanr(CRB[ok],ERR[ok]); print(f"   CRB   vs median epicentral error      rho={rho:+.3f}  p={p:.2e}")
print("\nLEVEL 3  geometry -> empirical error  (skipping gamma entirely)")
for name,q in [('directional coverage (BLR-like)',G),('azimuthal gap (classical)',-AZ),('station count',NS)]:
    rho,p=spearmanr(q[ok],ERR[ok]); print(f"   {name:34s} rho={rho:+.3f}  p={p:.2e}")

print("\nby geometry family (median values)")
print(f"{'family':>10} {'coverage':>9} {'az gap':>8} {'gamma':>10} {'CRB km':>8} {'err km':>8}")
for kind in kinds:
    m=[i for i,r in enumerate(rows) if r[0]==kind and ok[i]]
    print(f"{kind:>10} {np.median(G[m]):9.4f} {np.median(AZ[m]):8.1f} {np.median(GAM[m]):10.2e} {np.median(CRB[m]):8.3f} {np.median(ERR[m]):8.3f}")
np.save('/home/claude/ci_rows.npy', np.array([[r[2],r[3],r[4],r[5],r[6],r[7]] for r in rows]))
