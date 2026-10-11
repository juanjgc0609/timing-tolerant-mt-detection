"""C1, corrected test: origin time is a FREE parameter, and network density varies.
theta = (a1,a2, x1,x2, tau1,tau2).  Equal-energy unit-norm atoms => equal SNR.
The P-S differential decouples distance from origin time; that decoupling can
only matter where the network does not already over-determine it.
"""
import numpy as np
FC=2.0; VS=3.0; Z0=8.0; SNR=20.0; sigma=1/SNR
DT=1.0/(8*FC); T=np.arange(0.0,25.0,DT)

def make_atom(XS):
    def atom(x,tau,vp,ap,as_):
        d=np.sqrt((x-XS)**2+Z0**2); A=np.zeros((len(XS),len(T)))
        if ap: A+=(ap/d)[:,None]*np.sinc(2*FC*(T[None,:]-(d/vp+tau)[:,None]))
        if as_:A+=(as_/d)[:,None]*np.sinc(2*FC*(T[None,:]-(d/VS+tau)[:,None]))
        v=A.ravel(); return v/np.linalg.norm(v)
    return atom

def sd_delta(D,vp,ap,as_,atom,h=1e-3):
    x1,x2=-D/2,D/2; t1=t2=5.0
    cols=[atom(x1,t1,vp,ap,as_),atom(x2,t2,vp,ap,as_)]
    for (x,t) in [(x1,t1),(x2,t2)]:
        cols.append((atom(x+h,t,vp,ap,as_)-atom(x-h,t,vp,ap,as_))/(2*h))
    for (x,t) in [(x1,t1),(x2,t2)]:
        cols.append((atom(x,t+h,vp,ap,as_)-atom(x,t-h,vp,ap,as_))/(2*h))
    B=np.stack(cols,axis=1)
    J=B.T@B/sigma**2
    s=np.linalg.svd(J,compute_uv=False)
    if s[-1]/s[0]<1e-13: return np.inf           # structurally degenerate
    C=np.linalg.inv(J)
    return np.sqrt(max(C[2,2]+C[3,3]-2*C[2,3],0.0))

def delta_res(vp,ap,as_,atom,lo=0.005,hi=40.0,tol=1e-3):
    f=lambda D: sd_delta(D,vp,ap,as_,atom)-D
    if not np.isfinite(f(hi)) or f(hi)>0: return np.nan
    while hi-lo>tol:
        mid=.5*(lo+hi)
        if f(mid)>0: lo=mid
        else: hi=mid
    return hi

print("Origin time FREE.  Resolution limit Delta_res [km]; nan = not identifiable\n")
print(f"{'Ns':>3} | {'P only':>9} {'S only':>9} {'P+S':>9} | {'gain of P+S vs best single':>27}")
print("-"*70)
for ns in [1,2,3,4,6,12]:
    XS=np.linspace(-30,30,ns) if ns>1 else np.array([12.0])
    atom=make_atom(XS); vp=1.73*VS
    dP=delta_res(vp,1.,0.,atom); dS=delta_res(vp,0.,1.,atom); dPS=delta_res(vp,1.,1.,atom)
    singles=[d for d in (dP,dS) if np.isfinite(d) and not np.isnan(d)]
    if singles and not np.isnan(dPS):
        g=f"{100*(min(singles)-dPS)/min(singles):+.1f}%"
    else:
        g="single branches not identifiable" if not singles else "n/a"
    print(f"{ns:3d} | {dP:9.4f} {dS:9.4f} {dPS:9.4f} | {g:>27}")

print("\nVp/Vs sweep at the sparse-network regime (Ns=2), origin time free")
print(f"{'Vp/Vs':>6} | {'dRes P+S':>9} {'kappa':>8}")
print("-"*30)
XS=np.linspace(-30,30,2); atom=make_atom(XS); out=[]
for r in [1.50,1.60,1.67,1.73,1.80,1.90,1.95]:
    d=delta_res(r*VS,1.,1.,atom); out.append((r,d))
    print(f"{r:6.2f} | {d:9.4f} {d*FC/VS:8.4f}")
out=np.array(out)
print(f"\nslope d(Delta_res)/d(Vp/Vs) = {np.polyfit(out[:,0],out[:,1],1)[0]:+.4f} km per unit ratio")
print(f"relative change 1.50 -> 1.95 = {100*(out[-1,1]-out[0,1])/out[0,1]:+.2f}%")
