"""
B1: does the Proposition survive interfaces and multipathing?

Two-layer model: layer 1 velocity V1 down to depth H, halfspace V2 > V1.  Two arrivals per
phase: the DIRECT ray and the critically refracted HEAD WAVE.  Both are analytic.

  direct:  t = D/V1 ,             source slowness  p = n/V1     (n = straight-line direction)
  head  :  t = x/V2 + (2H - z) cos(thc)/V1 ,
           source slowness p has magnitude 1/V1 at the CRITICAL angle thc, azimuth toward
           the receiver, downgoing.   sin(thc) = V1/V2.

Eikonal check: |grad_source t| must equal 1/V1 for BOTH branches.
Certificate under test:  gamma_lead = lam_min( sum over ALL branches of w_b v_b v_b^T ),
                         v_b = (-p_b, 1).
"""
import numpy as np
from scipy.stats import spearmanr
rng=np.random.default_rng(404)
FC=2.0; DT=1/(8*FC); T=np.arange(0,40,DT); NT=len(T); SIGMA=0.02
V1S,V2S=2.6,4.2; VPR=1.73; H=14.0
SRC=np.array([0.,0.,7.,6.0])
def vel(b,layer): return (V1S if layer==1 else V2S)*(VPR if b=='P' else 1.)
def ricker(t):
    a=(np.pi*FC*t)**2; return (1-2*a)*np.exp(-a)

def branches(src,r,b):
    """returns list of (traveltime, source-slowness vector, amplitude weight)"""
    V1=vel(b,1); V2=vel(b,2); out=[]
    v=r-src; D=np.linalg.norm(v); n=v/D
    out.append((D/V1, n/V1, 1.0/D))                       # direct
    x=np.linalg.norm((r-src)[:2]); z=src[2]
    if V2>V1 and x>0:
        s=V1/V2; c=np.sqrt(max(1-s*s,0.))
        th=x/V2+(2*H-z)*c/V1
        h=(r-src)[:2]/max(x,1e-9)
        p=np.array([s*h[0]/V1, s*h[1]/V1, -c/V1])         # critical angle, downgoing
        xc=2*H*np.sqrt((V2-V1)/(V2+V1))                   # crossover
        if x>0.8*xc: out.append((th, p, 0.55/max(x,1.)))  # weaker amplitude
    return out

def atom(theta,STA):
    src=theta[:3]; tau=theta[3]; out=np.zeros((len(STA),3,NT))
    for i,r in enumerate(STA):
        for b,A in (('P',1.0),('S',1.5)):
            for tt,p,w in branches(src,r,b):
                n=p/np.linalg.norm(p)
                if b=='P': pol=n
                else:
                    ez=np.array([0.,0.,1.]); sv=ez-np.dot(ez,n)*n; nn=np.linalg.norm(sv)
                    pol=sv/nn if nn>1e-9 else np.array([1.,0.,0.])
                out[i]+=A*w*np.outer(pol,ricker(T-(tt+tau)))
    return out.ravel()
def gamma_num(STA,h=np.array([.02,.02,.02,.005])):
    cols=[]
    for k in range(4):
        e=np.zeros(4); e[k]=h[k]
        cols.append((atom(SRC+e,STA)-atom(SRC-e,STA))/(2*h[k]))
    J=np.stack(cols,axis=1); return np.linalg.eigvalsh(J.T@J/SIGMA**2)[0]
def gamma_lead(STA,use_head=True):
    M=np.zeros((4,4))
    for r in STA:
        for b,A in (('P',1.0),('S',1.5)):
            bl=branches(SRC[:3],r,b)
            if not use_head: bl=bl[:1]
            for tt,p,w in bl:
                v=np.concatenate([-p,[1.]]); M+=(A*w)**2*np.outer(v,v)
    return np.linalg.eigvalsh(M)[0]

print("eikonal check, |grad_src t| vs 1/V1 for both branches")
r0=np.array([34.,0.,0.])
for b in ('P','S'):
    V1=vel(b,1)
    for name,(tt,p,w) in zip(('direct','head'),branches(SRC[:3],r0,b)):
        print(f"   {b} {name:>6}: |p|={np.linalg.norm(p):.6f}  1/V1={1/V1:.6f}  "
              f"rel.err={abs(np.linalg.norm(p)-1/V1)*V1:.2e}   t={tt:.3f}s")
print()

def net(N,kind):
    if kind=='wide':  r=np.linspace(8,55,N); a=np.linspace(0,2*np.pi,N,endpoint=False)+.3
    elif kind=='far': r=np.linspace(30,60,N); a=np.linspace(0,2*np.pi,N,endpoint=False)+.7
    elif kind=='near':r=np.linspace(5,18,N);  a=np.linspace(0,2*np.pi,N,endpoint=False)+1.1
    else: a=rng.uniform(0,2*np.pi,N); r=55*np.sqrt(rng.uniform(.03,1,N))
    return np.stack([r*np.cos(a),r*np.sin(a),np.zeros(N)],axis=1)

L=[];Ld=[];G=[]
for kind in ['wide','far','near','random']:
    for N in [4,5,6,8,10]:
        for _ in range(2):
            S=net(N,kind)
            L.append(gamma_lead(S,True)); Ld.append(gamma_lead(S,False)); G.append(gamma_num(S))
L=np.array(L);Ld=np.array(Ld);G=np.array(G); ok=(L>1e-14)&(G>0)
print(f"{ok.sum()} geometries with a two-layer medium (direct + head wave)\n")
r1=spearmanr(L[ok],G[ok]).statistic; r2=spearmanr(Ld[ok],G[ok]).statistic
print(f"Spearman(leading with branch sum , gamma)   = {r1:+.4f}     <- B1 requires >= 0.95")
print(f"Spearman(leading, DIRECT branch only, gamma) = {r2:+.4f}     <- ignoring multipathing")
rat=G[ok]/L[ok]
print(f"ratio gamma/leading: median {np.median(rat):.3f}  IQR [{np.percentile(rat,25):.3f}, {np.percentile(rat,75):.3f}]")
res=rat/np.median(rat)
print(f"95% of geometries within {np.percentile(np.abs(res-1)*100,95):.1f}% of the closed-form prediction")
print(f"\nhead wave contribution: median gamma_lead ratio (with head / direct only) = {np.median(L[ok]/np.maximum(Ld[ok],1e-30)):.2f}x")
