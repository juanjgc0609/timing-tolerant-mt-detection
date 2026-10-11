"""
Heterogeneous extension of the Proposition.

By Fermat / eikonal, d T / d x_source = - p_source, the slowness vector at the SOURCE end
of the connecting ray, with |p| = 1/V(x_source), independent of the path.  Hence

    v_b^s = ( -p_b^s , 1 )        replaces      ( -n_s / V_b , 1 )

Two testable consequences:
  (i) the velocity entering is the LOCAL velocity at the source, not a path average;
  (ii) the direction is the ray TAKE-OFF direction, not the straight line to the station.

Medium: constant velocity gradient v(z) = v0 + k z, for which the two-point traveltime is
known in closed form (Slotnick):   T = (1/k) arccosh( 1 + k^2 D^2 / (2 v1 v2) ),
D = straight-line distance.  Slowness vectors obtained by differentiating that closed form,
validated against the eikonal identity |grad T| = 1/v.
"""
import numpy as np
from scipy.stats import spearmanr
rng=np.random.default_rng(5)
FC=2.0; DT=1/(8*FC); SIGMA=0.02
VS0, VPR = 2.2, 1.73          # S velocity at surface, Vp/Vs
SRC=np.array([0.,0.,8.,5.0]); AP,AS=1.0,1.5

def vel(z,k,branch): return (VS0+k*z)*(VPR if branch=='P' else 1.0)
def ttime(x1,x2,k,branch):
    v1=vel(x1[2],k,branch); v2=vel(x2[2],k,branch); D=np.linalg.norm(x2-x1)
    if k<1e-9: return D/v1
    return np.arccosh(1+ k**2*D**2/(2*v1*v2))/k
def grad_src(x1,x2,k,branch,h=1e-4):
    g=np.zeros(3)
    for j in range(3):
        e=np.zeros(3); e[j]=h
        g[j]=(ttime(x1+e,x2,k,branch)-ttime(x1-e,x2,k,branch))/(2*h)
    return g                                   # = -p_source
def grad_rec(x1,x2,k,branch,h=1e-4):
    g=np.zeros(3)
    for j in range(3):
        e=np.zeros(3); e[j]=h
        g[j]=(ttime(x1,x2+e,k,branch)-ttime(x1,x2-e,k,branch))/(2*h)
    return g

def ricker(t):
    a=(np.pi*FC*t)**2; return (1-2*a)*np.exp(-a)
def atom(theta,STA,T,k):
    src=theta[:3]; tau=theta[3]
    out=np.zeros((len(STA),3,len(T)))
    for i,r in enumerate(STA):
        D=np.linalg.norm(r-src)
        for b,A in (('P',AP),('S',AS)):
            tb=ttime(src,r,k,b)+tau
            gr=grad_rec(src,r,k,b); nr=gr/np.linalg.norm(gr)      # arrival direction at receiver
            if b=='P': pol=nr
            else:
                ez=np.array([0.,0.,1.]); sv=ez-np.dot(ez,nr)*nr
                pol=sv/np.linalg.norm(sv) if np.linalg.norm(sv)>1e-9 else np.array([1.,0.,0.])
            out[i]+=A*np.outer(pol,ricker(T-tb))/D
    return out.ravel()
def gamma_num(STA,T,k,h=np.array([.02,.02,.02,.005])):
    cols=[]
    for j in range(4):
        e=np.zeros(4); e[j]=h[j]
        cols.append((atom(SRC+e,STA,T,k)-atom(SRC-e,STA,T,k))/(2*h[j]))
    J=np.stack(cols,axis=1); return np.linalg.eigvalsh(J.T@J/SIGMA**2)[0]

def lead(STA,T,k,mode='takeoff'):
    M=np.zeros((4,4))
    for r in STA:
        D=np.linalg.norm(r-SRC[:3])
        for b,A in (('P',AP),('S',AS)):
            if not (T[0]+.5 < ttime(SRC[:3],r,k,b)+SRC[3] < T[-1]-.5): continue
            if mode=='takeoff':
                p=-grad_src(SRC[:3],r,k,b)                       # true source slowness
            else:                                                # naive straight-line + surface velocity
                n=(r-SRC[:3])/D; p=n/vel(SRC[2],k,b)
            v=np.concatenate([-p,[1.]]); M+=(A/D)**2*np.outer(v,v)
    return np.linalg.eigvalsh(M)[0]

# ---- eikonal validation ----
r0=np.array([25.,0.,0.])
for k in [0.0,0.05,0.12]:
    p=-grad_src(SRC[:3],r0,k,'S')
    print(f"k={k:.2f}  |p_src|={np.linalg.norm(p):.6f}  1/v(src)={1/vel(SRC[2],k,'S'):.6f}  "
          f"rel.err={abs(np.linalg.norm(p)-1/vel(SRC[2],k,'S'))/(1/vel(SRC[2],k,'S')):.2e}")
print()

def make(kind,N):
    if kind=='random':  a=rng.uniform(0,2*np.pi,N); rr=45*np.sqrt(rng.uniform(.05,1,N))
    elif kind=='ring':  a=np.linspace(0,2*np.pi,N,endpoint=False)+rng.uniform(0,1); rr=np.full(N,30.)
    elif kind=='onesided':a=rng.uniform(0,np.pi*.7,N); rr=45*np.sqrt(rng.uniform(.05,1,N))
    else:               a=rng.uniform(0,2*np.pi,N); rr=rng.uniform(4,14,N)
    return np.stack([rr*np.cos(a),rr*np.sin(a),np.zeros(N)],axis=1)

T=np.arange(0,40,DT)
print(f"{'grad k':>7} {'n':>4} | {'rho(takeoff, gamma)':>20} {'rho(straight, gamma)':>21} | {'med |dev| takeoff':>18} {'straight':>10}")
print("-"*94)
for k in [0.0,0.03,0.06,0.10,0.15]:
    Lt=[];Ls=[];G=[]
    for kind in ['random','ring','onesided','near']:
        for N in [4,6,8,12]:
            for _ in range(2):
                STA=make(kind,N)
                lt=lead(STA,T,k,'takeoff'); ls=lead(STA,T,k,'straight'); g=gamma_num(STA,T,k)
                if lt>1e-10 and ls>1e-10: Lt.append(lt);Ls.append(ls);G.append(g)
    Lt=np.array(Lt);Ls=np.array(Ls);G=np.array(G)
    rt=spearmanr(Lt,G).statistic; rs=spearmanr(Ls,G).statistic
    dt=np.median(np.abs(G/Lt/np.median(G/Lt)-1))*100
    ds=np.median(np.abs(G/Ls/np.median(G/Ls)-1))*100
    print(f"{k:7.2f} {len(G):4d} | {rt:20.4f} {rs:21.4f} | {dt:17.1f}% {ds:9.1f}%")
