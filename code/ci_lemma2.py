"""
Tight form.  Instead of chaining worst-case inequalities, evaluate the leading-order
Fisher matrix EXACTLY.  It is still purely geometric: built from station directions,
distances and velocities, with no wave simulation.

   gamma  =  (||g'||^2 / sigma^2) * lam_min( sum_{(s,b) in A(T)} w_b^s v_b^s v_b^s^T )  +  O(1/omega)
   v_b^s  = (-n_s / V_b , 1)  in R^4

Corollaries tested:
  C1 (time)   A(T) empty  ->  leading term vanishes
  C2 (single branch) leading term = 0  iff all n_s lie on a common cone a.n = c V
  C3 (two branches)  leading term > 0  iff {n_s} spans R^3
"""
import numpy as np
from scipy.stats import spearmanr, linregress
rng=np.random.default_rng(11)
FC=2.0; VS=3.0; VP=1.73*VS; DT=1/(8*FC); SIGMA=0.02
SRC=np.array([0.,0.,8.,5.0]); AP,AS=1.0,1.5
def ricker(t):
    a=(np.pi*FC*t)**2; return (1-2*a)*np.exp(-a)
tt=np.arange(-4,4,DT); GD2=np.sum(np.gradient(ricker(tt),DT)**2)*DT

def atom(theta,STA,T,br='PS'):
    x,y,z,tau=theta; src=np.array([x,y,z]); v=STA-src
    d=np.linalg.norm(v,axis=1); n=v/d[:,None]; ez=np.array([0.,0.,1.])
    out=np.zeros((len(STA),3,len(T)))
    for i in range(len(STA)):
        sv=ez-np.dot(ez,n[i])*n[i]; nn=np.linalg.norm(sv)
        sv=sv/nn if nn>1e-9 else np.array([1.,0.,0.])
        if 'P' in br: out[i]+=AP*np.outer(n[i],ricker(T-(d[i]/VP+tau)))/d[i]
        if 'S' in br: out[i]+=AS*np.outer(sv,ricker(T-(d[i]/VS+tau)))/d[i]
    return out.ravel()
def gamma_num(STA,T,br='PS',h=np.array([.02,.02,.02,.005])):
    cols=[]
    for k in range(4):
        e=np.zeros(4); e[k]=h[k]
        cols.append((atom(SRC+e,STA,T,br)-atom(SRC-e,STA,T,br))/(2*h[k]))
    J=np.stack(cols,axis=1); return np.linalg.eigvalsh(J.T@J/SIGMA**2)[0]

def gamma_lead(STA,T,br='PS'):
    v=STA-SRC[:3]; d=np.linalg.norm(v,axis=1); n=v/d[:,None]
    M=np.zeros((4,4))
    for i in range(len(STA)):
        for b,V,A in (('P',VP,AP),('S',VS,AS)):
            if b not in br: continue
            ta=d[i]/V+SRC[3]
            if not (T[0]+0.5 < ta < T[-1]-0.5): continue      # capture test = time half of GCC
            vec=np.concatenate([-n[i]/V,[1.]]); w=(A/d[i])**2
            M+=w*np.outer(vec,vec)
    return GD2/SIGMA**2*np.linalg.eigvalsh(M)[0], M

def make(kind,N):
    if kind=='random':   a=rng.uniform(0,2*np.pi,N); r=40*np.sqrt(rng.uniform(.02,1,N))
    elif kind=='onesided':a=rng.uniform(0,np.pi*.7,N); r=40*np.sqrt(rng.uniform(.05,1,N))
    elif kind=='ring':   a=np.linspace(0,2*np.pi,N,endpoint=False)+rng.uniform(0,1); r=np.full(N,28.)
    elif kind=='clustered':a=rng.normal(1.,.25,N); r=rng.normal(30,4,N)
    elif kind=='linear': a=np.where(rng.random(N)<.5,.35,.35+np.pi); r=rng.uniform(5,40,N)
    else:                a=rng.uniform(0,2*np.pi,N); r=rng.uniform(3,12,N)
    return np.stack([r*np.cos(a),r*np.sin(a),np.zeros(N)],axis=1)

T=np.arange(0,32,DT); L=[];G=[];K=[]
for kind in ['random','onesided','ring','clustered','linear','near']:
    for N in [3,4,5,6,8,10,12,16]:
        for _ in range(2):
            STA=make(kind,N); l,_M=gamma_lead(STA,T); g=gamma_num(STA,T)
            L.append(l);G.append(g);K.append(kind)
L=np.array(L);G=np.array(G);K=np.array(K); ok=L>1e-12
print("TIGHT LEADING-ORDER PREDICTOR vs measured coercivity\n")
print(f"usable (leading term non-degenerate): {ok.sum()}/{len(L)}")
r,p=spearmanr(L[ok],G[ok]); print(f"Spearman      = {r:+.3f}  p={p:.1e}")
sl=linregress(np.log10(L[ok]),np.log10(G[ok]))
print(f"log-log slope = {sl.slope:.3f}   R^2={sl.rvalue**2:.3f}   (theory: slope 1)")
rat=G[ok]/L[ok]; print(f"ratio gamma/leading: median {np.median(rat):.2f}  IQR [{np.percentile(rat,25):.2f}, {np.percentile(rat,75):.2f}]")
print(f"\n{'family':>10} {'median ratio':>13} {'n':>4}")
for kk in ['random','onesided','ring','clustered','near','linear']:
    m=(K==kk)&ok
    print(f"{kk:>10} {np.median(G[m]/L[m]) if m.sum() else float('nan'):13.2f} {m.sum():4d}")
mdeg=(K=='linear')&(~ok)
print(f"\nlinear arrays with degenerate leading term: {mdeg.sum()}  -> measured gamma "
      f"{np.median(G[mdeg]) if mdeg.sum() else float('nan'):.3f} comes from polarization alone")

print("\nC2  single-branch cone degeneracy (leading term, P only)")
print(f"{'family':>10} {'lead P':>11} {'lead S':>11} {'lead P+S':>11}")
for kind in ['ring','linear','random','onesided','near']:
    STA=make(kind,8)
    print(f"{kind:>10} {gamma_lead(STA,T,'P')[0]:11.3e} {gamma_lead(STA,T,'S')[0]:11.3e} {gamma_lead(STA,T,'PS')[0]:11.3e}")
