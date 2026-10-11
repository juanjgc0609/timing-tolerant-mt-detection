"""Appendix R at 30 replicates: pick-free advantage at low picker miss rate."""
import numpy as np
from scipy.stats import fisher_exact, beta
exec(open('/mnt/user-data/outputs/d1_coupling.py').read().split('SETS=[')[0])
def ci(k,n,a=0.05):
    lo=beta.ppf(a/2,k,n-k+1) if k>0 else 0.
    hi=beta.ppf(1-a/2,k+1,n-k) if k<n else 1.
    return lo,hi
S=net(6); DIJ=np.linalg.norm(S[:,None,:]-S[None,:,:],axis=2); R=30
print("Apendice R con 30 replicas: tasa de recuperacion\n")
print(f"{'perdida':>8} | {'TT':>16} | {'sin picks':>16} | {'dif':>6} {'p':>8}")
print("-"*60)
for m in [0.2,0.3,0.5]:
    a=b=0
    for _ in range(R):
        D,pk=gen(S,30,m); pre=precompute(D,pk)
        a+=recovered(S,DIJ,pre,set()); b+=recovered(S,DIJ,pre,{'WF'})
    la,ha=ci(a,R); lb,hb=ci(b,R)
    _,p=fisher_exact([[a,R-a],[b,R-b]])
    print(f"{m:8.1f} | {a:2d}/{R}={a/R:.2f} [{la:.2f},{ha:.2f}] | {b:2d}/{R}={b/R:.2f} [{lb:.2f},{hb:.2f}] | "
          f"{b/R-a/R:+6.2f} {p:8.4f}")
