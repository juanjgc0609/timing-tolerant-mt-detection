"""Appendix T at 30 replicates: coherent decoy, travel time vs polarization."""
import numpy as np
from scipy.stats import fisher_exact, beta
exec(open('/mnt/user-data/outputs/e1_critical.py').read().split('gx=np.arange')[0])
gx=np.arange(-18,18.1,3.0); gz=np.array([5.,8.,11.]); TAUS=np.arange(5.,8.51,0.25)
GP=[np.array([x,y,z]) for x in gx for y in gx for z in gz]
def outcome(dec,ndec,use,R=30):
    ok=0; spur=0
    for _ in range(R):
        D,pk=gen(dec,ndec); pre=pre_(D,pk)
        v=np.array([max(score(g,tt,pre,use) for tt in TAUS) for g in GP])
        b=GP[int(np.argmax(v))]
        if np.linalg.norm(b-SRC)<3.: ok+=1
        elif np.linalg.norm(b-dec)>3.: spur+=1
    return ok,spur,R
def ci(k,n,a=0.05):
    lo=beta.ppf(a/2,k,n-k+1) if k>0 else 0.
    hi=beta.ppf(1-a/2,k+1,n-k) if k<n else 1.
    return lo,hi
print("Apendice T con 30 replicas: tasa de solucion espuria\n")
print(f"{'sep':>5} {'sop':>6} | {'TT espurio':>12} {'IC95':>16} | {'C-IV espurio':>13} {'IC95':>16} | {'p':>8}")
print("-"*84)
for sep,nd in [(4.,6),(12.,12)]:
    dec=SRC+np.array([sep,0.,0.])
    o1,s1,R=outcome(dec,nd,{'CIII'}); o2,s2,_=outcome(dec,nd,{'CIII','CIV'})
    l1,h1=ci(s1,R); l2,h2=ci(s2,R)
    _,p=fisher_exact([[s1,R-s1],[s2,R-s2]])
    print(f"{sep:5.0f} {nd:4d}/16 | {s1:4d}/{R} = {s1/R:4.2f} [{l1:.2f},{h1:.2f}] | "
          f"{s2:4d}/{R} = {s2/R:4.2f} [{l2:.2f},{h2:.2f}] | {p:8.4f}")
