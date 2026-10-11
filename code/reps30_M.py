"""Appendix M at 30 replicates: channel separation in the spread geometry at high clutter."""
import numpy as np
from scipy.stats import fisher_exact, beta
exec(open('/mnt/user-data/outputs/a2_minimal.py').read().split("SETS=[('TT',set())")[0])
def ci(k,n,a=0.05):
    lo=beta.ppf(a/2,k,n-k+1) if k>0 else 0.
    hi=beta.ppf(1-a/2,k+1,n-k) if k<n else 1.
    return lo,hi
SETS=[('TT',set()),('C-IV',{'CIV'}),('MIN',{'CIII','CIV'}),('ALL',{'CII','CIII','CIV','CV'})]
R=30
print("Apendice M con 30 replicas: tasa de evento falso\n")
for kind in ['spread','random']:
    STA=geom(kind); DIJ=np.linalg.norm(STA[:,None,:]-STA[None,:,:],axis=2)
    for nc in [40]:
        fe={n:0 for n,_ in SETS}
        for _ in range(R):
            D,pk=gen(STA,nc); pre=precompute(D,pk)
            for n,u in SETS:
                f,_m=evaluate(STA,DIJ,pre,u); fe[n]+=f
        print(f"  geometria {kind}, clutter {nc}")
        for n,_ in SETS:
            l,h=ci(fe[n],R)
            print(f"    {n:>5}: {fe[n]:2d}/{R} = {fe[n]/R:.2f}  IC95 [{l:.2f},{h:.2f}]")
        _,p=fisher_exact([[fe['TT'],R-fe['TT']],[fe['C-IV'],R-fe['C-IV']]])
        _,p2=fisher_exact([[fe['C-IV'],R-fe['C-IV']],[fe['ALL'],R-fe['ALL']]])
        print(f"    TT vs C-IV: p={p:.4f}   C-IV vs ALL: p={p2:.4f}\n")
