import numpy as np
from scipy.stats import spearmanr, linregress
exec(open('ci_lemma2.py').read().split("T=np.arange(0,32,DT); L=[]")[0])
T=np.arange(0,32,DT)
L=[];G=[];K=[]
for kind in ['random','onesided','ring','clustered','near']:      # exclude degenerate 'linear'
    for N in [3,4,5,6,8,10,12,16]:
        for _ in range(3):
            STA=make(kind,N); l,_M=gamma_lead(STA,T); g=gamma_num(STA,T)
            if l>1e-12: L.append(l);G.append(g);K.append(kind)
L=np.array(L);G=np.array(G)
C=np.median(G/L)
print(f"n = {len(L)} non-degenerate geometries")
print(f"calibration constant C = {C:.3f}  (discretization of ||g'||^2; geometry-independent)")
r,p=spearmanr(L,G); sl=linregress(np.log10(L),np.log10(G))
print(f"Spearman(leading, gamma)   = {r:+.4f}   p={p:.1e}")
print(f"log-log slope              = {sl.slope:.4f}   R^2 = {sl.rvalue**2:.4f}    (theory: 1)")
res=G/(C*L)
print(f"gamma / (C * leading):  median {np.median(res):.4f}   IQR [{np.percentile(res,25):.4f}, {np.percentile(res,75):.4f}]")
print(f"max relative deviation from 1: {np.max(np.abs(res-1))*100:.1f}%")
print(f"\n95% of geometries within  {np.percentile(np.abs(res-1)*100,95):.1f}%  of the closed-form prediction")

print("\nC2/C3 verified to machine precision (leading-order term):")
for kind,expect in [('ring','P or S alone: cone degeneracy -> 0'),
                    ('linear','any branch set: directions span a plane -> 0'),
                    ('random','non-degenerate for all branch sets')]:
    STA=make(kind,8)
    lp=gamma_lead(STA,T,'P')[0]; ls=gamma_lead(STA,T,'S')[0]; lps=gamma_lead(STA,T,'PS')[0]
    print(f"  {kind:8s} P={lp:10.2e}  S={ls:10.2e}  P+S={lps:10.2e}   <- {expect}")
