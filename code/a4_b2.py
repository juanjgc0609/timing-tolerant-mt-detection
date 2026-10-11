"""
A4: cost-effectiveness of the minimal sufficient set.
B2: does a network designed by maximizing the observability constant beat one designed by
    minimizing azimuthal gap, at equal station count?
"""
import numpy as np, time, statistics
exec(open('/home/claude/a2_minimal.py').read().split("SETS=[('TT',set())")[0])

# ---------------- A4 : measured cost per objective evaluation ----------------
STA=geom('spread'); DIJ=np.linalg.norm(STA[:,None,:]-STA[None,:,:],axis=2)
D,pk=gen(STA,30); pre=precompute(D,pk)
CH=[('TT',set()),('C-III',{'CIII'}),('C-IV',{'CIV'}),('C-V',{'CV'}),('C-II',{'CII'}),
    ('MIN=III+IV',{'CIII','CIV'}),('ALL',{'CII','CIII','CIV','CV'})]
print("A4  measured cost per candidate evaluation (mean of 4000 calls)\n")
print(f"{'channel set':>12} {'us/eval':>9} {'x TT':>7}")
print("-"*32)
base=None; costs={}
for n,u in CH:
    t0=time.perf_counter()
    for _ in range(4000): score(SRC+np.array([1.,1.,0.]),6.1,STA,DIJ,pre,u)
    dt=(time.perf_counter()-t0)/4000*1e6
    costs[n]=dt
    if base is None: base=dt
    print(f"{n:>12} {dt:9.1f} {dt/base:7.2f}")
print(f"""
reference costs from Section 12 of the proposal, same unit:
   analytic beam synthesis + transport misfit, per objective evaluation   38000.0 us
   direct elastodynamic solve, per source evaluation                 576000000.0 us

A4 ratio, minimal set vs the beam+transport waveform condition:
   cost  {costs['MIN=III+IV']:.1f} us / 38000 us = {100*costs['MIN=III+IV']/38000:.4f}% of the waveform condition
   cost  {costs['C-IV']:.1f} us / 38000 us = {100*costs['C-IV']/38000:.4f}%  (polarization alone)
   benefit (Appendix M): matched the full constraint stack in 8/8 cells
   -> A4 registered margin was >=90% of benefit at <1% of cost: met by three to four orders
      of magnitude on cost, with the benefit condition already established in Appendix M.
   Against the direct-solver condition the cost ratio is {100*costs['MIN=III+IV']/5.76e8:.8f}%.
""")

# ---------------- B2 : network design ----------------
def gamma_lead(STA_):
    M=np.zeros((4,4))
    for r in STA_:
        d=np.linalg.norm(r-SRC); n=(r-SRC)/d
        for V,A in ((VP,AP),(VS,AS)):
            v=np.concatenate([-n/V,[1.]]); M+=(A/d)**2*np.outer(v,v)
    return np.linalg.eigvalsh(M)[0]
def azgap(STA_):
    az=np.sort(np.degrees(np.arctan2(STA_[:,1]-SRC[1],STA_[:,0]-SRC[0]))%360)
    return np.diff(np.concatenate([az,[az[0]+360]])).max()

def random_net(N,rmax=45.):
    a=rng.uniform(0,2*np.pi,N); r=rmax*np.sqrt(rng.uniform(.02,1,N))
    return np.stack([r*np.cos(a),r*np.sin(a),np.zeros(N)],axis=1)

print("B2  network design at equal station count\n")
print(f"{'N':>3} | {'criterion':>22} {'gamma':>10} {'az gap':>8} {'median epic. err km':>20} {'false-event':>12}")
print("-"*82)
for N in [5,7]:
    pool=[random_net(N) for _ in range(4000)]
    g=np.array([gamma_lead(s) for s in pool]); a=np.array([azgap(s) for s in pool])
    picks={'max observability':pool[int(np.argmax(g))],
           'min azimuthal gap':pool[int(np.argmin(a))],
           'random (median gamma)':pool[int(np.argsort(g)[len(g)//2])]}
    for name,S in picks.items():
        DJ=np.linalg.norm(S[:,None,:]-S[None,:,:],axis=2)
        errs=[]; fe=0; REPS=10
        for _ in range(REPS):
            Dx,pkx=gen(S,25); px=precompute(Dx,pkx)
            v=np.array([max(score(gp,tt,S,DJ,px,{'CIII','CIV'}) for tt in TAUS) for gp in GP])
            j=int(np.argmax(v)); errs.append(float(np.hypot(GP[j][0]-SRC[0],GP[j][1]-SRC[1])))
            st=max(score(SRC,tt,S,DJ,px,{'CIII','CIV'}) for tt in TAUS)
            fe+= int(v[FAR].max()>=st)
        print(f"{N:3d} | {name:>22} {gamma_lead(S):10.3f} {azgap(S):8.1f} {statistics.median(errs):20.2f} {fe/REPS:12.2f}")
    print()
