"""
C-III, reformulated as MEDIUM-INDEPENDENT kinematic admissibility.

From the eikonal equation |grad T| = 1/V <= 1/Vmin, the traveltime field is Lipschitz with
constant 1/Vmin in ANY medium whose velocity is bounded below.  Hence, for two arrivals of
the SAME phase belonging to the SAME event,

        | t_i - t_j |  <=  | r_i - r_j | / Vmin                    (L)

and for the two branches at one station,   t_S > t_P   and  dt_PS = t_S - t_P >= 0,   with
dt_PS increasing with t_P when Vp/Vs is bounded                      (O)

No velocity model, no source location, no ray tracing.  These prune which COMBINATIONS of
picks can be one event, before any hypothesis is posited - a different level of action from
C-IV/C-V, which score a given hypothesis.

Measured: rejection rate of clutter-contaminated pick sets at zero rejection of true sets.
"""
import numpy as np
from itertools import combinations
rng=np.random.default_rng(31)
VS=3.0; VP=1.73*VS; VMIN=2.0            # only a LOWER BOUND on velocity is assumed
SRC=np.array([0.,0.,8.]); TAU=5.0
th=np.linspace(0,2*np.pi,8,endpoint=False)+.3
RR=np.array([6.,11.,18.,27.,36.,45.,54.,62.])
STA=np.stack([RR*np.cos(th),RR*np.sin(th),np.zeros(8)],axis=1)
d=np.linalg.norm(STA-SRC,axis=1)
tP=d/VP+TAU; tS=d/VS+TAU
DIJ=np.linalg.norm(STA[:,None,:]-STA[None,:,:],axis=2)

def lipschitz_ok(times, idx):
    for a,b in combinations(range(len(idx)),2):
        if abs(times[a]-times[b]) > DIJ[idx[a],idx[b]]/VMIN + 1e-9: return False
    return True

def order_ok(tp,ts):
    if np.any(ts<=tp): return False
    dps=ts-tp
    if len(tp)<3: return True
    return np.corrcoef(np.argsort(np.argsort(tp)),np.argsort(np.argsort(dps)))[0,1] > 0.5

print("test (L): same-phase Lipschitz;  test (O): P-S ordering and monotonicity\n")
print("true event:")
print(f"   (L) on true P set: {lipschitz_ok(tP,list(range(8)))}   (L) on true S set: {lipschitz_ok(tS,list(range(8)))}")
print(f"   (O) on true event: {order_ok(tP,tS)}\n")

# --- discrimination: random pick sets drawn from a clutter pool ---
for WIN in [30.,60.]:
    for k in [3,4,5,6,8]:
        rejL=0; rejO=0; rejAny=0; N=4000
        for _ in range(N):
            idx=sorted(rng.choice(8,k,replace=False))
            tp=rng.uniform(0,WIN,k)                       # random ("clutter") P times
            ts=tp+rng.uniform(0,WIN/3,k)                  # random S times after them
            L=lipschitz_ok(tp,idx) and lipschitz_ok(ts,idx)
            O=order_ok(tp,ts)
            rejL+= (not L); rejO+= (not O); rejAny+= (not (L and O))
        print(f"window {WIN:4.0f}s  k={k}:  rejected by (L) {100*rejL/N:5.1f}%   by (O) {100*rejO/N:5.1f}%   by either {100*rejAny/N:5.1f}%")
    print()

# --- realistic mixed case: true picks + one clutter substitution ---
print("substitution test: true event with j of its picks replaced by clutter")
for j in range(0,5):
    rej=0; N=3000
    for _ in range(N):
        tp=tP.copy(); ts=tS.copy()
        if j>0:
            sw=rng.choice(8,j,replace=False)
            tp[sw]=rng.uniform(0,30,j); ts[sw]=tp[sw]+rng.uniform(.5,8,j)
        L=lipschitz_ok(tp,list(range(8))) and lipschitz_ok(ts,list(range(8)))
        O=order_ok(tp,ts)
        rej+= (not (L and O))
    print(f"   {j} substituted picks -> rejected {100*rej/N:5.1f}%")
