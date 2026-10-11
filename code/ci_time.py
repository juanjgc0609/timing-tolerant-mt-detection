"""
C-I verification, level 4: the TIME half of the geometric control condition.
BLR requires rays to reach the observation region WITHIN time T.  Truncating the
recording window should collapse the coercivity constant exactly when the arrivals
of the outermost stations fall outside the window.  No purely spatial heuristic
(azimuthal gap, station count, GDOP) can see this.
"""
import numpy as np
FC=2.0; VS=3.0; VP=1.73*VS; DT=1/(8*FC); SIGMA=0.02
SRC=np.array([0.,0.,8.,5.0])
rng=np.random.default_rng(3)
a=np.linspace(0,2*np.pi,8,endpoint=False)+0.2
r=np.array([6,10,15,20,26,32,38,44.])
STA=np.stack([r*np.cos(a),r*np.sin(a),np.zeros(8)],axis=1)
d=np.linalg.norm(STA-SRC[:3],axis=1)
tP=d/VP+SRC[3]; tS=d/VS+SRC[3]
print("station distances (km):", np.round(d,1))
print("P arrivals (s):        ", np.round(tP,2))
print("S arrivals (s):        ", np.round(tS,2))
print(f"last P at {tP.max():.2f} s,  last S at {tS.max():.2f} s\n")

def ricker(t):
    x=(np.pi*FC*t)**2; return (1-2*x)*np.exp(-x)
def atom(theta,T,As=1.5):
    x,y,z,tau=theta; src=np.array([x,y,z])
    v=STA-src; dd=np.linalg.norm(v,axis=1); n=v/dd[:,None]; ez=np.array([0.,0.,1.])
    out=np.zeros((len(STA),3,len(T)))
    for i in range(len(STA)):
        sv=ez-np.dot(ez,n[i])*n[i]; sv/=np.linalg.norm(sv)
        out[i]=np.outer(n[i],ricker(T-(dd[i]/VP+tau)))/dd[i] \
             + As*np.outer(sv,ricker(T-(dd[i]/VS+tau)))/dd[i]
    return out.ravel()
def gamma_of(Tend,h=np.array([0.02,0.02,0.02,0.005])):
    T=np.arange(0,Tend,DT); cols=[]
    for k in range(4):
        e=np.zeros(4); e[k]=h[k]
        cols.append((atom(SRC+e,T)-atom(SRC-e,T))/(2*h[k]))
    J=np.stack(cols,axis=1); F=J.T@J/SIGMA**2
    ev=np.linalg.eigvalsh(F); C=np.linalg.pinv(F,rcond=1e-15)
    return ev[0], np.sqrt(max(C[0,0]+C[1,1],0)), np.sqrt(max(C[3,3],0))

print(f"{'T end':>7} {'stations w/ P':>13} {'w/ S':>6} | {'gamma':>11} {'CRB epi km':>11} {'CRB tau s':>10}")
print("-"*66)
prev=None
for Tend in [5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,15.0,17.0,20.0,25.0,32.0]:
    nP=int((tP<Tend-0.3).sum()); nS=int((tS<Tend-0.3).sum())
    g,ce,ct=gamma_of(Tend)
    jump="" if prev is None else f"  x{g/prev:6.2f}"
    print(f"{Tend:7.1f} {nP:13d} {nS:6d} | {g:11.3e} {ce:11.4f} {ct:10.5f}{jump}")
    prev=g
