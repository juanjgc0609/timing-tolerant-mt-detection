import sys, numpy as np
sys.path.insert(0,'/home/claude/work'); sys.path.insert(0,'/home/claude/work/gamma_src')
import task1x_mtglrt as X
from task0_common import DEFAULT_TAUS, DT, NT, VP, VS, make_network
from task1o_consistent import draw_sources, matched
AX=np.arange(-21.0,21.0+1e-9,3.0)+0.37; DEP=np.array([5.,9.,13.]); TAUS=np.asarray(DEFAULT_TAUS)
SH=X.SH; K=X.CAPV
def T_node(R,STA,src,Q):
    v,tpk,_,_=X.node_vectors(R,STA,src,TAUS); vf=v.reshape(len(TAUS),-1)
    return np.sum((vf@Q)**2,axis=1), v, tpk      # T per tau (no 1/2)
def bound(R,STA,c,h,tc,dt):
    E=np.sum(R**2,axis=1)       # (ns,NT) 3C energy per sample
    d=np.linalg.norm(STA-c,axis=1); B=0.0
    for V in (VP,VS):
        L=h/V+dt; t1=d/V+tc-L; t2=d/V+tc+L
        k1=np.clip(np.floor(t1/DT).astype(int),0,NT-1); k2=np.clip(np.floor(t2/DT).astype(int),0,NT-1)
        for i in range(len(STA)):
            lo=max(0,k1[i]-SH); hi=min(NT-1,k2[i]+SH)
            B+=min(K**2, E[i,lo:hi+1].max())
    return B
rng=np.random.default_rng(3); rm=np.random.default_rng(4); rh=np.random.default_rng(5)
viol=0; ncell=0; worst=np.inf
for w in range(4):
    STA=X.make_network("onesided",8,seed=100+w) if w%2==0 else make_network("random",4,seed=100+w)
    srcs,taus=draw_sources(rng)
    D,_,_,_=X.generate_x(STA,srcs,taus,20,rng,rm,True,rh)
    MF,sd=matched(D); R0=MF/sd
    grid=[np.array([x,y,z]) for x in AX for y in AX for z in DEP]
    Gs=[X.radiation_operator(g,STA) for g in grid]; Qs=[X.basis(G) for G in Gs]
    # one deflation step at the global max to get a residual
    best=(-1,None,None,None)
    for gi,g in enumerate(grid):
        Tt,v,tpk=T_node(R0,STA,g,Qs[gi]); j=int(Tt.argmax())
        if Tt[j]>best[0]: best=(Tt[j],gi,v[j].reshape(-1),tpk[j])
    _,gi,vf,tpk=best; m,*_=np.linalg.lstsq(Gs[gi],vf,rcond=None); pred=(Gs[gi]@m).reshape(2,len(STA),3)
    R1=R0.copy()
    for b in range(2):
        for i in range(len(STA)):
            t0=int(tpk[b,i]); lo,hi=max(0,t0-X.AC_HALF),min(NT,t0+X.AC_HALF+1)
            R1[i,:,lo:hi]-=np.outer(pred[b,i],X.AC[X.AC_HALF-(t0-lo):X.AC_HALF+(hi-t0)])
    for R in (R0,R1):
        for _ in range(60):
            sx=[rng.integers(1,8) for _ in range(2)]+[rng.integers(1,4)]; st=rng.integers(1,32)
            ox=[rng.integers(0,len(AX)-sx[0]+1),rng.integers(0,len(AX)-sx[1]+1),rng.integers(0,3-sx[2]+1)]; ot=rng.integers(0,len(TAUS)-st+1)
            xs,ys,zs,ts=AX[ox[0]:ox[0]+sx[0]],AX[ox[1]:ox[1]+sx[1]],DEP[ox[2]:ox[2]+sx[2]],TAUS[ot:ot+st]
            c=np.array([(xs[0]+xs[-1])/2,(ys[0]+ys[-1])/2,(zs[0]+zs[-1])/2])
            h=np.sqrt(((xs[-1]-xs[0])/2)**2+((ys[-1]-ys[0])/2)**2+((zs[-1]-zs[0])/2)**2)
            tc,dt=(ts[0]+ts[-1])/2,(ts[-1]-ts[0])/2
            B=bound(R,STA,c,h,tc,dt); Tmax=-1
            for x in xs:
                for y in ys:
                    for z in zs:
                        g=np.array([x,y,z]); Q=X.basis(X.radiation_operator(g,STA))
                        v,tpk,_,_=X.node_vectors(R,STA,g,ts); Tt=np.sum((v.reshape(len(ts),-1)@Q)**2,axis=1)
                        Tmax=max(Tmax,Tt.max())
            ncell+=1; worst=min(worst,B-Tmax); viol+=int(Tmax>B+1e-9)
print('cells',ncell,'violations',viol,'min slack B-maxT',worst)
