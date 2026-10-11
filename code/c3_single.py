"""
C3: single three-component station, multiple overlapping sources.

Structural prediction: with ONE branch the score field is a ridge along the back-azimuth,
because polarization fixes direction but t_P cannot fix distance while tau is unknown.  With
TWO branches, t_S - t_P fixes distance independently of tau and the ridge collapses to a point.

Multiple sources are extracted by hierarchical subtraction, the strategy used for overlapping
gravitational-wave signals: take the strongest peak, subtract its best-fit contribution, repeat.
"""
import numpy as np, statistics
rng=np.random.default_rng(303)
FC=2.0; VS=3.0; VP=1.73*VS; DT=1/(8*FC); T=np.arange(0,45,DT); NT=len(T)
SNR=8.0; AP,AS=1.0,1.5
STA=np.array([[0.,0.,0.]])                      # ONE three-component station
DEPTH=8.0
def ricker(t):
    a=(np.pi*FC*t)**2; return (1-2*a)*np.exp(-a)
WAV=ricker(np.arange(-3,3,DT))
def pol(src,r,b):
    v=r-src; n=v/np.linalg.norm(v)
    if b=='P': return n
    ez=np.array([0.,0.,1.]); sv=ez-np.dot(ez,n)*n; nn=np.linalg.norm(sv)
    return sv/nn if nn>1e-9 else np.array([1.,0.,0.])
def wave(src,tau,branches,amp=1.0):
    d=np.linalg.norm(STA[0]-src); out=np.zeros((3,NT))
    for b,A,V in (('P',AP,VP),('S',AS,VS)):
        if b not in branches: continue
        out+=amp*A*np.outer(pol(src,STA[0],b),ricker(T-(d/V+tau)))/d
    return out
def gen(srcs,taus,branches,nclut):
    D=np.zeros((3,NT))
    for s,tt in zip(srcs,taus): D+=SNR*wave(s,tt,branches)
    for _ in range(nclut):
        tc=rng.uniform(2,42); q=rng.normal(size=3); q/=np.linalg.norm(q)
        D+=SNR*0.30*np.outer(q,ricker(T-tc))/25.
    D+=rng.normal(scale=SNR*0.02,size=D.shape)
    return D
def matched(D,branches,GX,GY,TAUS):
    R=np.stack([np.convolve(D[c],WAV,mode='same') for c in range(3)],axis=1)  # (NT,3)
    best=None
    for x in GX:
        for y in GY:
            src=np.array([x,y,DEPTH]); d=np.linalg.norm(STA[0]-src)
            ws=[];ps=[]
            for b,A,V in (('P',AP,VP),('S',AS,VS)):
                if b not in branches: continue
                ws.append((A/d,V)); ps.append(pol(src,STA[0],b))
            for tau in TAUS:
                num=0.;nrm=0.
                for (w,V),p in zip(ws,ps):
                    j=int((d/V+tau)/DT)
                    if not (0<=j<NT): continue
                    num+=w*np.dot(p,R[j]); nrm+=w*w
                s=max(num,0.)/np.sqrt(nrm+1e-15)
                if best is None or s>best[0]: best=(s,x,y,tau)
    return best
def extract(D,branches,nmax,GX,GY,TAUS):
    Dw=D.copy(); found=[]
    for _ in range(nmax):
        b=matched(Dw,branches,GX,GY,TAUS)
        if b is None: break
        s,x,y,tau=b; src=np.array([x,y,DEPTH])
        w=wave(src,tau,branches)
        a=np.sum(Dw*w)/max(np.sum(w*w),1e-12)
        if a<=0: break
        Dw=Dw-a*w; found.append((x,y,tau,s,a))
    return found

GX=np.arange(-30,30.1,2.5); GY=np.arange(-30,30.1,2.5); TAUS=np.arange(3.,10.01,0.05)
print("C3  single 3C station.  epicentral error of the recovered sources [km]\n")
print(f"{'sources':>8} {'branches':>9} | " + " ".join(f"{'src'+str(i+1):>8}" for i in range(3)) + "   spurious")
print("-"*62)
CONF=[([np.array([12.,9.,DEPTH])],[6.0]),
      ([np.array([12.,9.,DEPTH]),np.array([-15.,6.,DEPTH])],[6.0,7.2]),
      ([np.array([12.,9.,DEPTH]),np.array([-15.,6.,DEPTH]),np.array([4.,-18.,DEPTH])],[6.0,7.2,8.1])]
for srcs,taus in CONF:
    for br in ('P','PS'):
        D=gen(srcs,taus,'PS',nclut=12)                 # data always contain both branches
        f=extract(D,br,len(srcs)+1,GX,GY,TAUS)
        errs=[]
        for s in srcs:
            if f: errs.append(min(np.hypot(x-s[0],y-s[1]) for x,y,_,_,_ in f))
            else: errs.append(np.nan)
        matched_n=sum(1 for x,y,_,_,_ in f if min(np.hypot(x-s[0],y-s[1]) for s in srcs)<5.)
        print(f"{len(srcs):8d} {br:>9} | " + " ".join(f"{e:8.2f}" for e in errs) +
              " "*(8*(3-len(errs))+3*(3-len(errs))) + f"   {len(f)-matched_n}")
    print()
