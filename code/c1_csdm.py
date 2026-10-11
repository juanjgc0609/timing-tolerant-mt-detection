"""
C1: does the INTER-STATION cross-spectral matrix count sources, where the per-station 3x3
covariance does not?

For one source the array response across stations at frequency w is a rank-one outer product;
K mutually incoherent sources give rank K; noise that is independent between distant sites
adds a roughly isotropic floor.  So the number of eigenvalues above that floor is a source
counter that uses no picks, no arrival times, no origin time and no velocity model.

Model order is selected by MDL, the standard criterion in array processing, and compared
against the per-station 3x3 polarization rank that Appendix M showed cannot count sources.
"""
import numpy as np
rng=np.random.default_rng(4242)
FS=20.; DT=1/FS; T=np.arange(0,80,DT); NT=len(T)
VS=3.0; VP=1.73*VS; FC=2.0
a=np.linspace(0,2*np.pi,10,endpoint=False)+.3
r=np.array([9.,13.,18.,23.,29.,35.,41.,47.,52.,58.])
STA=np.stack([r*np.cos(a),r*np.sin(a),np.zeros(10)],axis=1); NS=len(STA)

def stf(n):
    """independent band-limited source-time function per event"""
    x=rng.normal(size=n)
    X=np.fft.rfft(x); f=np.fft.rfftfreq(n,DT)
    X*= np.exp(-((f-FC)/1.2)**2)
    return np.fft.irfft(X,n)

def pol(src,rr,b):
    v=rr-src; n=v/np.linalg.norm(v)
    if b=='P': return n
    ez=np.array([0.,0.,1.]); sv=ez-np.dot(ez,n)*n; nn=np.linalg.norm(sv)
    return sv/nn if nn>1e-9 else np.array([1.,0.,0.])

def gen(srcs,taus,snr=6.0,nclut=40,clut_coherent=False):
    """vertical component only for the N x N test; clutter either incoherent between stations
       or (stress test) radiated by a decoy so that it IS coherent between stations"""
    D=np.zeros((NS,NT))
    for src,tau in zip(srcs,taus):
        s=stf(NT); d=np.linalg.norm(STA-src,axis=1)
        for i in range(NS):
            for b,A,V in (('P',1.0,VP),('S',1.5,VS)):
                sh=int((d[i]/V+tau)/DT)
                w=snr*A*pol(src,STA[i],b)[2]/d[i]
                D[i]+=w*np.roll(s,sh)
    if clut_coherent:
        dec=np.array([18.,-15.,9.]); s=stf(NT); dd=np.linalg.norm(STA-dec,axis=1)
        for i in range(NS):
            sh=int((dd[i]/VP+12.)/DT)
            D[i]+=snr*0.5*pol(dec,STA[i],'P')[2]/dd[i]*np.roll(s,sh)
    for _ in range(nclut):
        i=rng.integers(NS); s=stf(NT)*0.0
        j=rng.integers(NT); amp=rng.uniform(.3,1.0)
        pulse=np.zeros(NT); w=int(1.0/DT)
        pulse[j:j+w]=amp*snr*0.25*rng.normal(size=min(w,NT-j))
        D[i]+=pulse
    D+=rng.normal(scale=snr*0.05,size=D.shape)
    return D

def csdm_eigs(D,fband=(1.2,3.0),nseg=24):
    """cross-spectral density matrix between stations, averaged over snapshots and band"""
    L=NT//nseg; R=np.zeros((NS,NS),dtype=complex); cnt=0
    win=np.hanning(L)
    for k in range(nseg):
        seg=D[:,k*L:(k+1)*L]*win
        X=np.fft.rfft(seg,axis=1); f=np.fft.rfftfreq(L,DT)
        sel=(f>=fband[0])&(f<=fband[1])
        for m in np.where(sel)[0]:
            x=X[:,m:m+1]; R+=x@x.conj().T; cnt+=1
    R/=cnt
    return np.sort(np.real(np.linalg.eigvalsh(R)))[::-1], nseg*int(((fband[1]-fband[0])*NT/nseg*DT))

def mdl(eigs,Nsnap):
    """minimum description length model order selection"""
    p=len(eigs); best=(np.inf,0)
    for k in range(p):
        lam=eigs[k:]
        if np.any(lam<=0): break
        g=np.exp(np.mean(np.log(lam))); ar=np.mean(lam)
        L=-Nsnap*(p-k)*np.log(g/ar)
        pen=0.5*k*(2*p-k)*np.log(Nsnap)
        if L+pen<best[0]: best=(L+pen,k)
    return best[1]

def rank3c(D):
    """the per-station 3x3 approach, for contrast: here only one component, so report n/a"""
    return None

print("source counting from the INTER-STATION cross-spectral matrix\n")
print(f"{'true K':>7} {'clutter':>10} | {'MDL estimate':>13} | {'top eigenvalues (normalized)':>34}")
print("-"*78)
CFG=[([np.array([5.,6.,9.])],[10.]),
     ([np.array([5.,6.,9.]),np.array([-12.,8.,7.])],[10.,22.]),
     ([np.array([5.,6.,9.]),np.array([-12.,8.,7.]),np.array([3.,-14.,11.])],[10.,22.,34.])]
for srcs,taus in CFG:
    for coh in (False,True):
        D=gen(srcs,taus,clut_coherent=coh)
        ev,nsnap=csdm_eigs(D)
        k=mdl(ev,max(nsnap,NS+1))
        top=" ".join(f"{v/ev[0]:5.2f}" for v in ev[:5])
        print(f"{len(srcs):7d} {'coherent' if coh else 'random':>10} | {k:13d} | {top:>34}")
