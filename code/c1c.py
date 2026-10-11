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

def per_bin_counts(D,fband=(1.0,3.2),nseg=40,overlap=0.5):
    """narrowband treatment: one CSDM PER FREQUENCY BIN, never averaged across frequency.
       At a fixed omega each source contributes one steering vector, so K sources give rank K."""
    L=int(NT/(nseg*(1-overlap)+overlap)); step=int(L*(1-overlap))
    win=np.hanning(L); segs=[]
    k=0
    while k+L<=NT:
        segs.append(D[:,k:k+L]*win); k+=step
    segs=np.array(segs)                                  # (Nsnap, NS, L)
    X=np.fft.rfft(segs,axis=2); f=np.fft.rfftfreq(L,DT)
    sel=np.where((f>=fband[0])&(f<=fband[1]))[0]
    Nsnap=len(segs); out=[]
    for m in sel:
        x=X[:,:,m]                                       # (Nsnap, NS)
        R=(x.conj().T@x)/Nsnap
        ev=np.sort(np.real(np.linalg.eigvalsh(R)))[::-1]
        out.append((mdl(ev,Nsnap),ev))
    return out,Nsnap

def mdl(eigs,Nsnap):
    p=len(eigs); best=(np.inf,0)
    for k in range(p):
        lam=eigs[k:]
        if np.any(lam<=1e-30): break
        g=np.exp(np.mean(np.log(lam))); ar=np.mean(lam)
        L=-Nsnap*(p-k)*np.log(g/ar); pen=0.5*k*(2*p-k)*np.log(Nsnap)
        if L+pen<best[0]: best=(L+pen,k)
    return best[1]

import collections
print("source counting vs TEMPORAL OVERLAP of the events\n")
print(f"{'true K':>7} {'spread s':>9} {'clutter':>10} | {'median K':>9} {'mode K':>7} | {'gap at K':>9}")
print("-"*62)
S3=[np.array([5.,6.,9.]),np.array([-12.,8.,7.]),np.array([3.,-14.,11.])]
for K in (2,3):
    for spread in (0.0,1.0,3.0,12.0):
        taus=[10.+j*spread for j in range(K)]
        for coh in (False,True):
            D=gen(S3[:K],taus,clut_coherent=coh)
            res,Nsnap=per_bin_counts(D)
            ks=[k for k,_ in res]
            gaps=[ev[K-1]/max(ev[K],1e-30) for _,ev in res]
            print(f"{K:7d} {spread:9.1f} {'coherent' if coh else 'random':>10} | "
                  f"{int(np.median(ks)):9d} {collections.Counter(ks).most_common(1)[0][0]:7d} | {np.median(gaps):9.2f}")
    print()
