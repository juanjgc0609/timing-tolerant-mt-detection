"""
C0 decisive test: is source separability governed by COHERENCE of the physical
records (polarization + amplitude + moveout) rather than by arrival-time separation?

Matched comparison, same geometry, same pulse, same normalization:
  mu_wf : coherence of full 3C physical atoms (polarization + 1/d spreading + P,S moveout)
  mu_tt : coherence of arrival-time-only atoms (same t_P,t_S; no polarization, no amplitude)
Lower coherence = more separable. If mu_wf << mu_tt the wavefield adds real
separating power over travel times, and by how much, and in which directions.
"""
import numpy as np
FC=2.0; VS=3.0; VP=1.73*VS; DT=1/(8*FC); T=np.arange(0,30,DT)

th=np.linspace(0,2*np.pi,12,endpoint=False)
STA=np.stack([25*np.cos(th),25*np.sin(th),np.zeros(12)],axis=1)   # ring, 25 km
SRC0=np.array([0.,0.,8.])

def ricker(t):
    a=(np.pi*FC*t)**2
    return (1-2*a)*np.exp(-a)

def atom(src, physical=True, As=1.5):
    out=[]
    for r in STA:
        v=r-src; d=np.linalg.norm(v); n=v/d
        ez=np.array([0.,0.,1.]); sv=ez-np.dot(ez,n)*n
        sv=sv/np.linalg.norm(sv) if np.linalg.norm(sv)>1e-9 else np.array([1.,0.,0.])
        gP=ricker(T-d/VP); gS=ricker(T-d/VS)
        if physical:
            tr=np.outer(n,gP)/d + As*np.outer(sv,gS)/d          # 3 x Nt
        else:
            e=np.array([0.,0.,1.])
            tr=np.outer(e,gP) + np.outer(e,gS)                   # no polariz., no spreading
        out.append(tr.ravel())
    a=np.concatenate(out); return a/np.linalg.norm(a)

def mu(s1,s2,physical): return abs(np.dot(atom(s1,physical),atom(s2,physical)))

dirs={'radial  (along x)':np.array([1.,0,0]),
      'transverse (y)   ':np.array([0,1.,0]),
      'vertical (depth) ':np.array([0,0,1.])}
seps=[0.25,0.5,1.0,2.0,4.0,8.0]

print(f"coherence mu  (lower = more separable).  fc={FC} Hz, Vs={VS}, Vp/Vs=1.73, 12 stations\n")
print(f"{'offset direction':18} {'sep km':>7} | {'mu_waveform':>11} {'mu_traveltime':>13} | {'ratio':>6}")
print("-"*66)
res={}
for name,u in dirs.items():
    for s in seps:
        s2=SRC0+s*u
        if s2[2]<1.0: continue
        mw=mu(SRC0,s2,True); mt=mu(SRC0,s2,False)
        res[(name,s)]=(mw,mt)
        print(f"{name:18} {s:7.2f} | {mw:11.4f} {mt:13.4f} | {mt/max(mw,1e-9):6.2f}x")
    print()

# separability threshold: smallest separation with mu < 0.5
print("separation at which mu drops below 0.5  (proxy for the demixing limit)")
print(f"{'direction':18} | {'waveform':>9} {'travel time':>12}")
print("-"*46)
fine=np.linspace(0.05,15,600)
for name,u in dirs.items():
    r={}
    for phys,key in [(True,'wf'),(False,'tt')]:
        val=np.nan
        for s in fine:
            s2=SRC0+s*u
            if s2[2]<1.0: continue
            if mu(SRC0,s2,phys)<0.5: val=s; break
        r[key]=val
    print(f"{name:18} | {r['wf']:9.2f} {r['tt']:12.2f}")
