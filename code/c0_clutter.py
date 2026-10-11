"""
The real failure mode: a FALSE source hypothesis supported by clutter.
Margin = (score at true source) - (best score among hypotheses >5 km away),
measured in units of the std of the score field.  Higher = harder to fool.
  travel-time score : fraction of predicted (tP,tS) matched by a pick within +-0.5 s
  waveform score    : normalized correlation with the physical 3C atom
"""
import numpy as np
rng=np.random.default_rng(7)
FC=2.0; VS=3.0; VP=1.73*VS; DT=1/(8*FC); T=np.arange(0,30,DT); NT=len(T)
SRC0=np.array([0.,0.,8.])
th=np.linspace(0,2*np.pi,6,endpoint=False)+.3
STA=np.stack([25*np.cos(th),25*np.sin(th),np.zeros(6)],axis=1)
def ricker(t):
    a=(np.pi*FC*t)**2; return (1-2*a)*np.exp(-a)
def arrivals(src):
    d=np.linalg.norm(STA-src,axis=1); return d/VP, d/VS, d
def atom(src,As=1.5):
    tp,ts,d=arrivals(src); out=[]
    for i,r in enumerate(STA):
        v=r-src; n=v/d[i]
        ez=np.array([0.,0.,1.]); sv=ez-np.dot(ez,n)*n; sv/=np.linalg.norm(sv)
        out.append((np.outer(n,ricker(T-tp[i]))/d[i]+As*np.outer(sv,ricker(T-ts[i]))/d[i]).ravel())
    a=np.concatenate(out); return a/np.linalg.norm(a)

def make_data(nclutter,snr=8.0):
    d=snr*atom(SRC0); picks=[[] for _ in STA]
    tp,ts,_=arrivals(SRC0)
    for i in range(len(STA)): picks[i]+= [tp[i],ts[i]]
    for _ in range(nclutter):
        i=rng.integers(len(STA)); tc=rng.uniform(2,25)
        v=rng.normal(size=3); v/=np.linalg.norm(v)       # random polarization
        amp=rng.uniform(.5,2.0)
        blk=np.zeros((len(STA),3,NT)); blk[i]=amp*np.outer(v,ricker(T-tc))
        d=d+blk.ravel()/np.linalg.norm(blk.ravel())*snr*0.35
        picks[i].append(tc)
    d=d+rng.normal(scale=np.linalg.norm(d)/np.sqrt(d.size)/3,size=d.size)
    return d,picks

def tt_score(src,picks,tol=0.5):
    tp,ts,_=arrivals(src); hit=0
    for i in range(len(STA)):
        for tt in (tp[i],ts[i]):
            if any(abs(tt-p)<tol for p in picks[i]): hit+=1
    return hit/(2*len(STA))
def wf_score(src,d): return abs(np.dot(atom(src),d))/np.linalg.norm(d)

gx=np.arange(-24,24.1,3.0); gz=np.array([4.,8.,12.,16.])
grid=[np.array([x,y,z]) for x in gx for y in gx for z in gz]
print(f"{'clutter':>8} | {'TRAVEL TIME':^30} | {'WAVEFORM':^30}")
print(f"{'events':>8} | {'true':>7} {'best false':>10} {'margin':>9} | {'true':>7} {'best false':>10} {'margin':>9}")
print("-"*78)
for nc in [0,5,15,30,60]:
    d,picks=make_data(nc)
    st=np.array([tt_score(g,picks) for g in grid]); sw=np.array([wf_score(g,d) for g in grid])
    far=np.array([np.linalg.norm(g-SRC0)>5.0 for g in grid])
    tt_true=tt_score(SRC0,picks); wf_true=wf_score(SRC0,d)
    mt=(tt_true-st[far].max())/max(st.std(),1e-9)
    mw=(wf_true-sw[far].max())/max(sw.std(),1e-9)
    print(f"{nc:8d} | {tt_true:7.3f} {st[far].max():10.3f} {mt:8.2f}s | {wf_true:7.3f} {sw[far].max():10.3f} {mw:8.2f}s")
