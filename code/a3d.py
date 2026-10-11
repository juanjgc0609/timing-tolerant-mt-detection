import numpy as np, statistics
from scipy.optimize import minimize
exec(open('/home/claude/a3b.py').read().split("sets=[('TT',set())")[0])

def score1(src,tau,AMP,PD,picks,use):
    d=np.linalg.norm(STA-src,axis=1)
    tt=np.stack([d/VPH,d/VSH])+tau; idx=np.clip((tt/DT).astype(int),0,NT-1)
    hit=np.array([[any(abs(tt[b,i]-pk)<TOL for pk in picks[i]) for i in range(8)] for b in range(2)])
    if hit.sum()==0: return 0.
    if 'CIII' in use:
        cl=[[None]*8,[None]*8]
        for b in range(2):
            for i in range(8):
                if hit[b,i]: cl[b][i]=min((pk for pk in picks[i] if abs(tt[b,i]-pk)<TOL),key=lambda p:abs(tt[b,i]-p))
        for b in range(2):
            ii=[i for i in range(8) if cl[b][i] is not None]
            for a_ in range(len(ii)):
                for c_ in range(a_+1,len(ii)):
                    if abs(cl[b][ii[a_]]-cl[b][ii[c_]])>DIJ[ii[a_],ii[c_]]/VMIN+1e-9: return 0.
        for i in range(8):
            if cl[0][i] is not None and cl[1][i] is not None and cl[1][i]<=cl[0][i]: return 0.
    s=hit.sum()/16.
    if 'CIV' in use:
        cs=[abs(np.dot(pol_from_dir((STA[i]-src)/d[i],b),PD[i,idx[bi,i]]))
            for bi,b in enumerate(('P','S')) for i in range(8) if hit[bi,i]]
        s*= np.mean(cs)**2 if cs else 0.
    if 'CV' in use:
        m=[AMP[i,idx[bi,i]]*d[i]/A for bi,(b,A) in enumerate((('P',AP),('S',AS))) for i in range(8) if hit[bi,i]]
        if len(m)>=2:
            m=np.array(m); s*=1./(1.+m.std()/max(m.mean(),1e-9))
        else: s*=0.5
    return s

gx=np.arange(-18,18.1,6.0); gz=np.array([4.,8.,12.]); TAUS=np.arange(4.5,8.01,0.5)
G0=[(np.array([x,y,z]),tt) for x in gx for y in gx for z in gz for tt in TAUS]
def best_continuous(D,picks,use):
    AMP,PD=feats(D)
    vals=[score1(s,t,AMP,PD,picks,use) for s,t in G0]
    s0,t0=G0[int(np.argmax(vals))]
    x0=np.array([s0[0],s0[1],s0[2],t0])
    f=lambda p: -score1(np.array([p[0],p[1],max(1.,min(20.,p[2]))]),p[3],AMP,PD,picks,use)
    r=minimize(f,x0,method='Nelder-Mead',options=dict(maxiter=400,xatol=.05,fatol=1e-4))
    return float(np.hypot(r.x[0]-SRC[0],r.x[1]-SRC[1]))

sets2=[('TT',set()),('C-IV',{'CIV'}),('MIN=III+IV',{'CIII','CIV'}),('ALL',{'CIII','CIV','CV'})]
REPS=12
print(f"CONTINUOUS optimization, median epicentral error over {REPS} realizations [km]\n")
print(f"{'k':>5} {'model err':>10} | " + " ".join(f"{n:>12}" for n,_ in sets2))
print("-"*70)
errs=[]; res={n:[] for n,_ in sets2}
for k in [0.0,0.06,0.12,0.20]:
    e=[]
    for i in range(8):
        for b,Vh in (('P',VPH),('S',VSH)):
            tt_=ttrue(SRC,STA[i],k,b); e.append((np.linalg.norm(STA[i]-SRC)/Vh-tt_)/tt_)
    err=100*np.sqrt(np.mean(np.array(e)**2)); errs.append(err)
    med={}
    for n,u in sets2:
        v=[]
        for _ in range(REPS):
            D,picks=gen(k); v.append(best_continuous(D,picks,u))
        med[n]=statistics.median(v); res[n].append(med[n])
    print(f"{k:5.2f} {err:9.1f}% | " + " ".join(f"{med[n]:12.2f}" for n,_ in sets2))
print("\ndegradation slope (km per 10% traveltime model error)")
for n,_ in sets2:
    print(f"   {n:>12}: {np.polyfit(np.array(errs),np.array(res[n]),1)[0]*10:+6.2f}")
