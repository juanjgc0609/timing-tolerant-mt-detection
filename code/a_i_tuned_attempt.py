import numpy as np, statistics
exec(open('/home/claude/a_i_tuned.py').read().split('print("A-i TUNED')[0])
def intersect(o1,o2):
    """intersect LINES, not rays: the eigenvector sign is arbitrary."""
    p1=STA[o1['sta']][:2]; d1=o1['az']
    p2=STA[o2['sta']][:2]; d2=o2['az']
    A=np.array([[d1[0],-d2[0]],[d1[1],-d2[1]]])
    if abs(np.linalg.det(A))<1e-6: return None
    ts=np.linalg.solve(A,p2-p1)
    x=p1+ts[0]*d1
    return x if np.hypot(*x)<70 else None
SRC=np.array([7.,-5.,8.])                       # OFF the origin, so failure cannot look like success
def locate2(obs,ntrial=900):
    best=(-1,None); n=len(obs)
    if n<3: return None
    for _ in range(ntrial):
        a_,b_=rng.choice(n,2,replace=False)
        if obs[a_]['sta']==obs[b_]['sta']: continue
        xy=intersect(obs[a_],obs[b_])
        if xy is None: continue
        sc,inl=support(xy,obs)
        if len(inl)<3 or not lip_ok(obs,inl): continue
        for _ in range(12):
            J=[];r=[]
            for q in inl:
                o=obs[q]; v=xy-STA[o['sta']][:2]; nn=np.linalg.norm(v)
                if nn<1e-6: continue
                u=v/nn; perp=np.array([-o['az'][1],o['az'][0]]); w=o['rect']/sigma_az(o['rect'])
                J.append(w*perp); r.append(-w*np.dot(perp,u)*nn)
            if len(J)<3: break
            J=np.array(J); r=np.array(r)
            try: dx=np.linalg.lstsq(J,r,rcond=None)[0]
            except Exception: break
            if not np.all(np.isfinite(dx)): break
            xy=xy+np.clip(dx,-3,3)
        sc,inl=support(xy,obs)
        if sc>best[0]: best=(sc,xy.copy())
    return best[1]
def err_of(xy):
    return np.nan if xy is None else float(np.hypot(xy[0]-SRC[0],xy[1]-SRC[1]))
TAUS=np.arange(5.,8.01,0.15)
GP2=[np.array([x,y,z]) for x in np.arange(-20,20.1,2.5) for y in np.arange(-20,20.1,2.5) for z in [5.,8.,11.]]
print("corrected: true source off the origin; failures reported as failures\n")
print(f"{'k':>5} {'model err':>10} | {'TT':>8} {'C-IV':>8} | {'A-i tuned':>10} {'fail rate':>10}")
print("-"*58)
for k in [0.0,0.06,0.12,0.20]:
    e=[]
    for i in range(8):
        for b,Vh in (('P',VPH),('S',VSH)):
            tt_=ttrue(SRC,STA[i],k,b); e.append((np.linalg.norm(STA[i]-SRC)/Vh-tt_)/tt_)
    merr=100*np.sqrt(np.mean(np.array(e)**2))
    ett=[];eciv=[];eaz=[];nf=0
    for _ in range(8):
        D,pk=gen(k); pre=pre_(D,pk)
        v=np.array([max(score_pred(g,tt,pre,set()) for tt in TAUS) for g in GP2])
        g=GP2[int(np.argmax(v))]; ett.append(np.hypot(g[0]-SRC[0],g[1]-SRC[1]))
        v=np.array([max(score_pred(g,tt,pre,{'CIV'}) for tt in TAUS) for g in GP2])
        g=GP2[int(np.argmax(v))]; eciv.append(np.hypot(g[0]-SRC[0],g[1]-SRC[1]))
        xy=locate2(obs_full(pk,D))
        if xy is None: nf+=1
        else: eaz.append(err_of(xy))
    m=statistics.median(eaz) if eaz else float('nan')
    print(f"{k:5.2f} {merr:9.1f}% | {statistics.median(ett):8.2f} {statistics.median(eciv):8.2f} | {m:10.2f} {nf/8:10.2f}")
print("\nlateral back-azimuth deviation, at k=0.12")
def obs_lat(pk,D,bias,jit):
    o=obs_full(pk,D); bb={i:np.radians(rng.normal(0,bias)) for i in range(8)}
    for e in o:
        th=bb[e['sta']]+np.radians(rng.normal(0,jit)); c,s=np.cos(th),np.sin(th)
        e['az']=np.array([c*e['az'][0]-s*e['az'][1], s*e['az'][0]+c*e['az'][1]])
    return o
print(f"{'bias':>6} {'jitter':>7} | {'median err km':>14} {'fail rate':>10}")
for bias,jit in [(0.,0.),(3.,2.),(6.,4.),(10.,6.),(15.,8.)]:
    ee=[];nf=0
    for _ in range(10):
        D,pk=gen(0.12); xy=locate2(obs_lat(pk,D,bias,jit))
        if xy is None: nf+=1
        else: ee.append(err_of(xy))
    print(f"{bias:5.0f}d {jit:6.0f}d | {statistics.median(ee) if ee else float('nan'):14.2f} {nf/10:10.2f}")
