"""Own first-arrival travel times for the registered JMA2001 points (linear between points), by ray shooting in fine
layers (dz = 0.02 km), up-going and diving (turning) rays; receivers at the surface. VS = VP/1.73."""
import numpy as np
PTS=[(0.0,4.800),(1.5,5.220),(2.5,5.450),(4.0,5.720),(5.0,5.820),(15.0,6.220),(22.0,6.544),(34.0,7.188),(41.0,7.494),(44.0,7.598),(48.0,7.688),(52.5,7.750),(65.0,7.845)]
DZ=0.02
ZB=np.arange(0.0,60.0+1e-9,DZ)
def vel(z,r=1.0): return np.interp(z,[p[0] for p in PTS],[p[1] for p in PTS])/r
def curve(zs,r=1.0):
    """return list of (X,T) branches sampled densely in p"""
    zt=np.arange(0.0,zs,DZ); h=np.minimum(zt+DZ,zs)-zt; v=vel(zt+h/2,r)
    vmax=v.max()
    XT=[]
    # up-going
    ps=np.linspace(0,1/vmax*(1-1e-9),4000)
    for p in ps:
        q=np.sqrt(1-(p*v)**2); XT.append((np.sum(h*p*v/q),np.sum(h/(v*q))))
    # down-going diving: turning below source
    zd=np.arange(zs,60.0,DZ); vd=vel(zd+DZ/2,r)
    pmax=1/vmax
    for j in range(1,len(zd)):
        vt=vd[j]
        if vt<=max(vmax,vd[:j].max()): continue
        p=1/vt*(1-1e-7)
        vv=np.r_[v,vd[:j]]; hh=np.r_[h,np.full(j,DZ)]
        q=np.sqrt(np.maximum(1-(p*vv)**2,1e-14))
        x=np.sum(h*p*v/q[:len(v)])+2*np.sum(DZ*p*vd[:j]/q[len(v):]); t=np.sum(h/(v*q[:len(v)]))+2*np.sum(DZ/(vd[:j]*q[len(v):]))
        XT.append((x,t))
    XT=np.array(XT); return XT
_cache={}
def tt(zs,X,r=1.0):
    k=(round(zs,3),r)
    if k not in _cache:
        XT=curve(zs,r); xs=np.arange(0,80.01,0.05); T=np.full(len(xs),np.inf)
        # first arrival: for each x, min over segments of the sampled curve (piecewise linear interp)
        o=XT
        for a,b in zip(o[:-1],o[1:]):
            lo,hi=min(a[0],b[0]),max(a[0],b[0])
            m=(xs>=lo)&(xs<=hi)
            if not m.any(): continue
            if hi-lo<1e-12: t=np.full(m.sum(),min(a[1],b[1]))
            else: t=a[1]+(xs[m]-a[0])*(b[1]-a[1])/(b[0]-a[0])
            T[m]=np.minimum(T[m],t)
        _cache[k]=(xs,T)
    xs,T=_cache[k]; return np.interp(X,xs,T)
if __name__=="__main__":
    import sys; sys.path.insert(0,"/home/claude/work/explore_contrib")
    import c18_tt as TT
    print("c18_tt JMA2001 points equal registered:", TT.JMA2001_POINTS==PTS)
    X=np.array([5.,10.,12.,15.,20.,25.,28.,30.,35.])
    for zs in (6.3,7.0,7.8,8.4,3.0,11.0):
        for r in (1.0,1.73):
            mine=tt(zs,X,r); theirs=TT.first_arrival(TT.JMA2001,zs,X,r)
            print("z %.1f %s max|diff| %.4f s  (mine %s)"%(zs,"P" if r==1 else "S",np.max(np.abs(mine-theirs)),np.round(mine,3)))
