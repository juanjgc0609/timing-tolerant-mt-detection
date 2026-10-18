import sys, math, json
import numpy as np
from obspy import read, UTCDateTime
D="/home/claude/realdata/MatchLocate2/example"
T0=UTCDateTime("2012-09-02T03:20:00")
jma=[]
for l in open(D+"/catalog.dat"):
    f=l.split()
    if not f: continue
    t=UTCDateTime(f[0].replace("/","-")+"T"+f[1])-T0
    jma.append(dict(t=t,lat=float(f[2]),lon=float(f[3]),dep=float(f[4]),M=float(f[6])))
ml=[]
for l in open(D+"/MultipleTemplate/DetectedFinal.dat"):
    if l.startswith("#") or not l.strip(): continue
    f=l.split()
    t=UTCDateTime(f[1].replace("/","-")+"T"+f[2])-T0
    ml.append(dict(no=int(f[0]),t=t,lat=float(f[3]),lon=float(f[4]),dep=float(f[5]),M=float(f[6]),coef=float(f[7]),mad=float(f[8]),ref=f[9]))
def km(a,b,c,d):
    # haversine
    R=6371.0; p1,p2=math.radians(a),math.radians(c); dp=p2-p1; dl=math.radians(d-b)
    h=math.sin(dp/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
    return 2*R*math.asin(math.sqrt(h))
for e in ml:
    e["jma"]=any(abs(e["t"]-j["t"])<=1.0 and km(e["lat"],e["lon"],j["lat"],j["lon"])<=5.0 for j in jma)
print("JMA",len(jma),"ML",len(ml),"jma-flag",sum(e["jma"] for e in ml))
# which JMA events matched, by which ML event
for j in jma:
    m=[e for e in ml if abs(e["t"]-j["t"])<=1.0 and km(e["lat"],e["lon"],j["lat"],j["lon"])<=5.0]
    print(" JMA t %.2f M %.1f -> %s"%(j["t"],j["M"],[(e["no"],round(e["t"]-j["t"],2),e["M"],round(km(e["lat"],e["lon"],j["lat"],j["lon"]),2),e["ref"]) for e in m]))
sys.path.insert(0,"/home/claude/work/r1")
import r1_lib as R
cat=R.read_catalogue()
assert len(cat)==len(ml)
dt=max(abs(a["t"]-b["t"]) for a,b in zip(cat,ml)); dM=max(abs(a["M"]-b["M"]) for a,b in zip(cat,ml))
dj=sum(a["jma"]!=b["jma"] for a,b in zip(cat,ml))
print("max |dt| vs r1_lib %.2e  max|dM| %.2e  jma flag mismatches %d"%(dt,dM,dj))
rj=R.read_jma(); print("read_jma max dt %.2e"%max(abs(a["t"]-b["t"]) for a,b in zip(rj,jma)))
cal=[j for j in jma if 0<=j["t"]<840]; c=(np.mean([j["lat"] for j in cal]),np.mean([j["lon"] for j in cal]))
print("CAL JMA",len(cal),"centre",c,"r1_lib.frame",R.frame())
segs={"CAL":(0,840),"EVAL":(840,1960),"TAIL":(1960,2000)}
for s,(a,b) in segs.items():
    E=[e for e in ml if a<=e["t"]<b]
    print(s,"truth",len(E),"M>=0.5",sum(e["M"]>=0.5 for e in E),"M>=1.0",sum(e["M"]>=1.0 for e in E),"jma",sum(e["jma"] for e in E),"JMA cat",sum(a<=j["t"]<b for j in jma))
print("calwin (CAL, M>=1, t>=6.5):",sum(1 for e in ml if 6.5<=e["t"]<840 and e["M"]>=1.0))
# stations
print("stations:")
for n in R.SUBNETS["H7"]:
    for comp in "ENU":
        tr=read(f"{D}/Trace/20120902/N.{n}.{comp}")[0]; h=tr.stats.sac
        if comp=="U":
            x,y=R.to_xy(h.stla,h.stlo,c)
            print(" %s lat %.4f lon %.4f stel %s stdp %s x %.2f y %.2f r %.2f"%(n,h.stla,h.stlo,h.get("stel"),h.get("stdp"),x,y,math.hypot(x,y)))
        print("   %s.%s start %s npts %d sr %.3f end %s cmpaz %s cmpinc %s kcmpnm %s"%(n,comp,tr.stats.starttime,tr.stats.npts,tr.stats.sampling_rate,tr.stats.endtime,h.get("cmpaz"),h.get("cmpinc"),h.get("kcmpnm")))
# event spread
xy=np.array([R.to_xy(e["lat"],e["lon"],c) for e in ml]); r=np.hypot(xy[:,0],xy[:,1])
print("truth epicentres: x range %.2f %.2f y range %.2f %.2f max r %.2f, depth range %.2f %.2f"%(xy[:,0].min(),xy[:,0].max(),xy[:,1].min(),xy[:,1].max(),r.max(),min(e["dep"] for e in ml),max(e["dep"] for e in ml)))
print("unique locations", len(set((e["lat"],e["lon"],e["dep"]) for e in ml)))
print("M range %.2f %.2f"%(min(e["M"] for e in ml),max(e["M"] for e in ml)))
from collections import Counter
print("refs",Counter(e["ref"] for e in ml).most_common(15))
json.dump(ml,open("/home/claude/work/r1_verif_indep/truth_own.json","w"))
json.dump(jma,open("/home/claude/work/r1_verif_indep/jma_own.json","w"))
