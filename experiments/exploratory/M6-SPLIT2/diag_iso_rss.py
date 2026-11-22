"""M6-SPLIT2 diagnostic (after the run, seed 999 outside the declared seeds): RSS of the isolated windows in the
two-event trial and in the one-event alternative of each test. Run from the repository root."""
import sys, importlib.util, numpy as np
sys.argv=["run.py","1","1","1"]
spec=importlib.util.spec_from_file_location("ms2","experiments/exploratory/M6-SPLIT2/run.py"); M=importlib.util.module_from_spec(spec); spec.loader.exec_module(M)
import ttmt.v2.m6split as MS
orig_free, orig_iso = MS.free_one_event, MS.iso_rss
log=[]
def iso_rss(iso,G,resp,frames):
    v=orig_iso(iso,G,resp,frames); log.append(v); return v
MS.iso_rss=iso_rss
def free(*a,**k):
    n=len(log); r=orig_free(*a,**k); print("   one-event iso_rss", round(log[-1],1) if len(log)>n else None, "total", round(r,1)); return r
MS.free_one_event=free
for c in range(2):
  srcs,Y,q,iso=M.scene((999,1,c),"C1","iii")
  log.clear()
  a=M.loops("iii",Y,q,-np.inf,iso)
  print("scene",c,"iso_rss calls",[round(x,1) for x in log],"gains",[x["pair_gain"] and round(x["pair_gain"],1) for x in a])
