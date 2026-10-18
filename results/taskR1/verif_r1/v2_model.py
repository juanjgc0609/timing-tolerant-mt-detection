import sys
sys.path[:0]=["/home/claude/work/hyb/pydeps_pyrocko"]
import numpy as np, pandas as pd, pyocto
P=[(0.0,4.800),(1.5,5.220),(2.5,5.450),(4.0,5.720),(5.0,5.820),(15.0,6.220),(22.0,6.544),(34.0,7.188),(41.0,7.494),(44.0,7.598),(48.0,7.688),(52.5,7.750),(65.0,7.845)]
z=[p[0] for p in P if p[0]<40]+[40.0]; vp=np.interp(z,[p[0] for p in P],[p[1] for p in P])
df=pd.DataFrame({"depth":z,"vp":vp,"vs":vp/1.73}); print(df)
pyocto.VelocityModel1D.create_model(df,0.25,80.0,40.0,"/home/claude/work/r1_verif_indep/pyocto_jma2001_own.model")
