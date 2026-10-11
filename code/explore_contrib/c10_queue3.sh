#!/bin/bash
# C10 post hoc: extra (true statics + catalogue level), then the first-iteration diagnostic on STD windows; chained
# after the sweep of the same regime.  (The first version of this queue, diagnostic on STD and DENSE, was stopped
# while still waiting, before it had run anything.)
cd /home/claude/work/explore_contrib
export OPENBLAS_CORETYPE=Prescott OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
R=$1
while [ ! -f out/c10_sweepdone_$R.flag ]; do sleep 20; done
python3 c10_extra.py DENSE $R 160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179 5 out/c10x_DENSE_$R.json >> out/c10_extra_$R.log 2>&1
python3 c10_extra.py STD $R 160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179 5 out/c10x_STD_$R.json >> out/c10_extra_$R.log 2>&1
python3 c10_diag.py STD $R 160,161,162,163,164,165,166,167,168,169 5 out/c10d_STD_$R.json >> out/c10_diag_$R.log 2>&1
date -u > out/c10_diagdone_$R.flag
