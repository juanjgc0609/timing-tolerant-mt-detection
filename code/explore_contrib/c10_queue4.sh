#!/bin/bash
# C10 post hoc: SC-R corrections in the MT-FW family, chained after queue3 of the same regime
cd /home/claude/work/explore_contrib
export OPENBLAS_CORETYPE=Prescott OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
R=$1
while [ ! -f out/c10_diagdone_$R.flag ]; do sleep 20; done
python3 c10_screloc_run.py DENSE $R 160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179 5 out/c10r_DENSE_$R.json >> out/c10_scr_$R.log 2>&1
python3 c10_screloc_run.py STD $R 160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179 5 out/c10r_STD_$R.json >> out/c10_scr_$R.log 2>&1
date -u > out/c10_scrdone_$R.flag
