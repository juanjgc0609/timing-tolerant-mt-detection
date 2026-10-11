#!/bin/bash
# C11 evaluation runs (launched after C11_PLAN.txt was hashed): one process per regime, networks 210-239 in order
cd /home/claude/work/explore_contrib
export OPENBLAS_CORETYPE=Prescott OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
R=$1
python3 c11_run.py $R 210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239 5 out/c11_eval_$R.json >> out/c11_eval_$R.log 2>&1
date -u > out/c11_evaldone_$R.flag
