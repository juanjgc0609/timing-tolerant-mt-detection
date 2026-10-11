#!/bin/bash
cd /home/claude/work
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1
for s in 2024 7331; do python3 task1o_baselines.py --seed $s --workers 2 > log_1ob_$s.txt 2> err_1ob_$s.txt; done
echo DONE > done_1ob.flag
