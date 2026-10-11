#!/bin/bash
cd /home/claude/work
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1
for s in 2024 7331; do python3 task1r_fair.py --seed $s --workers 2 > log_1r_$s.txt 2> err_1r_$s.txt; done
for s in 2024 7331; do python3 task1r_baselines.py --seed $s --workers 2 2>/dev/null | grep -a hecho > log_1rb_$s.txt; done
echo DONE > done_1r.flag
