#!/bin/bash
cd /home/claude/work
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for s in 2024 7331; do python3 task1o_consistent.py --seed $s --workers 2 > log_1o_$s.txt 2> err_1o_$s.txt; done
echo DONE > done_1o.flag
