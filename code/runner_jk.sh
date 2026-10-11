#!/bin/bash
cd /home/claude/work
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
for s in 2024 7331; do python3 task1jk_controls.py --part j --seed $s --workers 2 2>/dev/null | grep -av -E "Associating|Clustering" > log_j_$s.txt; done
for s in 2024 7331; do python3 task1jk_controls.py --part k --seed $s --workers 2 2>/dev/null | grep -av -E "Associating|Clustering" > log_k_$s.txt; done
echo DONE > done.flag
