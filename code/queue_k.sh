#!/bin/bash
cd /home/claude/work
while [ ! -f done_j2.flag ]; do sleep 30; done
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
run_one() { s=$1; c=$2; j=$3; TASK1K_JITTER=$j python3 task1jk_controls.py --part k --seed $s --scenarios $c --workers 1 --out k_${c}_${s}_j${j}.json > klog_${c}_${s}_j${j}.out 2> klog_${c}_${s}_j${j}.err; }
export -f run_one
printf "%s\n" "2024 unilateral_c40 0" "2024 ladder_4 0" "2024 dispersa_c40 0" "2024 ladder_6 0" "7331 unilateral_c40 0" "7331 ladder_4 0" "7331 dispersa_c40 0" "7331 ladder_6 0" "2024 unilateral_c40 0.2" "2024 ladder_4 0.2" "7331 unilateral_c40 0.2" "7331 ladder_4 0.2" | xargs -P 2 -L 1 bash -c 'run_one $0 $1 $2'
echo DONE > done_k.flag
