#!/bin/bash
cd /home/claude/work
while [ ! -f done.flag ]; do sleep 30; done
export OMP_NUM_THREADS=1
for s in 2024 7331; do python3 task1j2_noise.py --seed $s --workers 2 > log_j2_$s.txt 2>&1; done
echo DONE > done_j2.flag
