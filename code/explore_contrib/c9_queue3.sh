#!/bin/bash
# C9 queue 3 (post hoc in-loop polarization selection): waits for queue 2, then one process, sequential.
cd /home/claude/work/explore_contrib
while [ ! -f out/c9_queue2_done.txt ]; do sleep 20; done
J=$(seq -s, 130 149)
for D in DENSE STD DENSE-POL0 DENSE-POL20 DENSE-POL45; do
  for R in thin strong; do
    python3 c9_extra.py $D $R $J 5 out/c9x_${D}_${R}.json >> out/c9_run.log 2>&1
  done
done
date -u > out/c9_queue3_done.txt
