#!/bin/bash
# C9 queue 2 (post hoc worlds): waits for queue 1, then one process, sequential.
cd /home/claude/work/explore_contrib
while [ ! -f out/c9_queue_done.txt ]; do sleep 20; done
J=$(seq -s, 130 149)
for D in DENSE-POL20 DENSE-POL45; do
  for R in thin strong; do
    python3 c9_run.py $D $R $J 5 out/c9_${D}_${R}.json >> out/c9_run.log 2>&1
  done
done
date -u > out/c9_queue2_done.txt
