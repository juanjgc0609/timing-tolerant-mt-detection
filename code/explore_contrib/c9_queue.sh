#!/bin/bash
# C9 P4 queue: one process, sequential.  Seed 999, networks 130-149, 5 windows per network.
cd /home/claude/work/explore_contrib
J=$(seq -s, 130 149)
for D in DENSE STD DENSE-BOTH0 DENSE-POL0 DENSE-AMP0 STD-BOTH0; do
  for R in thin strong; do
    python3 c9_run.py $D $R $J 5 out/c9_${D}_${R}.json >> out/c9_run.log 2>&1
  done
done
date -u > out/c9_queue_done.txt
