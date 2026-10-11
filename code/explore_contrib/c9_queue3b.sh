#!/bin/bash
# C9 queue 3b (replaces queues 3 and 4, cut for time: extras only for DENSE, STD and DENSE-POL45), one process.
cd /home/claude/work/explore_contrib
J=$(seq -s, 130 149)
for X in "DENSE strong" "STD thin" "STD strong" "DENSE-POL45 thin" "DENSE-POL45 strong"; do
  set -- $X
  python3 c9_extra.py $1 $2 $J 5 out/c9x_${1}_${2}.json >> out/c9_run.log 2>&1
done
for D in DENSE STD; do NBOOT=1000 python3 c9_analyze.py $D > out/c9_analysis_$D.txt 2> out/c9_analysis_$D.err; done
for D in DENSE-BOTH0 DENSE-POL0 DENSE-AMP0 STD-BOTH0 DENSE-POL20 DENSE-POL45; do
  WORLD=1 NBOOT=500 python3 c9_analyze.py $D > out/c9_analysis_$D.txt 2> out/c9_analysis_$D.err
done
python3 c9_run.py validate > out/c9_validation_end.txt 2>&1
date -u > out/c9_queue4_done.txt
