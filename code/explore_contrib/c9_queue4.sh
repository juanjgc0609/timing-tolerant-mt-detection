#!/bin/bash
# C9 queue 4: analyses, after queue 3; one process, sequential.
cd /home/claude/work/explore_contrib
while [ ! -f out/c9_queue3_done.txt ]; do sleep 20; done
for D in DENSE STD; do NBOOT=1000 python3 c9_analyze.py $D > out/c9_analysis_$D.txt 2> out/c9_analysis_$D.err; done
for D in DENSE-BOTH0 DENSE-POL0 DENSE-AMP0 STD-BOTH0 DENSE-POL20 DENSE-POL45; do
  WORLD=1 NBOOT=500 python3 c9_analyze.py $D > out/c9_analysis_$D.txt 2> out/c9_analysis_$D.err
done
python3 c9_run.py validate > out/c9_validation_end.txt 2>&1
date -u > out/c9_queue4_done.txt
