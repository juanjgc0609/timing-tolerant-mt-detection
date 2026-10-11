#!/bin/bash
# C8 supplement (post hoc): Huber-constant sweep on the DENSE-IID windows, after the main queue of the same lane.
cd /home/claude/work/explore_contrib
export OPENBLAS_CORETYPE=Prescott OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
REG=$1
while ! grep -q ALLDONE out/c8_$REG.log; do sleep 30; done
python c8_run.py DENSE-IID-HSWEEP $REG 60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79 5 out/c8_DENSE-IID-HSWEEP_$REG.json
echo ALLDONE2 $REG $(date)
