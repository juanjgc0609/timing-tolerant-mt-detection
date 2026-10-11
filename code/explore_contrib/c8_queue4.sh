#!/bin/bash
cd /home/claude/work/explore_contrib
export OPENBLAS_CORETYPE=Prescott OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
REG=$1
while ! grep -q ALLDONE out/c8_${REG}3.log; do sleep 20; done
python c8_run.py DENSE-IID-HWIDE $REG 60,61,62,63,64,65,66,67,68,69 5 out/c8_DENSE-IID-HWIDE_$REG.json
echo ALLDONE4 $REG $(date)
