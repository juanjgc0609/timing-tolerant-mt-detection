#!/bin/bash
# C8 queue: lane = regime (thin | strong); designs in priority order.
cd /home/claude/work/explore_contrib
export OPENBLAS_CORETYPE=Prescott OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
REG=$1
M="60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79"
S="80,81,82,83,84,85,86,87,88,89"
python c8_run.py DENSE-IID $REG $M 5 out/c8_DENSE-IID_$REG.json
python c8_run.py DENSE-PERS $REG $M 4 out/c8_DENSE-PERS_$REG.json
python c8_run.py STD-IID $REG $M 4 out/c8_STD-IID_$REG.json
for D in DENSE-IID-S80 DENSE-IID-S25 DENSE-PERS-S80 DENSE-PERS-PSSAME; do
  python c8_run.py $D $REG $S 4 out/c8_${D}_$REG.json
done
python c8_run.py STD-PERS $REG 60,61,62,63,64,65,66,67,68,69 4 out/c8_STD-PERS_$REG.json
echo ALLDONE $REG $(date)
