#!/bin/bash
# C16 development: REAL arm on networks 2400-2401 (6 windows per design) and the validation of the injection code.
cd /home/claude/work/explore_contrib
export OPENBLAS_NUM_THREADS=1 OPENBLAS_CORETYPE=Prescott OMP_NUM_THREADS=1
mkdir -p out/c16_dev_parts
setsid nohup python3 c16_run.py REAL 2400,2401 6 out/c16_dev_parts > out/c16_dev_REAL.log 2>&1 < /dev/null &
setsid nohup python3 c16_validate.py > out/c16_validation.txt 2> out/c16_validation.err < /dev/null &
