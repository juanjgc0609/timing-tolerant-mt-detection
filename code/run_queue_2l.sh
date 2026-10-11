#!/bin/bash
# Tarea 2L: corrida de las semillas 7171 y 9393 (lanzada despues de congelar task2l_code_hashes.txt); reanudable.
cd /home/claude/work
export OPENBLAS_CORETYPE=Prescott OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
sha256sum -c --quiet task2l_code_hashes.txt || { echo "hashes distintos" >> task2l_FAIL.txt; exit 1; }
python3 task2l_run.py run --seeds 7171,9393 --workers 2 >> task2l_run.log 2>&1
date -u > task2l_run_done.flag
