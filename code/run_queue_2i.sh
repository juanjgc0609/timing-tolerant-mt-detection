#!/bin/bash
# Tarea 2I: primario (las dos semillas) y despues los brazos secundarios; reanudable (partes por red).
cd /home/claude/work
export OPENBLAS_CORETYPE=Prescott OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 task2i_run.py run --seeds 5151,6262 --arms prim --workers 2 >> log2i_run_stdout.txt 2>&1
python3 task2i_run.py run --seeds 5151,6262 --arms thr,s25,s80 --workers 2 >> log2i_run_stdout.txt 2>&1
date -u > done_2i.flag
