#!/bin/sh
# Tarea 2D: parte (ii) (semilla 999, 200 ventanas vacias por celda y regimen; Prescott) y parte (i) (redes de 2A;
# SkylakeX en las 18 redes de 7331 que 2A produjo con el, Prescott en las demas). Reanudable: las partes escritas en
# task2d_parts/ se saltan.
cd /home/claude/work
OPENBLAS_CORETYPE=Prescott python -u task2d_run.py run --part 2 --workers 2 || exit 1
OPENBLAS_CORETYPE=SkylakeX python -u task2d_run.py run --part 1 --workers 2 || exit 1
OPENBLAS_CORETYPE=Prescott python -u task2d_run.py run --part 1 --workers 2 || exit 1
echo done > done_2d.flag
