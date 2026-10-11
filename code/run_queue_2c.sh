#!/bin/sh
# Tarea 2C: parte (i) (redes de 2A; nucleo SkylakeX en las 18 redes de 7331 que 2A produjo con el, Prescott en las
# demas) y parte (ii) (25 redes nuevas por celda, semilla y regimen; Prescott). Reanudable: las partes escritas en
# task2c_parts/ se saltan.
cd /home/claude/work
OPENBLAS_CORETYPE=SkylakeX python -u task2c_run.py run --part 1 --workers 2 || exit 1
OPENBLAS_CORETYPE=Prescott python -u task2c_run.py run --part 1 --workers 2 || exit 1
OPENBLAS_CORETYPE=Prescott python -u task2c_run.py run --part 2 --workers 2 || exit 1
echo done > done_2c.flag
