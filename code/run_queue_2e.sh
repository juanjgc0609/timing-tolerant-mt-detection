#!/bin/sh
# Tarea 2E: evaluacion en las semillas 2024 y 7331 (30 redes por celda y regimen). Nucleo de OpenBLAS de 2C: Prescott
# en todas las redes salvo las 18 redes de 7331 de 2A que 2A produjo con SkylakeX. Reanudable: las partes escritas en
# task2e_parts/ se saltan.
cd /home/claude/work
OPENBLAS_CORETYPE=Prescott python3 -u task2e_run.py run --workers 2 || exit 1
OPENBLAS_CORETYPE=SkylakeX python3 -u task2e_run.py run --workers 2 || exit 1
echo done > done_2e.flag
