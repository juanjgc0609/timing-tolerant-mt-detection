#!/bin/bash
# Tarea 2A: evaluacion (semillas 2024 y 7331, ambos regimenes, 2 procesos). Los puntos de 2A-e ya estan
# congelados en task2a_frozen_2Ae.json (semilla 999).
cd /home/claude/work
export TASK2A_SCRATCH=/tmp/claude-0/-home-claude/d3c612af-4e2f-52d9-b554-c7fa18d867f4/scratchpad/t2a
for s in 2024 7331; do
  python task2a_confirm.py run --seed $s --workers 2 > log2a_seed$s.txt 2> err2a_seed$s.txt || { echo "fallo semilla $s" > fail2a_$s.flag; exit 1; }
done
date > done_2a.flag
