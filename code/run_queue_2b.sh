#!/bin/sh
# Tarea 2B: redes de 2A producidas con el nucleo AVX-512 de OpenBLAS (sesiones reanudadas de la semilla 7331) con
# OPENBLAS_CORETYPE=SkylakeX; el resto con el nucleo por defecto de esta maquina; al final se agregan los JSON.
cd /home/claude/work
OPENBLAS_CORETYPE=SkylakeX python task2b_ablation.py run --workers 2 --no-aggregate --nets 7331_strong_ladder_6_2,7331_strong_ladder_6_3,7331_strong_ladder_6_4,7331_thin_ladder_6_0,7331_thin_ladder_6_1,7331_thin_ladder_6_2,7331_thin_ladder_6_3,7331_thin_ladder_6_4,7331_strong_ladder_4_0,7331_strong_ladder_4_1,7331_strong_ladder_4_2,7331_strong_ladder_4_3,7331_strong_ladder_4_4,7331_thin_ladder_4_0,7331_thin_ladder_4_1,7331_thin_ladder_4_2,7331_thin_ladder_4_3,7331_thin_ladder_4_4 || exit 1
python task2b_ablation.py run --workers 2 || exit 1
echo done > done_2b.flag
