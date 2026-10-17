#!/bin/bash
# POST HOC: after m56 lane B, M2 on N8 with +SC from a +-2 s search (verif_widesc.py)
export PATH=/home/claude/venv311/bin:$PATH OPENBLAS_NUM_THREADS=1 OPENBLAS_CORETYPE=Prescott OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
V=/home/claude/work/r0_verif
until grep -q LANEB_DONE $V/m56/m56.log 2>/dev/null; do sleep 20; done
cd /home/claude/work/r0
[ -f $V/widesc/det_eval_N8.json ] || { python3 $V/verif_widesc.py N8 $V/widesc/chunks_N8 50 >> $V/widesc/widesc.log 2>&1 &&
  python3 r0_run_chunks.py merge eval N8 M2 $V/widesc/chunks_N8 $V/widesc/det_eval_N8.json >> $V/widesc/widesc.log 2>&1; }
echo WIDESC_DONE >> $V/widesc/widesc.log
