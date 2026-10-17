#!/bin/bash
# After the M2 lane finishes: frozen M5 and M6 on EVAL (r0_run_chunks.py, unchanged), outputs in r0_verif/m56;
# lane A then runs the POST HOC Z-flip sensitivity of M2 on N8 (verif_zflip.py).
cd /home/claude/work || exit 1
head -34 r0/R0_FREEZE.txt | grep -v -e gamma_src -e R0_DEV_LOG | sha256sum -c --quiet || { echo "freeze check failed"; exit 1; }
export PATH=/home/claude/venv311/bin:$PATH OPENBLAS_NUM_THREADS=1 OPENBLAS_CORETYPE=Prescott OMP_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
V=/home/claude/work/r0_verif; O=$V/m56
until grep -q LANE$1_DONE $V/m2/m2.log; do sleep 20; done
cd r0
one() { [ -f $O/det_eval_$1.json ] && return
  python3 r0_run_chunks.py eval $1 M5,M6 $O/chunks_M56_$1 50 >> $O/m56.log 2>&1 &&
  python3 r0_run_chunks.py merge eval $1 M5,M6 $O/chunks_M56_$1 $O/det_eval_$1.json >> $O/m56.log 2>&1; }
if [ "$1" = A ]; then one N8; one N4; echo LANEA_DONE >> $O/m56.log
  mkdir -p $V/zflip
  [ -f $V/zflip/det_eval_N8.json ] || { python3 $V/verif_zflip.py N8 $V/zflip/chunks_N8 50 >> $V/zflip/zflip.log 2>&1 &&
    python3 r0_run_chunks.py merge eval N8 M2 $V/zflip/chunks_N8 $V/zflip/det_eval_N8.json >> $V/zflip/zflip.log 2>&1; }
  echo ZFLIP_DONE >> $V/zflip/zflip.log; fi
if [ "$1" = B ]; then one N6; one U8; echo LANEB_DONE >> $O/m56.log; fi
