#!/bin/bash
# Re-run frozen M2 on EVAL for 4 subnets via r0_run_chunks.py (unchanged), outputs in r0_verif/m2
cd /home/claude/work || exit 1
head -34 r0/R0_FREEZE.txt | grep -v -e gamma_src -e R0_DEV_LOG | sha256sum -c --quiet || { echo "freeze check failed"; exit 1; }
export PATH=/home/claude/venv311/bin:$PATH OPENBLAS_NUM_THREADS=1 OPENBLAS_CORETYPE=Prescott OMP_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
cd r0
O=/home/claude/work/r0_verif/m2
one() { [ -f $O/det_eval_$1.json ] && return
  python3 r0_run_chunks.py eval $1 M2 $O/chunks_M2_$1 50 >> $O/m2.log 2>&1 &&
  python3 r0_run_chunks.py merge eval $1 M2 $O/chunks_M2_$1 $O/det_eval_$1.json >> $O/m2.log 2>&1; }
if [ "$1" = A ]; then one N8; one N4; echo LANEA_DONE >> $O/m2.log; fi
if [ "$1" = B ]; then one N6; one U8; echo LANEB_DONE >> $O/m2.log; fi
