#!/bin/bash
export PATH=/home/claude/venv311/bin:$PATH OPENBLAS_NUM_THREADS=1 OPENBLAS_CORETYPE=Prescott OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
V=/home/claude/work/r0_verif; s=$1
cd /home/claude/work/r0
[ -f $V/wides2/det_eval_$s.json ] || { python3 $V/verif_wides2.py $s $V/wides2/chunks_$s 50 >> $V/wides2/wides2.log 2>&1 &&
  python3 r0_run_chunks.py merge eval $s M2 $V/wides2/chunks_$s $V/wides2/det_eval_$s.json >> $V/wides2/wides2.log 2>&1; }
echo WIDES2_DONE_$s >> $V/wides2/wides2.log
