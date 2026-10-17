#!/bin/bash
# R0 closure (post hoc, declared in the log before launching): M1 and M3 on EVAL of the four subnetworks with the
# frozen R0 code (sha256sum -c against R0_FREEZE.txt first), run through r0_run_chunks.py (imports r0_run.py unchanged;
# resumable chunks of 50 windows), merged into run_close/det_eval_<s>.json and detM3_eval_<s>.json, then
# r0_analyze.py on these files only (registered secondary contrast M1 - M3, and a re-run check of M1).
cd /home/claude/work || exit 1
head -34 r0/R0_FREEZE.txt | grep -v -e gamma_src -e R0_DEV_LOG | sha256sum -c --quiet || { echo "R0 freeze check failed"; exit 1; }
export PATH=/home/claude/venv311/bin:$PATH OPENBLAS_NUM_THREADS=1 OPENBLAS_CORETYPE=Prescott OMP_NUM_THREADS=1
cd r0 && mkdir -p run_close
one() {  # method subnet output
  [ -f "$3" ] && return
  nice -n 10 python3 r0_run_chunks.py eval $2 $1 run_close/chunks_$1_$2 50 >> run_close/close.log 2>&1 &&
  python3 r0_run_chunks.py merge eval $2 $1 run_close/chunks_$1_$2 "$3" >> run_close/close.log 2>&1
}
laneA() { one M3 N8 run_close/detM3_eval_N8.json; one M3 U8 run_close/detM3_eval_U8.json
          one M1 N8 run_close/det_eval_N8.json; one M1 N6 run_close/det_eval_N6.json; echo LANEA_DONE >> run_close/close.log; }
laneB() { one M3 N6 run_close/detM3_eval_N6.json; one M3 N4 run_close/detM3_eval_N4.json
          one M1 N4 run_close/det_eval_N4.json; one M1 U8 run_close/det_eval_U8.json; echo LANEB_DONE >> run_close/close.log; }
laneA & laneB & wait
echo CLOSE_DONE >> run_close/close.log
