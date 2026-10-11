#!/bin/bash
# R0 frozen run, resumed after the container restart of 9 Oct 21:07 (-05): same frozen code, outputs per method subset
# in run/parts, merged into run/det_eval_<subnet>.json by merge_parts.py.  Primary contrast inputs (M2, picks) first.
cd /home/claude/work/r0
export OPENBLAS_NUM_THREADS=1 OPENBLAS_CORETYPE=Prescott OMP_NUM_THREADS=1
run() { [ -f "$4" ] || python3 r0_run.py "$1" "$2" "$3" "$4" >> run/resume.log 2>&1; }
pk() { [ -f run/picks_eval_$1.json ] || python3 r0_picks.py eval $1 run/picks_eval_$1.json >> run/resume.log 2>&1; }
p1() {
  for s in N8 N6 N4 U8; do pk $s; done
  run eval U8 M2,M6 run/parts/det_eval_U8_M2M6.json
  run eval U8 M1,M5 run/parts/det_eval_U8_M1M5.json
  run eval N4 M3 run/detM3_eval_N4.json
  run eval U8 M3 run/detM3_eval_U8.json
  echo P1_DONE >> run/resume.log
}
p2() {
  run eval N6 M2,M6 run/parts/det_eval_N6_M2M6.json
  run eval N4 M2,M6 run/parts/det_eval_N4_M2M6.json
  run eval N6 M1,M5 run/parts/det_eval_N6_M1M5.json
  run eval N4 M1,M5 run/parts/det_eval_N4_M1M5.json
  run eval N6 M3 run/detM3_eval_N6.json
  echo P2_DONE >> run/resume.log
}
p1 & p2 & wait
