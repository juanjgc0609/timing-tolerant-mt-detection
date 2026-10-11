#!/bin/bash
# R0 frozen run.  Two lanes; finished outputs are skipped.  Usage: bash queue_eval.sh
cd /home/claude/work/r0
export OPENBLAS_NUM_THREADS=1 OPENBLAS_CORETYPE=Prescott OMP_NUM_THREADS=1
mkdir -p run
laneA() {
  for s in N8 N6 N4 U8; do
    [ -f run/det_eval_$s.json ] || python3 r0_run.py eval $s M1,M2,M5,M6 run/det_eval_$s.json >> run/laneA.log 2>&1
    [ -f run/picks_eval_$s.json ] || python3 r0_picks.py eval $s run/picks_eval_$s.json >> run/laneA.log 2>&1
  done
  for s in N8 N6 N4 U8; do
    [ -f run/picks_cal_$s.json ] || python3 r0_picks.py cal $s run/picks_cal_$s.json >> run/laneA.log 2>&1
    [ -f run/det_cal_$s.json ] || python3 r0_run.py cal $s M1,M2,M5,M6 run/det_cal_$s.json >> run/laneA.log 2>&1
  done
  echo LANEA_DONE >> run/laneA.log
}
laneB() {
  for s in N8 N6 N4 U8; do
    [ -f run/detM3_eval_$s.json ] || python3 r0_run.py eval $s M3 run/detM3_eval_$s.json >> run/laneB.log 2>&1
  done
  for s in N8 N6 N4 U8; do
    [ -f run/detM3_cal_$s.json ] || python3 r0_run.py cal $s M3 run/detM3_cal_$s.json >> run/laneB.log 2>&1
  done
  echo LANEB_DONE >> run/laneB.log
}
laneA &
laneB &
wait
