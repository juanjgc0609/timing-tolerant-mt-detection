#!/bin/bash
# Task R1 single run (frozen code): PhaseNet and associators on the full record, P1 on CAL and EVAL, hybrid lists,
# analysis of EVAL.  Each step is skipped when its output exists (resumable after a container restart).
cd /home/claude/work/r1 || exit 1
. ./env.sh
O=${1:-/home/claude/work/r1/run}
mkdir -p $O
[ -f $O/picks_all.json ] && grep -q '"done": true' $O/picks_all.json || python3 r1_picks.py eval $O/picks_all.json >> $O/picks.log 2>&1
grep -q '"done": true' $O/picks_all.json || { echo PICKS_FAILED >> $O/queue.log; exit 1; }
for s in cal eval; do
  [ -f $O/p1_$s.json ] || { python3 r1_run.py run $s $O/picks_all.json $O/p1_${s}_chunks 40 >> $O/p1_$s.log 2>&1 &&
                            python3 r1_run.py merge $s $O/p1_${s}_chunks $O/p1_$s.json >> $O/p1_$s.log 2>&1; }
  [ -f $O/hyb_$s.json ] || python3 r1_hybrid.py $s $O/picks_all.json $O/p1_$s.json $O/hyb_$s.json >> $O/hyb.log 2>&1
done
[ -f $O/analysis_eval.txt ] || python3 r1_analyze.py eval $O/hyb_eval.json $O/hyb_cal.json $O/analysis_eval.txt >> $O/queue.log 2>&1
echo R1_DONE $(date -u) >> $O/queue.log
