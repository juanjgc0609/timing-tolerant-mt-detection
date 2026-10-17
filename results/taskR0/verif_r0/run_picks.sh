#!/bin/bash
# Re-run frozen r0_picks.py (unchanged) on EVAL for 4 subnets; PhaseNet at PHASENET_HOME inside r0_verif
cd /home/claude/work || exit 1
head -34 r0/R0_FREEZE.txt | grep -v -e gamma_src -e R0_DEV_LOG | sha256sum -c --quiet || { echo "freeze check failed"; exit 1; }
grep gamma_src r0/R0_FREEZE.txt | sha256sum -c --quiet || { echo "gamma check failed"; exit 1; }
export PATH=/home/claude/venv311/bin:$PATH OPENBLAS_NUM_THREADS=1 OPENBLAS_CORETYPE=Prescott OMP_NUM_THREADS=1
export PHASENET_HOME=/home/claude/work/r0_verif/phasenet PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=/home/claude/work/r0_verif/pydeps
cd r0
O=/home/claude/work/r0_verif/picks
for s in N8 N6 N4 U8; do
  [ -f $O/picks_eval_$s.json ] && continue
  echo "START $s $(date -u)" >> $O/picks.log
  nice -n 5 python3 r0_picks.py eval $s $O/picks_eval_$s.json >> $O/picks.log 2>&1
done
echo PICKS_DONE >> $O/picks.log
