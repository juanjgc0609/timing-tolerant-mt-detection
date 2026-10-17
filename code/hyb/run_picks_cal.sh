#!/bin/bash
# PhaseNet + GaMMA + PyOcto (frozen r0_picks.py, unchanged) on the Ridgecrest CAL segment, for hybrid development.
cd /home/claude/work/r0 || exit 1
export PATH=/home/claude/venv311/bin:$PATH OPENBLAS_NUM_THREADS=1 OPENBLAS_CORETYPE=Prescott OMP_NUM_THREADS=1
export PHASENET_HOME=/home/claude/work/r0_verif/phasenet PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=/home/claude/work/r0_verif/pydeps
O=/home/claude/work/hyb
for s in $@; do
  [ -f $O/picks_cal_$s.json ] || nice -n 5 python3 r0_picks.py cal $s $O/picks_cal_$s.json >> $O/picks_cal.log 2>&1
done
echo PICKS_CAL_DONE $@ >> $O/picks_cal.log
