#!/bin/bash
cd /home/claude/work
while ps -eo args | grep -q "^python task1r_fair"; do sleep 30; done
python task1r_baselines.py --seed 2024 --reps 200 --workers 1 --methods PYOCTO > log1rb_2024_r200.txt 2>&1 &
python task1r_baselines.py --seed 7331 --reps 200 --workers 1 --methods PYOCTO > log1rb_7331_r200.txt 2>&1 &
wait
