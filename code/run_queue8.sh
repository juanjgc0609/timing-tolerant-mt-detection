#!/bin/bash
cd /home/claude/work
while ps -eo args | grep -q "^python task1w_mech.py --seed 7331 --reps 120 --workers 1$"; do sleep 20; done
python task1w_mech.py --seed 7331 --reps 120 --workers 1 --wave-only --evbase 9000 --moff 131 > log1wR_7331.txt 2>&1
