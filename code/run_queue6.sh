#!/bin/bash
cd /home/claude/work
while ps -eo args | grep -q "^python task1t_wave.py --seed 2024"; do sleep 20; done
python task1t_wave.py --seed 7331 --reps 200 > log1t_7331_exact.txt 2>&1
python task1u_hetero.py --seed 2024 --reps 100 --wave-only > log1u_2024_exact.txt 2>&1
python task1u_hetero.py --seed 7331 --reps 100 --wave-only > log1u_7331_exact.txt 2>&1
