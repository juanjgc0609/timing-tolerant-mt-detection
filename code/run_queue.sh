#!/bin/bash
# cola: cuando terminan las corridas de 1s, 1r (consistente) y baselines con 200 ventanas
cd /home/claude/work
while ps -eo args | grep -q "^python task1s_gamma_pol"; do sleep 30; done
python task1r_fair.py --seed 2024 --reps 200 --workers 1 --cond consistente > log1r_2024_r200.txt 2>&1 &
python task1r_fair.py --seed 7331 --reps 200 --workers 1 --cond consistente > log1r_7331_r200.txt 2>&1 &
wait
python task1r_baselines.py --seed 2024 --reps 200 --workers 1 > log1rb_2024_r200.txt 2>&1 &
python task1r_baselines.py --seed 7331 --reps 200 --workers 1 > log1rb_7331_r200.txt 2>&1 &
wait
python task1t_wave.py --seed 2024 --reps 200 > log1t_2024_r200.txt 2>&1 &
python task1t_wave.py --seed 7331 --reps 200 > log1t_7331_r200.txt 2>&1 &
wait
