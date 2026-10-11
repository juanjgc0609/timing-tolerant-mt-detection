#!/bin/bash
cd /home/claude/work
python task1x_mtglrt.py --seed 2024 --reps 100 --cond ideal > log1x_2024_ideal.txt 2>&1 &
python task1x_mtglrt.py --seed 7331 --reps 100 --cond ideal > log1x_7331_ideal.txt 2>&1 &
wait
python task1x_mtglrt.py --seed 2024 --reps 100 --cond pert > log1x_2024_pert.txt 2>&1 &
python task1x_mtglrt.py --seed 7331 --reps 100 --cond pert > log1x_7331_pert.txt 2>&1 &
wait
