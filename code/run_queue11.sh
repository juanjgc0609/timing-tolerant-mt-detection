#!/bin/bash
cd /home/claude/work
python task1y_clutter.py --seed 2024 --reps 100 --clutter thin > log1y_2024_thin.txt 2>&1 &
python task1y_clutter.py --seed 7331 --reps 100 --clutter thin > log1y_7331_thin.txt 2>&1 &
wait
python task1y_clutter.py --seed 2024 --reps 100 --clutter strong > log1y_2024_strong.txt 2>&1 &
python task1y_clutter.py --seed 7331 --reps 100 --clutter strong > log1y_7331_strong.txt 2>&1 &
wait
