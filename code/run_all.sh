#!/bin/bash
cd /home/claude/work
for s in 2024 7331; do python3 task1jk_controls.py --part j --seed $s --workers 2 > log_j_$s.txt 2>&1; done
for s in 2024 7331; do python3 task1jk_controls.py --part k --seed $s --workers 2 > log_k_$s.txt 2>&1; done
echo DONE > done.flag
