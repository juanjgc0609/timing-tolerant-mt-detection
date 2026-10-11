#!/bin/bash
cd /home/claude/work
while ps -eo args | grep -q "[r]un_queue9"; do sleep 30; done
python chk_nms.py unilateral_c40 pert 2024 100 > log_chk_nms_uni_pert.txt 2>&1 &
python chk_nms.py ladder_6 pert 2024 100 > log_chk_nms_l6_pert.txt 2>&1 &
wait
python chk_nms.py ladder_4 ideal 2024 100 > log_chk_nms_l4.txt 2>&1 &
python chk_nms.py ladder_4 pert 2024 100 > log_chk_nms_l4_pert.txt 2>&1 &
wait
