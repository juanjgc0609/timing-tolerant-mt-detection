#!/bin/bash
cd /home/claude/work/verif_2i
for a in "prim 6262 thin unilateral_c40 45" "prim 6262 thin ladder_6 91" "prim 5151 thin ladder_6 95" "prim 6262 strong ladder_4 98" "prim 5151 strong unilateral_c40 109" "prim 6262 strong ladder_6 98"; do
  nice -n 5 python3 regen_verif.py $a >> regen_out_verif.txt 2>&1
done
echo DONE >> regen_out_verif.txt
