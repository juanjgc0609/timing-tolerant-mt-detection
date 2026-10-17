#!/bin/bash
# Task 2M: run every (network, cell) of a list as its own call of the frozen runner, so that one network that exhausts
# its catalogue pool (declared as failed) does not stop the others.  Usage: ./t2m_resume_each.sh SEED ARM NETS
SEED=$1; ARM=$2
for j in $(echo $3 | tr ',' ' '); do
  for c in ladder_4 ladder_6 unilateral_c40; do
    ./t2m_resume.sh $SEED $ARM $j $c
  done
done
echo "RESUME_EACH_DONE $SEED $ARM $3 $(date)" >> t2m_out/run_${SEED}_${ARM}_resume.log
