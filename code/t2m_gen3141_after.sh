#!/bin/bash
# after the 3141 REAL resumption: the descriptive GEN arm of 3141 (networks 2900-2907, finished parts skipped), only if
# the REAL arm ended within 3 h of the first launch (registered budget rule).
cd /home/claude/work
while ! grep -q "RESUME_EACH_DONE 3141 REAL" t2m_out/run_3141_REAL_resume.log; do sleep 20; done
if [ $(( $(date +%s) - $(cat t2m_out/queue_start) )) -le 10800 ]; then
  ./t2m_resume_each.sh 3141 GEN 2900,2901,2902,2903,2904,2905,2906,2907
else
  echo "GEN_SKIPPED_BUDGET 3141 $(date)" >> t2m_out/run_3141_GEN_resume.log
fi
