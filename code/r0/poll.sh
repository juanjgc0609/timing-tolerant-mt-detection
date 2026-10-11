#!/bin/bash
# block up to 115 s until a lane log or the C16 part count changes; print a one-line status
cd /home/claude/work/r0
s0="$(wc -l run/laneA.log run/laneB.log 2>/dev/null | tail -1)$(ls /home/claude/work/explore_contrib/out/c16_eval_parts 2>/dev/null | wc -l)"
for i in $(seq 1 23); do
  s1="$(wc -l run/laneA.log run/laneB.log 2>/dev/null | tail -1)$(ls /home/claude/work/explore_contrib/out/c16_eval_parts 2>/dev/null | wc -l)"
  [ "$s1" != "$s0" ] && break
  sleep 5
done
echo "$(date +%H:%M) A: $(tail -n 1 run/laneA.log | cut -c1-90) | B: $(tail -n 1 run/laneB.log | cut -c1-60) | C16 parts $(ls /home/claude/work/explore_contrib/out/c16_eval_parts | wc -l) | files: $(ls run | grep json | tr '\n' ' ')"
