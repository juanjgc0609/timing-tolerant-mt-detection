#!/bin/bash
cd /home/claude/work/r0
st() { echo "$(grep -cE ' eval [0-9]+ 887|DONE|^phasenet' run/resume.log)$(ls run run/parts | wc -l)"; }
s0=$(st); for i in $(seq 1 22); do [ "$(st)" != "$s0" ] && break; sleep 5; done
echo "$(date +%H:%M) $(grep -E ' eval [0-9]+ 887|DONE|Trace|Error' run/resume.log | tail -n 2 | cut -c1-100 | tr '\n' '|') files: $(ls run/*.json run/parts/*.json 2>/dev/null | xargs -n1 basename | tr '\n' ' ')"
