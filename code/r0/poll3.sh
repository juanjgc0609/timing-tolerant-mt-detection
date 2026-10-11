#!/bin/bash
cd /home/claude/work/r0
st() { ls run run/parts | wc -l; }
s0=$(st); for i in $(seq 1 23); do [ "$(st)" != "$s0" ] && break; sleep 5; done
echo "$(date +%H:%M) assoc $(grep -cE '^gamma|^pyocto' run/resume.log) | $(grep -E ' eval [0-9]+ 887|Trace|Error' run/resume.log | tail -n 1 | cut -c1-90) | $(ls run/*.json run/parts/*.json 2>/dev/null | xargs -n1 basename | tr '\n' ' ')"
