#!/bin/bash
# C13 post-run queue (one process at a time): waits for the D3c run, then twins of the event windows, analysis, tables
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OPENBLAS_CORETYPE=Prescott
while pgrep -f c13_run_d3c.py > /dev/null; do sleep 3; done
python c13_twin.py > out/c13_twin.log 2>&1
python c13_analyze.py 400 1000 > out/c13_analyze.log 2>&1
python c13_tables.py >> out/c13_analyze.log 2>&1
date > out/c13_queue_done.flag
