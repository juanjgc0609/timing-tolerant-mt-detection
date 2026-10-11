#!/bin/sh
cd /home/claude/work
OPENBLAS_CORETYPE=SkylakeX python -u task2d_explore.py e1 --workers 2 || exit 1
OPENBLAS_CORETYPE=Prescott python -u task2d_explore.py e1 --workers 2 || exit 1
OPENBLAS_CORETYPE=Prescott python -u task2d_explore.py e2 --workers 2 || exit 1
echo done > done_2d_explore.flag
