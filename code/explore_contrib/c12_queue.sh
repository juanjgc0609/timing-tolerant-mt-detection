#!/bin/bash
cd /home/claude/work/explore_contrib
export OPENBLAS_NUM_THREADS=1 OPENBLAS_CORETYPE=Prescott OMP_NUM_THREADS=1
python c12_run.py eval 1200,1201,1202,1203,1204,1205,1206,1207,1208,1209,1210,1211 20
python c12_run.py nulln 1226,1227,1228,1229 50
python c12_run.py nullc 1220,1221,1222,1223,1224,1225 40
python c12_run.py eval 1212,1213,1214,1215 20
date > out/c12_done.flag
