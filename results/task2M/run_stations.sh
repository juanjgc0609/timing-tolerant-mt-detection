#!/bin/bash
# Task 2M station-set rule on every candidate day (data quality only; no injection, no detector run).
# Writes ${OUT:-/home/claude/work/t2m_out}/stations_<day>.json.
W=/home/claude/work
D0=/home/claude/realdata/pympa37_altotiberina_100723
DN=/home/claude/realdata/pympa37_days
O=${OUT:-$W/t2m_out}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python3 -I $W/t2m_stations.py 101231 ATBU,ATCA,ATCC,ATPC,ATPI,ATVO,FRON,MURB,PARC,PIEI,SSFR $DN $D0/templates.zmap $O/stations_101231.json
python3 -I $W/t2m_stations.py 100210 ATBU,ATCC,ATLO,ATPC,ATVO,BADI,MURB,PARC,PIEI $DN $D0/templates.zmap $O/stations_100210.json
python3 -I $W/t2m_stations.py 100211 ATBU,ATCC,ATLO,ATPC,ATVO,BADI,MURB,PARC,PIEI $DN $D0/templates.zmap $O/stations_100211.json
python3 -I $W/t2m_stations.py 100212 ATCC,ATFO,ATLO,ATVO,BADI,MURB,PARC,PIEI $DN $D0/templates.zmap $O/stations_100212.json
python3 -I $W/t2m_stations.py 100723 ATBU,ATCA,ATCC,ATFO,ATLO,ATPC,ATPI,ATVO,BADI,FRON,MURB,PARC,PIEI,SSFR $D0 $D0/templates.zmap $O/stations_100723.json $D0/200.100723.cat
