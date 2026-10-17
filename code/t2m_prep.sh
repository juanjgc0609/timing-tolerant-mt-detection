#!/bin/bash
# Task 2M noise preparation (no injection, no detector run).
#  1. day 100723 without QC, 14 stations: reproduction of C16's noise (t2m_validate.py noise: bit-identical).
#  2. station sets of the new days by the rule of t2m_stations.py (QC qc1; outputs t2m_out/stations_<day>.json):
#     101231 -> 10 stations (ATCC dropped), 100210 and 100211 -> 8 stations (MURB dropped), 100212 -> not used (0 slots
#     with 8 stations); 100723 is also evaluated, for the report only.
#  3. noise of the used new days with QC qc1 and the selected station sets; day 100723 with qc1, for the report only.
SP=/tmp/claude-0/-home-claude/dd3b46bc-0f2e-5360-b938-9624bcf101de/scratchpad
cd $SP
# C16 prepared its noise with the default OpenBLAS kernel (DYNAMIC_ARCH, CORETYPE unset); with Prescott the trend fit
# differs in the last bits.  Same environment here, so that day 100723 is bit-identical to c16_noise.npy.
unset OPENBLAS_CORETYPE
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
# Python 3.11, NumPy 1.26.4, SciPy 1.17.1, ObsPy 1.5.1 (the environment of the C16 preparation)
export PATH=/home/claude/venv311/bin:$PATH
W=/home/claude/work
D0=/home/claude/realdata/pympa37_altotiberina_100723
DN=/home/claude/realdata/pympa37_days
O=$W/t2m_out/noise
STEP=${1:-all}
if [ "$STEP" = all ] || [ "$STEP" = c16 ]; then
nice -n 10 python3 -I $W/t2m_data.py 100723 ATBU,ATCA,ATCC,ATFO,ATLO,ATPC,ATPI,ATVO,BADI,FRON,MURB,PARC,PIEI,SSFR $D0 $D0/inv.IV $D0/templates.zmap $O/100723 none $D0/200.100723.cat > $O/100723.log 2>&1
fi
if [ "$STEP" = all ] || [ "$STEP" = stations ]; then
bash $W/t2m_out/run_stations.sh > $W/t2m_out/stations.log 2>&1
fi
if [ "$STEP" = all ] || [ "$STEP" = new ]; then
nice -n 10 python3 -I $W/t2m_data.py 101231 ATBU,ATCA,ATPC,ATPI,ATVO,FRON,MURB,PARC,PIEI,SSFR $DN $D0/inv.IV $D0/templates.zmap $O/101231 qc1 > $O/101231.log 2>&1
nice -n 10 python3 -I $W/t2m_data.py 100210 ATBU,ATCC,ATLO,ATPC,ATVO,BADI,PARC,PIEI $DN $D0/inv.IV $D0/templates.zmap $O/100210 qc1 > $O/100210.log 2>&1
nice -n 10 python3 -I $W/t2m_data.py 100211 ATBU,ATCC,ATLO,ATPC,ATVO,BADI,PARC,PIEI $DN $D0/inv.IV $D0/templates.zmap $O/100211 qc1 > $O/100211.log 2>&1
nice -n 10 python3 -I $W/t2m_data.py 100723 ATBU,ATCA,ATCC,ATFO,ATLO,ATPC,ATPI,ATVO,BADI,FRON,MURB,PARC,PIEI,SSFR $D0 $D0/inv.IV $D0/templates.zmap $O/100723_qc1 qc1 $D0/200.100723.cat > $O/100723_qc1.log 2>&1
echo PREP_DONE > $O/PREP_DONE
fi
