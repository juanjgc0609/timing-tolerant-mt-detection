Results of Tasks 1r, 1s, 1t and 1u (Appendix Y.6 to Y.9 of the proposal).

Decision rules and dated amendments: PREREGISTRO_1j_1k.txt (section 1r onwards).

1r  fair repetition of 1o with negative control PHYS-SH
    task1r_seed{2024,7331}.json                60 windows, registered sweep, conditions consistente and exacto
    task1r_seed{s}_r200.json                   200 windows, condition consistente, registered sweep
    task1r_seed{s}_r200_exact.json             200 windows, exact FROC curve (primary reading)
    task1r_baselines_seed{s}.json              GaMMA (without seed, see E22) and PyOcto, 60 windows
    task1r_baselines_seed{s}_r200_pyocto.json  PyOcto, 200 windows
    analyze_1r_out.txt                         60 windows; analyze_1r200_{exact,coarse}_out.txt: 200 windows
1s  polarisation as a filter on the GaMMA events (GaMMA with seed; its GAMMA column is the
    reference GaMMA of 200 windows for 1r and 1t)
    task1s_seed{s}_r200.json, analyze_1s_r200_out.txt
    log_check_gamma.txt                        GaMMA without seed changes 15 of 80 calls; with seed, 0
1t  pick-free waveform migration (WAVE-E, WAVE-P, WAVE-PSH)
    task1t_seed{s}_r60.json                    reproduction of WAVE-P of 1r (identical)
    task1t_seed{s}_r200.json                   registered sweep; task1t_seed{s}_r200_exact.json: exact curve
    analyze_1t_r*_out.txt
1u  unmodelled medium: delays N(0, 0.15 s) and ~10 degrees of polarisation deviation
    task1u_seed{s}_r100.json                   WAVE-* (registered sweep) and GaMMA with seed
    task1u_seed{s}_r100_wave_exact.json        WAVE-* with exact curve, same windows
    analyze_1u_out.txt
tables_y.md    tables Y.6 to Y.10 generated from the JSON files by make_tables_y.py
log*.txt       reports printed by each run

Reproduction (from the scripts folder, with gamma_src/ in the same folder):
  python task1r_fair.py --seed 2024 --reps 200 --workers 1 --cond consistente --exact
  python task1s_gamma_pol.py --seed 2024 --reps 200 --workers 1
  python task1r_baselines.py --seed 2024 --reps 200 --workers 1 --methods PYOCTO
  python task1t_wave.py --seed 2024 --reps 200
  python task1u_hetero.py --seed 2024 --reps 100            (GaMMA and WAVE-* with registered sweep)
  python task1u_hetero.py --seed 2024 --reps 100 --wave-only
  python check_gamma_det.py
  python analyze_1r200.py ; python analyze_1s.py 200 ; python analyze_1t.py 200 ; python analyze_1u.py
  python make_tables_y.py
Times: 15 to 40 minutes per cell on one core (1r and 1s), 3 to 10 minutes (1t and 1u without GaMMA).
task0_common.generate now accepts statics and dir_noise (1u); with both set to None it reproduces the previous
generator bit for bit (verified).
GaMMA from its source code (github.com/AI4EPS/GaMMA, revision 80394dd), with np.random.seed(42)
before each call; PyOcto 0.2.0.

1w  double-couple sources and controls that take the geometry of another node (Appendix Y.10)
    task1w_seed{s}_r120.json                         KIN-LR, PHYS-LR, PHYS-NODE, WAVE-E, WAVE-P, WAVE-PT,
                                                     WAVE-PT-NODE, GaMMA; 60 calibration + 120 evaluation
    task1w_seed{s}_r120_ev9000_m131_wave.json        replication with 120 new windows (waveform and GaMMA)
    analyze_1w_out.txt, analyze_1w_ev9000_m131_wave_out.txt, table_y13.md
    log_chk_gamma_loc.txt   picks detected per source and GaMMA failures by distance (6 stations)
    log_chk_gamma_tr.txt    GaMMA with the origin-time filter widened to [0, 15] s
Reproduction:
  python task1w_mech.py --seed 2024 --reps 120 --workers 1
  python task1w_mech.py --seed 2024 --reps 120 --workers 1 --wave-only --evbase 9000 --moff 131
  python analyze_1w.py ; python analyze_1w.py _ev9000_m131_wave
Note: the tables of make_tables_y.py re-choose the operating point in each bootstrap replicate;
the files analyze_1r200/1t/1u_out.txt keep the intervals with the fixed point.

1x  coherent moment-tensor detector and Frank-Wolfe with deflation (Appendix Y.11)
    task1x_seed{s}_{ideal,pert}_r100.json, analyze_1x_out.txt, log_chk_nms_*.txt
    The GAMMA column of 1x uses the clutter of E26 (undetectable picks delivered to GaMMA): it is not a
    valid comparison.
1y  consistent clutter, regimes thin and strong, condition pert (Appendix Y.12)
    task1y_seed{s}_{thin,strong}_r100.json, analyze_1y_out.txt
Reproduction:
  python cal_capv.py
  python task1x_mtglrt.py --seed 2024 --reps 100 --cond ideal     (and --cond pert)
  python task1y_clutter.py --seed 2024 --reps 100 --clutter thin   (and --clutter strong)
  python chk_nms.py unilateral_c40 ideal 2024 100
  python analyze_1x.py ; python analyze_1y.py
