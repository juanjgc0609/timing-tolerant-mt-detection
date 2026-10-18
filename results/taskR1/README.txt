Task R1 - hybrid of pick-free moment-tensor candidates and PhaseNet association on the Hi-net swarm of
Geolab-USTC/MatchLocate2 (2 September 2012, 03:20:00-03:53:20 UTC). Appendix Y.33 of the evidence volume and
Section 6.6 of docs/thesis/Thesis.md.

log_entries/   registration (r1_reg_*), amendment 1 (amend_*), result (result_*) and post hoc entry (posthoc_*), Spanish
               and English, as appended to the sealed log (preregistration/)
dev/           development readings on CAL (in-sample), before the freeze; terms_check.py
sim_r1_rules_out.txt   simulation of the decision rules (code/hyb/sim_r1_rules.py)
run/           outputs of the single frozen run (code/r1/r1_queue.sh): picks_all.json, p1_cal.json, p1_eval.json and
               their chunk folders, hyb_cal.json, hyb_eval.json, analysis_eval.txt, logs
verif_r1/      independent re-derivation by a separate AI session (r1_verification.txt, scripts v*.py)
posthoc/       timing residuals of the CAL readings (code/r1/posthoc_timing_cal.py), declared post hoc

Reproduce the analysis:  python code/r1/r1_analyze.py eval results/taskR1/run/hyb_eval.json \
                            results/taskR1/run/hyb_cal.json OUT.txt
Freeze list: code/r1/FREEZE_R1.sha256 (its own SHA-256 is 1669903f22b1e6d4f1b9656ffdc78f5c7a3bc9b4ed4145ce52e41cefcd28a2a3).
