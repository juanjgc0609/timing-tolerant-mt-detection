Results of Tasks 1j, 1j-b and 1k (Appendix Y of the proposal).

task1j_seed{2024,7331}.json    1j: tie control, protocol P3, 6 cells, 30+30 realizations.
task1j2_seed{2024,7331}.json   1j-b: pick time error and velocity error, 9 cells.
k_<cell>_<seed>_j<error>.json   1k: 0-2 sources per window, absolute threshold, 60+60 windows;
                               j0 = exact picks, j0.2 = Gaussian error of 0.2 s.
PREREGISTRO_1j_1k.txt          decision rules and dated amendments.
log_*.txt, klog_*.txt          reports printed by each run.

Reproduction: python task1jk_controls.py --part j --seed 2024 --workers 2
              TASK1K_JITTER=0.2 python task1jk_controls.py --part k --seed 2024 --scenarios ladder_4 --workers 1
              python task1j2_noise.py --seed 2024 --workers 2
GaMMA was run from its source code (github.com/AI4EPS/GaMMA, revision 80394dd), with numba and tqdm;
PyOcto 0.2.0. Both can only improve with a larger hyperparameter grid.

Version note: part j was run before adding the PHYS-K method and the correction for
stations without picks (which only occurs in empty windows of part k). Running part j today adds
the PHYS-K column and does not change the others.

Task 1o (Appendix Y.4-Y.6):
task1o_seed{2024,7331}.json           consistent generator; own methods; endpoint per budget
task1o_baselines_seed{2024,7331}.json GaMMA and PyOcto tuned on the same windows (7331: consistent cells only)
log_1o_2024_recalltarget.txt          first run with recall target 0.90 (amended, see PREREGISTRO)
analyze_1o_out.txt                    tables Y.5 and Y.6
Reproduction: python task1o_consistent.py --seed 2024 ; python task1o_baselines.py --seed 2024 ; python analyze_1o.py
Certified pruning (Appendix Z): python bnb_octree.py ; python bnb_scaling.py
