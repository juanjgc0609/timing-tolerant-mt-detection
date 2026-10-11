# Code

The scripts keep the flat layout of the original working directory (`/home/claude/work`), because they import one
another by module name and the exploration and real-data scripts add the parent folder to `sys.path`. Run them from
this folder, or recreate that path.

## Core modules

| Module | Role |
|:--|:--|
| `task0_common.py` | constants (sampling, window length, velocities, origins), networks, pick formats |
| `gen_realistic.py` | synthetic generator: gradient medium, curved rays, random double couples, Ricker pulses with phase rotation, coda, site factors, statics, clutter |
| `task1x_mtglrt.py`, `certified_fw.py` | moment-tensor radiation operators, subspace bases, certified Frank-Wolfe search |
| `task2c_methods.py`, `task2d_methods.py`, `task2e_methods.py` | block dictionaries, BOMP-PK, the convex and off-grid (TS-SFW) solvers |
| `station_corrections.py` | catalogue station corrections (+SC) from matched-filter peak delays |
| `real_picker.py` | PhaseNet (original model) and STA/LTA pickers |
| `task1r_baselines.py`, `qm_wrapper.py` | GaMMA, PyOcto and QuakeMigrate wrappers |
| `froc_exact.py` | exact FROC at fixed realized false-event rates |

## Tasks

Each confirmatory task has a runner, an analysis script and, from Task 2D on, a checksum file of its frozen code:
`task2a_confirm.py` + `analyze_2a.py` (2A), `task2b_ablation.py` (2B), `task2c_run.py` (2C), `task2d_run.py` (2D),
`task2e_run.py` (2E), `task2f_run.py` (2F), `task2i_run.py` (2I), `task2l_run.py` + `task2l_control.py` (2L), each with
its `analyze_*.py`; the Task 1 programme scripts are `task1*.py` and `analyze_1*.py`. The prepared Task 2M is
`t2m_data.py`, `t2m_stations.py`, `t2m_lib.py`, `t2m_run.py`, `t2m_validate.py`, `t2m_analyze.py`, `t2m_power.py` and
`t2m_*.sh` (its configurations are in `results/task2M_prepared/`, expected under `t2m_out/`).

`explore_contrib/` holds the explorations C1-C16 on the calibration seed (`cN_lib.py`, `cN_run.py`, `cN_analyze.py`),
with the shared libraries `xc.py`, `c8_lib.py`, `c9_lib.py`, `c10_lib.py` (the MT-FW loop with timing modes),
`c11_lib.py` (pick-free JHD station terms S2) and the freeze files. `r0/` holds the real-data task R0.

## Notes

Comments, docstrings and printed labels inside the code are partly in Spanish; the documents and the analysis outputs
in `results/` are in English. GaMMA is not included: clone `AI4EPS/GaMMA` at commit
`80394dd4a4c29450da599f507916ecf1db6203dd` into `gamma_src/`. Scripts in this folder that are not referenced by any
task or appendix are earlier development scripts kept for completeness.
