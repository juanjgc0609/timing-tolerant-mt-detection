# Deviations of the core environment from the frozen runs

Reference of the frozen runs: `requirements-frozen.txt` (Python 3.11) and the full list of the rebuilt Task 2M
environment, `results/task2M/env_check/env_venv311.txt` (Python 3.11.17, Linux).

| Item | Frozen runs | Core environment | Why |
|---|---|---|---|
| Python | 3.11.17 | 3.12.15 | decision of the Phase 1 brief (numpy 1.26.4 has 3.12 wheels) |
| Direct pins of `requirements-frozen.txt` | numpy 1.26.4, scipy 1.17.1, obspy 1.5.1, pandas 3.0.6, scikit-learn 1.8.0, pyocto 0.2.0, matplotlib 3.11.2, pyproj 3.7.2 | identical | |
| Transitive pins | `env_venv311.txt` | identical to `env_venv311.txt` for every package it lists | |
| Added | — | numba 0.61.2 (lazy import in `task2c_methods.py`), pytest 8.4.2 with iniconfig, pluggy, Pygments; colorama on Windows only; tzdata 2026.5 (pandas on Windows); llvmlite 0.44.0 (numba) | needed by imports or by the tests |
| OpenBLAS | the build bundled in the numpy 1.26.4 wheel, kernel Prescott (E31) | same wheel family; `ttmt.runtime.pin_blas()` sets `OPENBLAS_CORETYPE=Prescott` and one thread before numpy loads | on Windows the wheel's OpenBLAS and the C runtime's libm differ from the Linux build, so bit-for-bit equality is only claimed on Linux (CI); this machine uses a declared tolerance |
| PhaseNet | tensorflow-cpu 2.15.1 in a separate target directory | not installed | tensorflow 2.15 has no Python 3.12 wheel; the picks of R1 are read from the stored `picks_all.json` (see the R1 regression) |
| GaMMA | commit 80394dd4 cloned into `code/gamma_src/` | not installed | only the pick stage of R1 and R0 needs it |

No frozen pin had to be replaced by a nearby version.
