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

## Picks environment (`picks.lock.txt`, `picks.yml`)

Scenario: the pick stage of R0 and R1 (PhaseNet, GaMMA, PyOcto) and, later, the comparators of the atlas.

| Item | Frozen runs | Picks environment | Why |
|---|---|---|---|
| Python | 3.11.17 | 3.11.17 | TensorFlow 2.15 has no wheel for Python 3.12 |
| TensorFlow | tensorflow-cpu 2.15.1, tensorflow-estimator 2.15.0, tensorflow-io-gcs-filesystem 0.37.1 (FREEZE_R1) | the same, except tensorflow-io-gcs-filesystem 0.31.0 on Windows (no newer Windows wheel), and tensorflow-intel 2.15.1, which tensorflow-cpu requires on Windows | platform; neither package takes part in the PhaseNet computation |
| PhaseNet | AI4EPS/PhaseNet 62005c6, model 190703-214543 | the same commit; the six model files match the sha256 of FREEZE_R1 | |
| GaMMA | 80394dd, cloned into code/gamma_src/ | the same commit, cloned outside the checkout (GAMMA_SRC on sys.path); the six files of R0_FREEZE match | keeps code/ untouched |
| numba, llvmlite, tqdm (GaMMA) | r0_verif/pydeps, not frozen; recorded by the R0 verification as numba 0.60.0, llvmlite 0.43.0, tqdm 4.66.5 | the same versions | |
| pyrocko, pyyaml (PyOcto 1-D model) | 2026.6.2, 6.0.3 (data/README.md) | the same | |
| Core pins | requirements-frozen.txt and env_venv311.txt | identical | |

## Numerical tolerance off Linux

The frozen runs were made on Linux. On Windows the C runtime's libm and the OpenBLAS bundled in the numpy wheel round
differently. The difference is already visible upstream of the detector: the R1 station scale after filtfilt and
resample_poly differs by 4.6e-12 (relative). The R1 P1 scores reach 3.6e-9 absolute, 4.8e-13 relative, on scores up
to 19 335. The legacy code run on the same Windows machine gives values identical to the package's on both
segments, bit for bit, so the difference comes from the platform, not from the port. Rule (`ttmt.testing`): discrete
outputs exact on every platform; continuous detector scores and station terms bit for bit on Linux and within
relative 1e-12 elsewhere.

## Comparator environments of the atlas (`gpu.lock.txt`, `gpu_weights.sha256`, `eqc.conda-explicit.txt`)

`ttmt-gpu` is a conda clone of the local `mt-gpu-test` env, which stays unmodified: Python 3.12.14, torch
2.13.0+cu130, CUDA on the RTX 4060 Laptop GPU. Installed into it:
- SeisBench 0.12.6, with the PhaseNet and EQTransformer weights `original`, `stead`, `instance`, `scedc` and `ethz`
  (`SEISBENCH_CACHE_ROOT=C:\Users\juanj\ttmt\seisbench_cache`; sha256 in `gpu_weights.sha256`);
- PyOcto 0.2.0 with pyrocko 2026.6.2;
- GaMMA at 80394dd (editable, no dependencies), with numba 0.68.0 and scikit-learn 1.8.0, the frozen version (GaMMA's
  own metadata asks for 1.7);
- PhaseNO v1.0.0 (sun-hongyu/PhaseNO, commit 1c6d3f0, archived as Zenodo 10.5281/zenodo.10224301). The repository
  given in the brief, AI4EPS/PhaseNO, does not exist. PhaseNO's code is a notebook written for torch 1.12, PyG 2.2 and
  pytorch-lightning 1.9; here it runs with torch_geometric 2.8.1 and pytorch-lightning 2.6.6. Its example files have
  ':' in their names, so the tree is exported from git blobs with ':' replaced by '-';
- EQNet at 7dd24ed (editable). The current repository serves PhaseNet, PhaseNet+ and PhaseNet-DAS. The 2022 EQNet
  model (shift-and-stack of learned features) has no released weights. An `EQNet2-v0` weight release (AI4EPS/models,
  11 Oct 2026) has no public inference code. So EQNet cannot run in the atlas core set until weights or inference
  code for it are public.

numpy is 2.5.2 here, against 1.26.4 in the core env. The comparators are not part of the frozen pipeline.

EQcorrscan 0.5.2 has no Windows wheel on PyPI, and its source build fails. conda-forge has a Windows build, but it
would replace numpy, scipy and MKL under the pip-installed torch stack of `ttmt-gpu`. It therefore lives in a separate
env, `ttmt-eqc` (conda-forge only, Python 3.12, numpy 2.5.4, obspy 1.5.1; `eqc.conda-explicit.txt`).
