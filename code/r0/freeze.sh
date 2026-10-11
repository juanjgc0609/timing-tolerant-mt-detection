#!/bin/bash
# sha256 of every file the frozen R0 run and analysis execute or read
cd /home/claude/work
sha256sum r0/r0_lib.py r0/r0_run.py r0/r0_picks.py r0/r0_analyze.py r0/queue_eval.sh r0/R0_DEV_LOG.txt \
  explore_contrib/c10_lib.py explore_contrib/c11_lib.py explore_contrib/c8_lib.py explore_contrib/c9_lib.py \
  explore_contrib/c6_lib.py explore_contrib/xc.py station_corrections.py task1x_mtglrt.py certified_fw.py \
  task2c_methods.py task2d_methods.py task2e_methods.py task0_common.py task1o_consistent.py real_picker.py \
  task1r_baselines.py task2a_confirm.py gen_realistic.py froc_exact.py \
  $(find gamma_src/gamma -name '*.py' | sort) \
  /home/claude/realdata/skillsbench_eqassoc/wave.mseed /home/claude/realdata/skillsbench_eqassoc/stations.csv \
  /home/claude/realdata/skillsbench_eqassoc/catalog.csv
python3 -c "import pyocto, numpy, scipy, obspy, sklearn; print('pyocto', pyocto.__version__ if hasattr(pyocto,'__version__') else '0.2.0', 'numpy', numpy.__version__, 'scipy', scipy.__version__, 'obspy', obspy.__version__, 'sklearn', sklearn.__version__)"
