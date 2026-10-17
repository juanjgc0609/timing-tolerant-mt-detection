#!/bin/bash
# Task R1 freeze: sha256 of every Python module that the run can import (top level of the work folder, r0/, r1/, hyb/,
# explore_contrib/), the run scripts, the data files, the GaMMA commit, the PyOcto version and the PhaseNet model.
cd /home/claude/work || exit 1
. r1/env.sh
F=r1/FREEZE_R1.sha256
{ for d in . r0 r1 hyb explore_contrib; do ls $d/*.py; done | sort | xargs sha256sum
  sha256sum r1/r1_queue.sh r1/env.sh
  find /home/claude/realdata/MatchLocate2/example/Trace /home/claude/realdata/MatchLocate2/example/catalog.dat \
       /home/claude/realdata/MatchLocate2/example/MultipleTemplate/DetectedFinal.dat -type f | sort | xargs sha256sum
  find /home/claude/work/r0_verif/phasenet/PhaseNet/model/190703-214543 -type f | sort | xargs sha256sum
  ls -d /home/claude/work/r0_verif/phasenet/tf215/tensorflow*-*.dist-info 2>/dev/null | xargs -n1 basename
  echo "gamma_src commit $(cd gamma_src && git rev-parse HEAD 2>/dev/null || cat COMMIT 2>/dev/null)"
  /home/claude/venv311/bin/python -c "from importlib.metadata import version as v; print(*[f'{p} {v(p)}' for p in ('pyocto', 'obspy', 'numpy', 'scipy', 'pandas', 'scikit-learn')])"
} > $F
sha256sum $F
