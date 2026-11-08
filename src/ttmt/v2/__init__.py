"""Model v2 components (specification ``especificacion-modelo-v2.md``, sections M0, M1, M4, M5; plan T4.13, T4.15).

Exploratory: no registered run uses this subpackage. Every function takes its conventions as arguments and leaves the
registered detector (``ttmt.detector``) unchanged.

  m0     noise covariance 3x3 per station, whitening of traces and operators, orientation audit
  m1     onset evidence q (picker probability on the detector grid), isotonic calibration, onset reader,
         onset-to-peak lag, sub-sample amplitude reading
  m4     joint reading of one hypothesis (all blocks chosen together by consistency with the moment-tensor subspace)
  loop   the v2 greedy loop: onset (or joint) reading, ACAP ranking (M5), deflation at the sub-sample read time
"""
