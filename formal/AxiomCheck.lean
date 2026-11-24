import TTMT

/-! Prints the axioms used by each formal result. Expected: only `propext`,
`Classical.choice`, `Quot.sound` (never `sorryAx`). Run: `lake env lean AxiomCheck.lean`. -/

#print axioms TTMT.pAmplitude_cross_eq_zero
#print axioms TTMT.kappa_P_eq_zero
#print axioms TTMT.gain_null_eq
#print axioms TTMT.gain_null_pivotal
#print axioms TTMT.reading_contraction
#print axioms TTMT.lipschitz_comp_one
#print axioms TTMT.lipschitz_sup
#print axioms TTMT.coincident_pair_mem_oneEventClass
#print axioms TTMT.Controls.generic_cross_ne_zero
#print axioms TTMT.Controls.gain_depends_on_signal_outside_class
#print axioms TTMT.Controls.overlapping_rows_not_contraction
#print axioms TTMT.Controls.gain_needs_idempotent
#print axioms TTMT.Controls.gain_needs_symmetric
#print axioms TTMT.Controls.gain_positive_control
#print axioms TTMT.Controls.large_row_not_contraction
#print axioms TTMT.Controls.nonlinear_pair_not_in_class
