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
#print axioms TTMT.radial_lipschitz
#print axioms TTMT.goodRadius_id
#print axioms TTMT.goodRadius_cap
#print axioms TTMT.goodRadius_log
#print axioms TTMT.readingMap_lipschitz
#print axioms TTMT.norm_adjoint_le
#print axioms TTMT.norm_adjoint_sub_le
#print axioms TTMT.norm_sum_smul_le
#print axioms TTMT.blockRadial_lipschitz
#print axioms TTMT.readBlocks_lipschitz
#print axioms TTMT.tailFunction_lipschitz
#print axioms TTMT.read_le_tailFunction
#print axioms TTMT.score_sub_le_blocks
#print axioms TTMT.influence_le
#print axioms TTMT.cap_contamination
#print axioms TTMT.cap_score_le
#print axioms TTMT.coincident_pair_margin_zero
#print axioms TTMT.integral_quadForm
#print axioms TTMT.integral_gain_null
#print axioms TTMT.trace_mul_sum_vecMulVec
#print axioms TTMT.trace_eq_lambda_mul_card
#print axioms TTMT.Controls2.radial_needs_le_self
#print axioms TTMT.Controls2.radial_needs_lipschitz
#print axioms TTMT.Controls2.adjoint_needs_contraction
#print axioms TTMT.Controls2.influence_needs_one_block
#print axioms TTMT.Controls2.cap_score_needs_contraction
#print axioms TTMT.Controls2.contamination_needs_equal_outside
#print axioms TTMT.Controls2.mean_needs_signal_in_class
