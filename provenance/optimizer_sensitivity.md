# ExGR optimizer sensitivity

**Decision, 11 November 2026:** retain the existing manuscript inputs and results. Keep fresh model regeneration separate from the frozen reference. No manuscript reanalysis was performed for this audit.

All 35,900 ExGR models were refitted. All 21,326,300 genotype/model keys matched, but 2,700 models exceeded an absolute BLUE difference of 10^-7. Seventeen exceeded 10^-4; the maximum was 0.001395. There were 815 convergence/gradient/Hessian messages despite zero optimizer exit codes.

For those 17 cases plus one control, using archived variance ratios on current raw data reproduced archived BLUEs within 4.50 × 10^-14. Using rebuilt ratios reproduced rebuilt BLUEs within 6.55 × 10^-14. These checks identify covariance estimation, rather than changed data or formula, as the source of the largest discrepancies tested.

Two tighter optimizers were run on the six largest cases and the control. Neither the archived nor rebuilt solution was uniformly preferable. The tightened estimates agreed within 1.28 × 10^-7. Across the selected cases, archived/rebuilt genotype BLUE Pearson correlations were at least 0.9997305; this does not establish identical downstream predictions after retraining.

Nine of the 17 largest discrepancies were outside the final 11 prediction environments. Eight were inside and were used in some saved models. Accordingly, the audit does not claim zero downstream impact or exact first-stage refit equivalence. Original inputs remain the reference for reported results; the comparison failure remains documented.

Small numerical evidence tables are in `results/reference/`. Script lineage and hashes are in `source_map.csv`. Original pixel extraction, public data hosting, and other repository workflows remain separate audit tasks.
