# ExGR optimizer sensitivity

On 11 November 2026, we kept the archived inputs as the reference for the manuscript. The refits below were used to investigate differences in BLUEs.

We refit all 35,900 ExGR models and matched all 21,326,300 genotype and model keys. For 2,700 models, at least one BLUE differed by more than 10^-7. Seventeen exceeded 10^-4, and the largest absolute difference was 0.001395. The fits produced 815 messages about convergence, gradients, or Hessians despite zero optimizer exit codes.

For those 17 models and one control, the archived variance ratios reproduced the archived BLUEs on the current raw data within 4.50 × 10^-14. The rebuilt ratios reproduced the rebuilt BLUEs within 6.55 × 10^-14. This points to covariance estimation as the source of the largest discrepancies tested.

We used two tighter optimizers on the six largest cases and the control. The tightened estimates agreed within 1.28 × 10^-7. The checks gave no consistent reason to prefer the archived or rebuilt solution. The Pearson correlation between archived and rebuilt genotype BLUEs was at least 0.9997305 in the selected cases. Downstream predictions after retraining were outside the scope of this comparison.

Nine of the 17 largest discrepancies were outside the final 11 prediction environments. The remaining eight occurred inside those environments and were used in some saved models. The archived inputs remain the reference for the reported results. See the [reference tables](../results/reference/) and [script lineage](source_map.csv) for the numerical evidence and checksums. The audit starts from the existing image-level quantiles.
