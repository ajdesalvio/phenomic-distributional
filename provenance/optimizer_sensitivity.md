# ExGR optimizer sensitivity

Use the saved ExGR BLUEs to reproduce the manuscript. Fresh fits can return slightly different estimates because the fitted covariance components can vary with the optimizer.

Refitting 35,900 models recovered all 21,326,300 genotype and model keys. Absolute differences in BLUEs exceeded 10^-7 in 2,700 models and 10^-4 in 17 models. The largest difference was 0.001395. The fits produced 815 messages about convergence, gradients, or Hessians despite zero optimizer exit codes.

For the 17 largest discrepancies and one control, applying the saved and refitted variance ratios to the same raw data reproduced the corresponding BLUEs within 6.55 × 10^-14. This points to covariance estimation as the source of the discrepancies tested. Eight of the 17 cases occurred in environments used for prediction. The effect of refitting on downstream predictions remains unverified. The [reference tables](../results/reference/) give the numerical evidence.
