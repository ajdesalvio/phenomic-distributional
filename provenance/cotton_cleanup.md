# Cotton cleanup verification

The cotton workflow draws on the V6 full-season heritability models, their V5 formula and input functions, and the V7 functional comparisons. The numbered originals and manuscript outputs provide the reference for the checks below.

## Changes

- The current scripts share seven heritability functions and fifteen functional functions from the original code.
- The heritability command runs in sequence and saves a checkpoint for each date. It uses the V5 helpers directly.
- The functional command uses installed R libraries. It retains the original package checks, random number generator, defaults, and inference routines. Its provenance records checksums for the runner, model code, and input helper.
- The exporter uses saved estimates to rebuild the tables with the original heritability summary calculations.

## Evidence

- All 22 extracted function definitions match the originals exactly.
- The heritability smoke test checks the long and wide inputs for all three vegetation indices and fits nine models. The estimates of Vg, Vb, Ve, and heritability agree with the production results within 1e-10. The cohort checksums match exactly.
- All 30,300 saved heritability estimates satisfy the original formula. All six summary tables match the manuscript tables within 1e-12.
- The functional defaults match production for all three indices. The ExG2 preparation and inspection recover the original design with 47,807 plant-date observations, 479 plants, 101 dates, 16 genotypes, ten blocks, and a fixed design rank of 125.
- The functional contrasts, bands, and adjusted means match the saved outputs within 1e-12. Smoothing, propagated covariance, and standard errors calculated from the 300 saved point fits agree within 1e-10.
- The 4,999 bootstrap weight vectors per index match exactly. Replaying the first three draws per index reproduces the saved bootstrap-t statistics within 1e-10. Joint and per-curve critical values agree with the saved draw summaries.
- A synthetic repeated-measures fit gives identical coefficients, variance components, CR2 calculations, and diagnostics in the original and current code. Elapsed time was excluded from the comparison.

The checks used saved results and selected refits or replays. The full 30,300-model and 4,999-draw runs were outside the scope of this review.

## Interpretation retained

Cotton heritability is `Vg / (Vg + Ve)`, conditional on block. The reported estimates concern separate within-date models. The analysis excludes heritability confidence intervals and genotype-by-date variance components. Test identifies sections of the field. The planned irrigation treatment was never applied.

The functional bands are approximate and condition on the fitted covariance weights and smoothers, so they exclude uncertainty from estimating those quantities. The bootstrap accounts for dependence within a plant and assumes independence across plants after adjustment. The joint adjustment covers 15 contrasts against TM-1 over the quantile grid for each vegetation index separately. The additive model estimates genotype contrasts averaged across flights. Date-specific differences in senescence would require a model that allows genotype effects to vary by date.

The RData inputs contain prepared image quantiles. Validation starts there. See the [release guide](../data/DOWNLOADS.md) for the cotton inputs and saved models.
