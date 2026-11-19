# Cotton workflow validation

The [cotton workflow](../scripts/cotton/README.md) covers heritability for PSRI, CVI, and ExG2 over the full season, along with functional comparisons across flights. Use the saved results to reproduce the manuscript tables and bands.

- The saved heritability summaries, functional contrasts, confidence bands, and adjusted means match the manuscript results within 1e-12.
- The heritability smoke test fits nine models across the three vegetation indices. Variance components and heritability match the reference fits within 1e-10.
- Recomputing the smoothing and propagated uncertainty from all 300 saved point fits agrees within 1e-10. The bootstrap weights match for all 4,999 draws per index. Replaying three draws per index reproduces the saved bootstrap-t statistics within 1e-10.

## Interpretation

Cotton heritability is `Vg / (Vg + Ve)`, conditional on block. The estimates come from separate models for each date, vegetation index, and quantile. Test identifies sections of the field. The planned irrigation treatment was never applied.

The functional bands are approximate and condition on the fitted covariance weights and smoothers. The bootstrap accounts for dependence within a plant and assumes independence across plants after adjustment. Joint adjustment covers 15 contrasts against TM-1 across the quantile grid separately for each vegetation index. The additive model estimates genotype contrasts across flights. Date-specific differences in senescence would require genotype effects that vary by date.

Validation uses prepared image quantiles, saved fits, and selected refits and bootstrap replays. See the [release guide](../data/DOWNLOADS.md) for the inputs and saved models.
