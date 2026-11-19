# Within-environment nested CV validation

The [nested CV workflow](../scripts/maize/README.md) compares the median quantile with all 100 quantiles. It runs each environment separately in the C locale to reproduce the fold assignments.

- The outer and inner fold assignments match the manuscript runs for all 11 environments.
- For TXH1.2020, one complete outer split per predictor set was refit with the full tuning grid. The results match the saved fits within 1e-12.
- Recalculating the combined summaries from saved predictions reproduces the manuscript tables exactly.
- Synthetic checks cover training-only preprocessing and tuning, constant-feature removal, the intercept-only fallback, checkpoint resume, and rejection of changed inputs.

These checks cover representative refits and all saved summaries. A complete refit across all environments remains unverified. Run the synthetic check from the repository root.

```sh
Rscript validation/check_nested_cv.R output/nested_synthetic_check
```

See the [release guide](../data/DOWNLOADS.md) for the inputs and saved predictions and the [software guide](../environment/README.md) for the tested setup.
