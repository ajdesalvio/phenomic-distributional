# Across-environment prediction validation

The manuscript compares the median quantile with all 100 quantiles using within-environment scaling. The saved archive contains 22 models, one for each predictor set in each of the 11 held-out environments. Use the [LOEO workflow](../scripts/maize/LOEO.md) to export predictions and recalculate the reported metrics.

- Predictions from the saved models match the archived predictions within 1e-10. Recalculated metrics and coefficient counts match the manuscript tables within 1e-12.
- For TXH1.2020, predictor matrices rebuilt from the saved projections match the archived scaled matrices within 1e-12 for both predictor sets.
- Assembling the five projection parts reproduces all 7,135,700 rows of the combined table. Synthetic data checks also cover FPCA and model tuning.

Use the saved models and projections for manuscript reproduction. Exact agreement from a complete rerun of FPCA and LOEO fitting remains unverified. BLUE estimation and vegetation index selection took place before cross-validation. Inner tuning reuses FPCA estimated in the outer training set, and within-environment scaling uses predictors from the unlabeled held-out cohort. Caret can interpolate or clamp a requested penalty to its fitted path.

The [release guide](../data/DOWNLOADS.md) lists the inputs. The [software guide](../environment/README.md) lists the tested environment.
