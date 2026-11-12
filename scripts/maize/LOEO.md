# Across-environment prediction

Install `data.table`, `jsonlite`, `digest`, `caret`, and `glmnet`. Set paths in `config/loeo.json`; relative paths start from the configuration directory. Large input files are not bundled.

## Reproduce the manuscript results

From the package root:

```sh
Rscript scripts/maize/export_loeo_reference.R config/loeo.json output/loeo_reference
```

This verifies the frozen model archive, exports its primary predictions and coefficients, and recalculates the metrics for Median versus all quantiles across 11 held-out environments. It does not fit models. Use a new output directory.

## Optional refitting

Use the frozen all-quantile projection table and the same yield input:

```sh
python scripts/maize/partition_projections.py data/2020_2021_G2F_FPC_Scores_VI_BLUEs_Projections_AllQuantiles_ExGR.csv output/projections
Rscript scripts/maize/run_loeo.R config/loeo.json output/projections output/loeo_refit
```

The runner fits only Median/all with within-environment scaling, sequentially. Append an environment such as `TXH1.2020` to fit just that environment. Repeat the command to resume matching checkpoints. Completion requires 22 models; partial summaries show their completed held-out counts. Keep refits separate from manuscript reference outputs.

The pipeline retains 95% FVE per quantile, seed 123, environment-grouped inner folds, and caret's original `tuneLength=25` search. Saved outer-training FPCA is reused during inner tuning. Scaling uses the unlabeled held-out predictor cohort. This is not fully nested FPCA or single-plant deployment. Caret's historical lambda interpolation/clamping is retained and reported.

## Optional FPCA preparation

To assemble the five archived projection parts, in their original order:

```sh
Rscript scripts/maize/assemble_projections.R output/projections.csv PART1.csv PART2.csv PART3.csv PART4.csv PART5.csv
```

For fresh FPCA, also install `fdapace`, `dplyr`, and `tidyr`:

```sh
Rscript scripts/maize/project_fpca.R config/loeo.json output/fresh_projections.csv
```

The default covers all 100 quantiles. Append `TXH1.2020 0.5` for one held-out environment and quantile. FPCA uses outer-training curves only; held-out curves are restricted to the training time range and require two observations. Diagnostic plots and R versions are saved alongside the CSV.

Fresh FPCA is an optional reconstruction of the historical algorithm, not a certified recreation of every archived basis. Its output is deliberately not accepted by the frozen projection hash unless the configuration is explicitly changed for a new analysis. See [verification and limitations](../../provenance/loeo_cleanup.md).
