# Across-environment prediction

Install `data.table`, `jsonlite`, `digest`, `caret`, and `glmnet`. See the [release guide](../../data/DOWNLOADS.md) for the large inputs. If you store the files outside their default locations, edit the paths in `config/loeo.json`. Relative paths start from that file’s directory.

## Reproduce the manuscript results

Run this from the repository root with a new output directory.

```sh
Rscript scripts/maize/export_loeo_reference.R config/loeo.json output/loeo_reference
```

The command checks the saved model archive, exports predictions and coefficients, and recalculates the comparison of Median and all quantiles across the 11 held-out environments.

## Optional refitting

Use the original projection table and yield file.

```sh
python scripts/maize/partition_projections.py data/2020_2021_G2F_FPC_Scores_VI_BLUEs_Projections_AllQuantiles_ExGR.csv output/projections
Rscript scripts/maize/run_loeo.R config/loeo.json output/projections output/loeo_refit
```

The runner fits Median and all quantiles with within-environment scaling, one model at a time. Append an environment such as `TXH1.2020` to run that environment alone. Repeat the command to resume from matching checkpoints. A complete run has 22 models. Partial summaries show how many held-out environments have finished. Save refits separately from the manuscript results.

The workflow retains 95% of the variance for each quantile, uses seed 123, groups inner folds by environment, and runs caret’s original `tuneLength=25` search. Inner tuning reuses the FPCA from the outer training set. Scaling uses predictors from the unlabeled held-out cohort. The output records caret’s historical interpolation or clamping of lambda values.

## Optional FPCA preparation

List all five archived projection parts in their original order to assemble them.

```sh
Rscript scripts/maize/assemble_projections.R output/projections.csv PART1.csv PART2.csv PART3.csv PART4.csv PART5.csv
```

To recompute FPCA, also install `fdapace`, `dplyr`, and `tidyr`.

```sh
Rscript scripts/maize/project_fpca.R config/loeo.json output/fresh_projections.csv
```

The default includes all 100 quantiles. Append `TXH1.2020 0.5` to run a single environment and quantile. FPCA uses curves from the outer training set. Held-out curves need at least two observations within the training time range. The script saves diagnostic plots and R package versions alongside the CSV.

Fresh FPCA follows the historical algorithm. For a new analysis with the resulting projection table, record its checksum in the configuration. Keep the published projection table and checksum for manuscript reproduction. See the [validation notes](../../provenance/loeo_cleanup.md).
