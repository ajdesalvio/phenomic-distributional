# Remaining figures and tables

Set paths in `config/figures.json`; relative paths start from the configuration file. Inputs are hash-pinned saved results. Install the R packages listed in `environment/R_versions.csv`. The original prediction figures use Arial; retain that font for matching exports.

From the package root:

```sh
Rscript scripts/figures/render_results.R config/figures.json output/figures
python scripts/tables/prediction_table.py config/figures.json output/table_s3
```

The renderer exports Figures 3, 4, 6, 7, 8 and S1–S25: 30 plots, each as 600-dpi PNG/JPEG plus PDF. Append `3`, `4`, `6`, `7`, or `8` to render one group. Group 4 includes S1–S14; group 7 includes S15–S25. Use a new output directory. Existing results are never overwritten.

The current numbering replaces historical S16–S26 filenames with S15–S25 and the old prediction-table S6 filenames with S3. Plot data, colors, layouts, uncertainty bands, and physical dimensions are retained. Figure 6 remains 7.25 inches square.

## Optional numerical reconstruction

These commands use saved inputs, not new model fits:

```sh
Rscript scripts/maize/quantile_correlations.R config/nested_cv.json output/correlations
Rscript scripts/maize/prepare_feature_figures.R V15_RUN_ROOT output/feature_data
```

The correlation calculation uses genotype BLUEs and raw support within each environment/date. It retains the V15 screening rules, averages on the Fisher-z scale, and back-transforms. Combined correlations weight valid environment-date matrices equally, not environments equally. Correlations are descriptive and are not a cross-validation procedure.

Feature summaries average model-relative absolute standardized coefficients over the 25 outer fits per scenario, with zero/excluded coefficients contributing zero. Update the figure configuration paths and hashes explicitly when using regenerated data.

## Other outputs

| Item | Source or command |
|---|---|
| Figure 5 | [HAM workflow](HAM.md) |
| Figures 1–2 | [Editable authored artwork](../../figures/source/README.md) |
| Table 1 and S1 | Curated CSVs in `metadata/` |
| Table 2 | [Cotton saved-result export](../cotton/README.md) |
| Table S2 | [Maize heritability](../maize/HERITABILITY.md) |
| Table S3A/B | `prediction_table.py`, above |

S3 reports equal-weight environment means; percentage RMSE reductions are calculated per environment before averaging. The overall N cell is blank. These CSVs supply values, not Word formatting or tracked changes.

See [verification and remaining limits](../../provenance/figures_cleanup.md).
