# Remaining figures and tables

Install the R packages in `environment/R_versions.csv` and use Arial to match the prediction figures. The renderer checks the saved input files against the checksums in `config/figures.json`. If you move the inputs, update the paths relative to that file.

From the package root:

```sh
Rscript scripts/figures/render_results.R config/figures.json output/figures
python scripts/tables/prediction_table.py config/figures.json output/table_s3
```

The renderer writes 30 plots as 600 dpi PNG and JPEG files and as PDFs. Append `3`, `4`, `6`, `7`, or `8` to render one figure group. Figure 4’s group includes S1–S14, and Figure 7’s includes S15–S25. Choose a new output directory.

The former S16–S26 figures are now S15–S25, and the former prediction Table S6 is now S3. The plots retain their original data, colors, layouts, uncertainty bands, and dimensions. Figure 6 is 7.25 inches square.

## Optional numerical reconstruction

Use these commands to recalculate correlations and feature summaries from saved inputs.

```sh
Rscript scripts/maize/quantile_correlations.R config/nested_cv.json output/correlations
Rscript scripts/maize/prepare_feature_figures.R V15_RUN_ROOT output/feature_data
```

The correlation script screens genotype BLUEs with the corresponding raw support within each environment and date. It calculates descriptive correlations, averages them on the Fisher z scale, and transforms the averages back. In the combined matrix, each valid environment and date contributes equally.

For feature summaries, the script scales absolute standardized coefficients by the largest one within each model, then averages across the 25 outer fits per scenario. Excluded predictors and zero coefficients contribute zero. If you plot regenerated data, check their values before updating the figure configuration paths and checksums.

## Other outputs

| Item | Source or command |
|---|---|
| Figure 5 | [HAM workflow](HAM.md) |
| Figures 1–2 | [Editable authored artwork](../../figures/source/README.md) |
| Table 1 and S1 | Curated CSVs in `metadata/` |
| Table 2 | [Cotton saved-result export](../cotton/README.md) |
| Table S2 | [Maize heritability](../maize/HERITABILITY.md) |
| Table S3A/B | `prediction_table.py`, above |

Table S3 gives each environment equal weight. The script calculates the percentage RMSE reduction for each environment before averaging those reductions. It leaves the overall N cell blank. Use the CSVs as the source for the table values.

See [verification and remaining limits](../../provenance/figures_cleanup.md).
