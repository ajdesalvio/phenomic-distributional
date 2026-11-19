# Figure and table validation

The [figure and table workflow](../scripts/figures/README.md) uses the saved analysis results to produce Figures 3, 4, and 6–8, Figures S1–S25, and Table S3. Figure 5 uses the [heatmap workflow](../scripts/figures/HAM.md). Figures 1 and 2 have [editable PowerPoint files](../figures/source/README.md).

- The correlation script reproduces all 192 environment and date matrices, including pair counts, screened quantiles, and all 120,000 plotted correlations. The renderer includes the p = 1.00 row and column.
- Feature summaries rebuilt from the saved outer-fold coefficients match Figure 7 and the supporting figures.
- The values plotted in Figures 3, 4, and 6 and the supporting flight-time plots match the saved inputs within an absolute tolerance of 1e-12.
- All 30 rendered PNGs match the 600 dpi reference images pixel for pixel in the tested environment.
- The two CSVs for Table S3 match the reference files byte for byte. The calculation gives each environment equal weight and computes percentage RMSE reductions within environments before averaging.

Table 1 summarizes model scope and Table S1 records genotype information. The [artwork and table provenance](authored_artifacts.json) identifies their sources. See the [release guide](../data/DOWNLOADS.md) for the inputs to the renderers.
