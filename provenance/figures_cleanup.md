# Correlation, figure, and table cleanup

The current renderers use the final numerical calculations and plot settings. Input paths are configurable, and the workflow uses the saved inputs and results for CV, FPCA, mixed models, and bootstrap inference.

## Verification

- The correlation script reproduces the saved V5 calculations exactly, including all 192 environment and date matrices, pair counts, retained quantiles, screening diagnostics, ranges, and 120,000 plotted correlations. The V6 renderer includes the full quantile 1 row and column in each panel.
- Feature importance from the archived V15 outer coefficients matches Figure 7 and the supporting data exactly. The script uses the original normalization and averaging.
- The values plotted for Figure 3 contrasts and bands, Figure 4 components and season summaries, supporting flight-time statistics, and Figure 6 match the saved inputs within an absolute tolerance of 1e-12.
- Both Table S3 CSVs match the final files byte for byte. The calculations give each environment equal weight and compute percentage RMSE reductions before averaging across environments.
- All 30 rendered PNGs matched the final 600 dpi reference PNGs pixel for pixel. The renderer uses the original PDF graphics device for the Figure 4 layout and the PNG device for the prediction layouts to preserve text and sizing. The review also checked resolution, dimensions, and a visual contact sheet.

The renderer covers Figures 3, 4, 6–8 and S1–S25. Figure 5 has a separate [heatmap workflow](../scripts/figures/HAM.md). The editable [Figures 1 and 2](../figures/source/README.md) are PowerPoint artwork. Figure 1 uses the author’s corrected 17 April flight start date. Validation of these files did not extend to reconstructing the photographs or raw quantile artwork.

## Metadata and limitations

Table 1 describes model scope and Table S1 records genotype information. Both were transcribed from the accepted manuscript and supplement. The transcription skips deleted revisions and repeats labels from vertically merged cells. Source checksums are in [artwork and table provenance](authored_artifacts.json).

The quantile correlations are descriptive. The functional bands are approximate and conditional on fitted weights and smoothers. The maize heritability calculations use saved components, and the primary LOEO results use the historical fits. See the [release guide](../data/DOWNLOADS.md) for the inputs to each renderer.
