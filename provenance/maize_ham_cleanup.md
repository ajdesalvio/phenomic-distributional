# Maize heritability and heatmap validation

The [heritability workflow](../scripts/maize/HERITABILITY.md) recalculates heritability from saved variance components. The 11 manuscript environments contain 629,983 fits across 33 index labels. After grouping equivalent indices, the summaries cover 592,119 available fits across 31 index representations and 192 environment-specific dates.

Figure inputs, Table S2, and the remaining manuscript summaries match the reference results within an absolute tolerance of 1e-12. The code rejects duplicate components, flags invalid variances or denominators, and keeps missing quantiles as missing. Heatmaps require complete curves. Other summaries use the available valid quantiles. The checks cover calculations from saved components. Assessing model convergence requires the underlying fits.

The [heatmap workflow](../scripts/figures/HAM.md) reproduced all eleven maps from the source TIFFs and the included curve tables. Vegetation index arrays, ranks, masks, curves, mapped heritability arrays, and rendered RGB and heatmap pixels agree exactly. The 1,100-row curve table matches byte for byte. One reported mean differs by 5.55e-17 under the tested Python libraries.

Each heatmap uses ranks calculated within an image to display values from the matching quantile heritability curve. See the [input provenance](ham_input_excerpts.json) for source checksums and row counts and the [release guide](../data/DOWNLOADS.md) for the TIFFs.
