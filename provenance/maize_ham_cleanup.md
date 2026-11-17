# Maize heritability and HAM cleanup

The workflow recalculates heritability from saved variance components and uses the original heatmap functions. Validation used the saved fits and source TIFFs.

## Maize

The heritability function matches V1 exactly. The runner uses V2’s strict parsing of decimals and scientific notation and the final manuscript calculations. It reads 5,748,865 component rows. The 11 manuscript environments contain 629,983 fits across 33 index labels. After grouping equivalent indices, the summaries cover 592,119 available fits across 31 representations and 192 environment-specific dates.

The saved-fit evidence, fit coverage, figure components, complete map curves, Table S2, quantile summaries, affine comparison, and flight-time results match the manuscript inputs within an absolute tolerance of 1e-12. Differences around 1e-16 in the affine summaries trace to an intermediate CSV export and reimport. The checks use absolute tolerance because some of the values compared are near zero.

The code rejects duplicate components and flags negative variances or zero denominators. Missing quantiles stay missing. Heatmaps require a complete curve, while summaries use any valid individual quantiles. The comparison starts from saved components, with historical convergence and raw-image processing outside its scope.

## HAM

Seven helper functions and the two focal vegetation index formulas match V17 exactly. The runner uses explicit input paths and uses the original calculations for ranks within each image. Figure 5 uses the original 600 dpi layout.

We regenerated all eleven maps from their original TIFFs. The vegetation index arrays, ranks, masks, heritability curves, and mapped arrays match exactly, including invalid pixels. RGB, heatmap, and color bar PNGs have identical decoded pixels. The 1,100-row heritability curve table matches byte for byte. One summary mean differs by 5.55e-17 under the tested Python libraries. Its underlying arrays match exactly.

The bundled excerpts reproduce the same maps as the full source tables. They contain 500 cotton and 600 maize curve rows, five cotton and six maize date records, and six maize image membership records. The [input provenance](ham_input_excerpts.json) lists the source checksums and row counts.

The workflow uses right-continuous ranks within each image and excludes nodata, black, and nonfinite pixels. It keeps the original interpolation, band formulas and scales, rotations, and colors. The resulting maps display heritability from the quantile curve at each pixel’s rank. The estimates and their interpretation stay at the quantile level.

The detailed checks and scripts are in the local review directory. See the [release guide](../data/DOWNLOADS.md) for the saved components and eleven TIFFs.
