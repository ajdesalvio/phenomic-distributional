# Maize heritability and HAM cleanup

This package extracts the original saved-component recalculation function and HAM helpers, with explicit configuration and short entry points. Originals and manuscript outputs remain untouched. No models were fitted.

## Maize

The recalculation function has identical arguments and body to V1. The runner retains V2's strict decimal/scientific parsing of saved variance text and the final evidence/summary calculations. It reads 5,748,865 component rows; the selected 11-environment evidence contains 629,983 fits across 33 labels. Canonical summaries contain 592,119 available fits across 31 index representations and 192 environment-specific dates.

The saved-fit evidence, coverage, figure components, complete map curves, Table S2, quantile/VI/date summaries, affine comparison, and flight-time results reproduce the saved manuscript inputs within an absolute tolerance of 1e-12. Tiny differences in affine-difference summaries arise from the historical intermediate CSV round trip (approximately 1e-16), not changed selections or formulas. The comparison uses absolute tolerance because some reported differences are themselves near zero.

Duplicate components are rejected; negative variances and zero denominators are flagged. No missing quantiles are imputed. Complete curves are required for HAM, while valid individual fitted quantiles remain in manuscript summaries. Historical convergence and raw-image processing are not newly certified.

## HAM

Seven helper function definitions and both focal VI formula expressions match V17 exactly. The cleaned map runner retains the original map calculations and uses explicit input paths. Unused VI formulas and the legacy pooled-flight command interface are omitted from this workflow. The 600-dpi layout renderer is unchanged.

All eleven maps were regenerated from their original TIFFs. Every saved VI array, rank array, mask, H2 curve, and mapped-H2 array matches exactly, including invalid pixels. Every RGB/HAM/colorbar PNG has identical decoded pixels. The 1,100-row H2-curve table matches byte-for-byte. One summary mean differs by 5.55e-17 under the tested Python libraries; the underlying arrays are identical.

Bundled curve/date/filename excerpts reproduce the same maps as the full source tables. They contain the 500 cotton and 600 maize curve rows, five and six date records, and six maize image-membership records required by Figure 5. Source hashes and row counts are recorded separately. Image inputs remain external and hash-pinned.

Within-image right-continuous ranks, nodata/black/nonfinite exclusions, interpolation, band formulas/scales, rotations, and colors are preserved. HAM remains a descriptive display of quantile-level heritability, not pixel-level heritability. It does not refit heritability or quantify pixel-level uncertainty.

Detailed checks and versioned verification scripts remain beside this local review package. Public hosting/citations are still required for large saved component inputs and the eleven TIFFs. Nothing was pushed to GitHub.
