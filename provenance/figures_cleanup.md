# Correlation, figure, and table cleanup

This numbered package preserves the final numerical expressions and plotting specifications while removing local-drive paths, old manuscript/image dependencies, retired scenario panels, and document-editing steps. It does not rerun CV, FPCA, mixed models, or bootstrap inference.

## Verification

- All 192 environment-date correlation matrices, pair counts, retained-quantile masks, screening diagnostics, ranges, and 120,000 plotted correlations reproduce the saved V5 calculation exactly. The V6 renderer retains the full q=1 row and column in every panel.
- Feature importance rebuilt from the archived V15 outer coefficients matches the saved Figure 7/supporting data exactly. Its normalization and averaging are unchanged.
- Figure 3 contrasts/bands, Figure 4 component/season summaries, supporting flight-time statistics, and Figure 6 panel values match their saved plotting inputs within absolute 1e-12.
- Both S3 CSVs reproduce the final files byte-for-byte, including equal-environment averaging and per-environment percentage reductions.
- All 30 rendered figures were compared with the final 600-dpi reference PNGs. The original Figure 4 PDF null device and prediction PNG null device are selected separately to preserve text/layout metrics. Pixel comparison, resolution/dimension checks, and visual contact-sheet inspection are recorded in the local verification directory.

The renderer includes Figures 3, 4, 6–8 and S1–S25. Figure 5 remains in the previously verified HAM workflow. Figure 1/2 source presentations are preserved authored artifacts, not regenerated statistical results.

## Metadata and limitations

Table 1 model scope and Table S1 genotype metadata were transcribed from the identified accepted manuscript and supplement. Deleted revisions were excluded and vertically merged labels repeated. They are curated metadata, not new calculations; source hashes are recorded in `authored_artifacts.json`.

The available Figure 1 reference source still needs the author's date-label correction before release. Figure 2's corrected editable source is preserved. No claims are made about a fresh reconstruction of their original photography or raw quantile artwork.

Scientific qualifications remain unchanged: descriptive quantile correlations, approximate conditional functional bands, saved-component maize heritability, and historical primary LOEO fits. Physical figure sizes, uncertainty meanings, colors, and final numbering are preserved. Large inputs still require public hosting/release citations. No GitHub push or manuscript changes were made.
