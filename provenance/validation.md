# Validation

The current preparation scripts preserve the original statistical calculations. We checked that formatting edits produced identical parsed R expressions.

- Fifteen tests with small inputs checked duplicate and missing keys, overwrite protection, and other input safeguards.
- The functions for field joins and raw support match the verified originals.
- Rebuilding raw support produced 11,684,000 cells and matched the manuscript table byte for byte (MD5 `a0c8cad4ecfffcc3fb184407f57efe20`).
- The original and current yield scripts produced identical output on synthetic data. The original script also reproduced the reference yield estimates for all 19 historical environments.
- We checked ExGR checkpoint resume, merge, and comparison using saved fit chunks covering four models and 1,350 rows.
- The manuscript inputs passed checksum and raw support provenance checks.

The preparation review used saved results and small tests. See the [optimizer sensitivity audit](optimizer_sensitivity.md) for the separate full ExGR refit and the [release guide](../data/DOWNLOADS.md) for its inputs.

From the repository root, run `Rscript validation/check_preprocessing.R output/checks` for the small tests.
