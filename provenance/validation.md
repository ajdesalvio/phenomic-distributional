# Data preparation validation

- The field join and input validation checks cover duplicate and missing keys, overwrite protection, and the raw support provenance.
- Rebuilding raw support produced 11,684,000 cells and matched the manuscript input byte for byte (MD5 `a0c8cad4ecfffcc3fb184407f57efe20`).
- The yield refit reproduced the reference estimates for all 19 environments. Small tests also cover saving, resuming, merging, and comparing ExGR fits.

Run the preparation checks with `Rscript validation/check_preprocessing.R output/checks` from the repository root. See the [release guide](../data/DOWNLOADS.md) for the inputs and the [optimizer sensitivity note](optimizer_sensitivity.md) for what to expect from ExGR refits.
