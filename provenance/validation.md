# Validation

The review copies preserve the statistical calculations. Formatting was checked for identical parsed R expressions.

- Fifteen fixture checks passed, including duplicate/missing-key and overwrite guards.
- Field-join and raw-support calculation bodies match the verified originals.
- Full raw-support regeneration produced 11,684,000 cells byte-identical to the manuscript support table (MD5 `a0c8cad4ecfffcc3fb184407f57efe20`).
- Cleaned/original yield code produced byte-identical synthetic estimates; the source version previously reproduced all 19 historical environments.
- ExGR checkpoint resume, merge, and comparison passed on saved raw-fit chunks (four models, 1,350 rows).
- The current frozen-input configuration passed hash and support-lineage validation.

No manuscript models were reanalysed during cleanup. The separate full ExGR refit and its numerical qualification are described in [the audit](optimizer_sensitivity.md). Public download locations remain to be supplied.

Run the portable fixture checks with `Rscript validation/check_preprocessing.R output/checks` from the package root.
