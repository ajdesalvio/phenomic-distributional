# Cotton models

Set the paths in `config/cotton.json`. Paths are relative to that file; keep the hashes to use the manuscript inputs. Required packages and tested versions are listed in `environment/R_versions.csv`. The functional calculation requires `clubSandwich` 0.7.0.

## Export the manuscript results

From the package root:

```sh
Rscript scripts/cotton/export_reference.R config/cotton.json output/cotton_reference
```

This checks the saved input hashes and regenerates the heritability summary tables and functional contrast/band tables without refitting. It requires the saved heritability metrics and three functional model files listed in the configuration.

## Fit heritability models

```sh
Rscript scripts/cotton/run_heritability.R config/cotton.json output/h2_smoke smoke
Rscript scripts/cotton/run_heritability.R config/cotton.json output/h2_full full
```

The smoke check fits nine models. The full run fits PSRI, CVI, and ExG2 across 101 dates and 100 quantiles: 30,300 models, sequentially. Each fits `response ~ 1 + (1 | Pedigree) + (1 | Block)` by REML. Block is Test × Tier. Heritability is `Vg / (Vg + Ve)`, conditional on block; no replication divisor is used. Singular fits are recorded, while unresolved convergence warnings stop the run.

Use a new output directory. Full outputs retain the original `_V6` filenames, including complete curves for figure preparation. Smoke results are explicitly incomplete. Heritability summaries aggregate separate within-date fits; their medians across quantiles are not heritability estimates for the median VI trait.

## Fit functional comparisons

```sh
Rscript scripts/cotton/run_functional.R config/cotton.json ExG2 output/fmm_ExG2 inspect
Rscript scripts/cotton/run_functional.R config/cotton.json ExG2 output/fmm_ExG2 full
```

Repeat for PSRI and CVI, each with its own output directory. `inspect` checks inputs and design without fitting or creating the output directory. Optional `smoke` runs the original reduced five-date/ten-quantile check with 199 draws; use a separate directory.

The full model is `response ~ Pedigree + Flight + Block + (1 | Plant)`. It retains 4,999 plant wild bootstrap-t draws and exports pointwise, per-curve, and joint 15-contrast bands for each VI. Covariance weights and smoothers stay fixed during bootstrap draws.

Append `resume` to the same full-run command to resume matching checkpoints. Use `validate` instead of `full` to check a completed cleaned run. Historical checkpoints cannot be adopted because source fingerprints differ. Keep the code, inputs, paths, configuration, and package versions unchanged while resuming.

These are additive genotype contrasts adjusted for flight and block, not genotype-by-date effects. See [verification and scope](../../provenance/cotton_cleanup.md) for the checks and inference qualifications.
