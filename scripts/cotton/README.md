# Cotton models

The [data release](../../data/DOWNLOADS.md) places the inputs at the paths in `config/cotton.json`. If you store them elsewhere, change the paths relative to that file and keep the manuscript checksums. See `environment/R_versions.csv` for the tested package versions. The functional workflow requires `clubSandwich` 0.7.0.

## Export the manuscript results

From the package root:

```sh
Rscript scripts/cotton/export_reference.R config/cotton.json output/cotton_reference
```

The exporter checks the saved heritability metrics and three functional model files against their recorded checksums. It uses those files to rebuild the heritability summaries, functional contrasts, and confidence bands.

## Fit heritability models

```sh
Rscript scripts/cotton/run_heritability.R config/cotton.json output/h2_smoke smoke
Rscript scripts/cotton/run_heritability.R config/cotton.json output/h2_full full
```

The smoke check fits nine models. The full run fits 30,300 models in sequence across PSRI, CVI, and ExG2, 101 dates, and 100 quantiles. Each uses REML with `response ~ 1 + (1 | Pedigree) + (1 | Block)`, where Block is Test × Tier. Heritability is `Vg / (Vg + Ve)`, conditional on block. The script records singular fits and stops on unresolved convergence warnings.

Use a new output directory. Use the full run’s outputs, including the complete curves, to prepare figures. When a summary reports the median across quantiles, it takes the median of the heritability estimates from those separate within-date fits. Heritability of the median vegetation index trait comes from the fit at quantile 0.50.

## Fit functional comparisons

```sh
Rscript scripts/cotton/run_functional.R config/cotton.json ExG2 output/fmm_ExG2 inspect
Rscript scripts/cotton/run_functional.R config/cotton.json ExG2 output/fmm_ExG2 full
```

Repeat for PSRI and CVI, giving each index its own output directory. `inspect` checks the inputs and design before fitting. For a smaller test, run `smoke` in a separate directory. It uses five dates, ten quantiles, and 199 draws.

The full model is `response ~ Pedigree + Flight + Block + (1 | Plant)`. It uses 4,999 plant wild bootstrap-t draws and exports pointwise bands, bands for each curve, and bands adjusted jointly over 15 contrasts for each vegetation index. Covariance weights and smoothers stay fixed during bootstrap draws.

Append `resume` to the full-run command to continue from compatible checkpoints created with the current code. Keep the code, inputs, paths, configuration, and package versions the same. Use `validate` in place of `full` to check a completed run.

The additive model estimates genotype contrasts across flights, with adjustment for flight and block. See the [validation and interpretation notes](../../provenance/cotton_cleanup.md).
