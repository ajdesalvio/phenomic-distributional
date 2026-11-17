# Within-environment nested CV

Install `glmnet`, `data.table`, and `jsonlite`. The release places the three inputs at the paths in `config/nested_cv.json`. If you store them elsewhere, edit the paths relative to that configuration file and keep the manuscript checksums.

From the package root:

```sh
Rscript scripts/maize/run_nested_cv.R config/nested_cv.json output/nested_cv
Rscript scripts/maize/aggregate_nested_cv.R output/nested_cv output/nested_summary
```

The first command runs all 11 environments in order. To run one environment, append its name, such as `TXH1.2020`. Start in a new output directory. To resume, repeat the command with the same directory, scripts, and R environment. Checkpoints from the original production runs belong to a different version of the code.

The primary comparison is **Median versus all 100 quantiles**, with adaptive selection disabled. Within each environment, the runner repeats five-fold outer CV five times and uses five inner folds for tuning. It uses seed 123 and the original alpha and lambda grid. Each split estimates preprocessing and screens raw support using its training genotypes. Running each environment separately in the C locale reproduces the manuscript’s fold assignments.

Outputs retain their `_V15` names for compatibility. Once all 11 environments finish, the aggregator checks genotype coverage, fold assignments, and support diagnostics, then recalculates the summaries for each repeat. Look in each environment’s directory for tuning results, coefficients, and preprocessing details.

The CV starts from saved genotype BLUEs. Spatial adjustment and the historical selection of ExGR took place before CV. The standard deviations describe variation across repeats. See the [validation notes](../../provenance/nested_cv_cleanup.md).
