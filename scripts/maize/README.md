# Within-environment nested CV

Install `glmnet`, `data.table`, and `jsonlite`. Set the three input paths in `config/nested_cv.json`; paths are relative to that file. Keep the recorded hashes to use the manuscript inputs.

From the package root:

```sh
Rscript scripts/maize/run_nested_cv.R config/nested_cv.json output/nested_cv
Rscript scripts/maize/aggregate_nested_cv.R output/nested_cv output/nested_summary
```

The first command runs the 11 environments sequentially. To run one environment, append its name, for example `TXH1.2020`. Use a new output directory; repeating the same command resumes matching checkpoints. Original production checkpoints are not reused by the cleaned code because source fingerprints differ. Keep the scripts and R environment unchanged while resuming.

The primary comparison is **Median versus all 100 quantiles**, with adaptive selection disabled. Each environment uses five outer folds repeated five times, five inner folds, seed 123, the original alpha/lambda grid, and training-only preprocessing/support screening. Environments remain separate runs in the C locale, preserving the manuscript's fold assignments.

Outputs retain their `_V15` names for compatibility. The aggregator requires all 11 completed environments, validates cohort/fold/support coverage, and reproduces the original repeat-wise summaries. Detailed tuning, coefficients, and preprocessing remain in each environment's directory.

This analysis starts from frozen genotype BLUEs. First-stage spatial estimation and historical ExGR selection are outside these CV folds. Repeated-CV standard deviations are descriptive. See `provenance/nested_cv_cleanup.md` for cleanup verification.
