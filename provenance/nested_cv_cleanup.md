# V15 cleanup verification

The runner compares Median and all 100 quantiles with adaptive selection disabled. It uses five repeats of five-fold outer CV, five inner folds, and the original tuning grid. It runs each environment separately in the C locale to reproduce the manuscript fold assignments. Running the original engine over all environments at once would change the seed offsets.

- The arguments and bodies of all 21 retained engine functions match the originals. The wrapper supplies explicit paths for the frozen inputs and uses the production settings.
- The settings and every saved outer and inner fold assignment match across all 11 environments.
- In TXH1.2020, we refit one complete outer split for each scenario using the full tuning grid. The original and current code agree exactly. The fits match the saved production results within 1e-12.
- All six combined summary and prediction tables match exactly when recalculated from saved predictions.
- Synthetic checks covered training-only preprocessing and tuning, constant-feature removal, intercept-only fallback, the command-line runner, checkpoint resume, source code checksums, and rejection of modified inputs.

The aggregator uses the original validation and summary calculations. It loads the engine from the first entry in the source path list. The run records now include checksums for both the engine and the runner. Use a separate output directory for new runs because the original checkpoints record different source files.

The checks cover the function definitions, fold assignments, summaries, and representative fits. A full production rerun was outside the scope of the review. See the [tested R versions](../environment/R_versions.csv), [source mappings](source_map.csv), and [release guide](../data/DOWNLOADS.md).

Run this small synthetic check from the repository root.

```sh
Rscript validation/check_nested_cv.R output/nested_synthetic_check
```
