# V15 cleanup verification

The public runner preserves the manuscript primary comparison: Median versus all quantiles, no adaptive selection, five repeated five-fold outer splits, five inner folds, and the original tuning grid. It runs each environment separately in the C locale. Both choices preserve the original production fold assignments; running the old engine once over all environments would change its environment-index seed offset.

Verification against preserved V15 code and production outputs:

- All 21 retained engine function bodies and arguments are identical to the originals. The wrapper replaces machine-specific paths and generic defaults with explicit frozen inputs and production settings.
- Settings and every saved outer/inner fold assignment match for all 11 environments.
- One complete outer split for each primary scenario in TXH1.2020 was refitted using the full tuning grid. Original and cleaned outputs matched exactly; saved production fits matched within 1e-12.
- All six combined summary/prediction tables recalculated from saved predictions matched exactly.
- Synthetic checks passed for training-only preprocessing/tuning, constant-feature exclusion, intercept fallback, the public CLI, unchanged resumed checkpoints, source fingerprints, and rejection of modified inputs.

The aggregator retains the original validation and summary calculations. Its source loader reads the engine from the first entry of the source-path vector; checkpoint provenance now fingerprints both engine and runner. Original production checkpoints therefore require a separate output directory from cleaned runs.

No full manuscript CV rerun was performed. These checks establish unchanged definitions, folds, summaries, and representative fits, rather than claiming a fresh complete rerun. Original scripts, manuscript inputs, and production outputs were preserved. Tested package versions are in `environment/R_versions.csv`; source mappings are in `source_map.csv`.

For a small synthetic check from the package root:

```sh
Rscript validation/check_nested_cv.R output/nested_synthetic_check
```

Frozen input files are not bundled. Their public hosting and release citation remain necessary before this package is ready for turnkey public reproduction.
