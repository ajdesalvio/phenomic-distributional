# FPCA/LOEO cleanup

The manuscript’s primary across-environment results come from 22 historical models comparing Median and all quantiles with WITHIN scaling. The 88 V4 reruns cover trimmed scenarios. The review checked saved production models and used synthetic fits to test changes to the code.

## Preserved analysis

The four V4 functions parse quantiles, assemble predictors, scale within environments, and fit models. Their arguments and bodies match V4 exactly. The current runner exposes the two primary scenarios with WITHIN scaling and one worker. It takes explicit input paths, checks file contents, and records checksums for source code and settings. Both primary quantile sets were unaffected by the historical trimming parser error.

The FPCA code keeps V2’s sparse smoothing of the training mean and covariance, GCV bandwidths, score estimation, and restriction to the training time range. The wrapper runs all 100 quantiles, exports its results, and stops if a basis fails. It matches input keys to recover full genotype and environment IDs. The projection assembly command restores the environment suffixes without truncating genotype names that contain periods.

## Verification

- All 22 historical models reproduce archived primary predictions within 1e-10. Recalculated correlations, errors, sample sizes, observed and predicted standard deviations, and nonzero coefficient counts match the manuscript tables within 1e-12. The reference export produces the same verified metrics.
- The Median and all-quantile held-out matrices for TXH1.2020, rebuilt from the saved projections, match the archived scaled matrices within 1e-12.
- All 7,135,700 projection rows reproduce the combined table exactly after environment IDs are restored and missing columns filled.
- The original and current FPCA code give identical scores, fractions of variance explained, and IDs in a synthetic three-environment example.
- A synthetic Median fit with the full historical `tuneLength=25` search gives identical tuning, predictions, coefficients, folds, and metrics. Elapsed time was excluded.
- The review checked the parser and the safeguards for environment scaling. Partition creation and reuse check file contents by checksum.

The detailed checks and scripts are kept in the local review directory. See [source mappings](source_map.csv) and [tested R versions](../environment/R_versions.csv) for the recorded software.

## Limits

The projection table traces back to PART1–PART5. The saved V2 script filtered to the median and had exports disabled, so the original execution record is incomplete for the full set of bases. The FPCA command reconstructs the historical method, and the tests cover representative fits and recovery of the archived projections. Exact reproduction by a full fresh FPCA and LOEO run remains unverified.

BLUE estimation and vegetation index selection took place before CV. Inner tuning reuses FPCA estimated in the outer training set, which limits the degree of nesting. WITHIN scaling uses predictors from the unlabeled held-out cohort. Caret can interpolate or clamp a requested penalty to its fitted path. See the [release guide](../data/DOWNLOADS.md) for the projection table, the five parts, and the saved model archive.
