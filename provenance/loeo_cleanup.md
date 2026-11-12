# FPCA/LOEO cleanup

The manuscript primary reference comprises 22 historical Median/all WITHIN models. The 88 V4 trimmed-scenario reruns are not its source. No production models were refitted or replaced during cleanup.

## Preserved analysis

The four LOEO functions for quantile parsing, predictor assembly, within-environment scaling, and fitting have identical arguments and bodies to V4. The public runner exposes only Median/all and WITHIN, with one worker. It adds explicit input paths, content checks, and source/configuration fingerprints. The primary quantile sets are unaffected by the historical trimming-parser error.

FPCA retains the V2 training-only sparse smooth mean/covariance model, GCV bandwidths, score estimation, and training-time-range restriction. The wrapper enables exports and all 100 quantiles, removes machine paths/manual part loading, and fails rather than silently skipping a failed basis. Full genotype/environment IDs are restored by matching input keys. The separate archived-part assembler restores environment suffixes without truncating dotted genotype names.

## Verification

- All 22 historical primary models reproduce their archived predictions within 1e-10. Recomputed correlations, errors, sample sizes, outcome/prediction SDs, and raw nonzero-coefficient counts match manuscript source tables within 1e-12.
- The public reference-export command produces the same verified metrics.
- For TXH1.2020, both Median and all-quantile held-out predictor matrices reconstructed from frozen projections match the archived scaled matrices within 1e-12.
- All 7,135,700 archived projection rows reproduce the combined table exactly after the documented environment-ID restoration and missing-column fill.
- Original and cleaned FPCA scores, FVE values, and IDs match exactly on a synthetic three-environment example.
- A synthetic Median fit using the full historical tuneLength=25 search gives identical tuning, predictions, coefficients, folds, and metrics, excluding elapsed time.
- Parser and environment-scaling safeguards were checked. Partition creation/reuse checks file content hashes.

The local numbered review directory contains the detailed checks and verification scripts. Source mappings and tested package versions accompany this package.

## Limits

The archived projection table is traced to PART1-PART5, but the saved historical V2 script had a median-only filter and disabled exports. We have not established a complete original execution record for all bases, or refitted every basis from recovered BLUEs. Fresh FPCA is therefore an optional reconstruction, not a certified byte-identical source of the frozen table. Numerical parity of a full fresh primary LOEO refit is likewise not claimed.

The original statistical qualifications remain: first-stage BLUEs and historical VI selection precede CV; outer-training FPCA is reused within tuning; WITHIN scaling uses the unlabeled held-out predictor cohort; caret can interpolate/clamp the requested penalty along its fitted path. The cleanup does not change those choices or promote this analysis to fully nested FPCA validation.

The projection table, five parts, and historical model archive need public hosting/release citations before turnkey public reproduction can be claimed. No GitHub push was performed as part of this local review.
