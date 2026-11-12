# Cotton cleanup verification

The source models are full-season block-only heritability V6, its V5 formula/input helpers, and functional comparison V7. Their numbered originals and all manuscript outputs remain untouched.

## Changes

- Shared numerical functions were extracted without changing their bodies or arguments: seven heritability functions and fifteen functional functions.
- Small command-line runners replace local paths and superseded orchestration. Heritability runs sequentially and writes separate date checkpoints; the original V5 replay is no longer a runtime dependency.
- Functional setup uses installed R libraries rather than a drive-specific dependency directory. The original package checks, random-number generator, defaults, and model/inference routines remain. Provenance now fingerprints the engine and input helper as well as the runner.
- A reference exporter regenerates tables from pinned saved estimates without fitting. Heritability summary expressions preserve the original manuscript table calculations.

## Evidence

- All 22 extracted function definitions match their source exactly.
- The public heritability smoke command validates all three original long/wide inputs and fits nine models. Vg, Vb, Ve, and heritability agree with production within 1e-10; cohort hashes match exactly.
- All 30,300 saved heritability estimates satisfy the original formula. All six heritability summary tables reproduce manuscript tables within 1e-12.
- Functional scientific defaults match production for all three VIs. The ExG2 input preparation reproduces the complete cohort and design: 47,807 plant-date observations, 479 plants, 101 dates, 16 genotypes, ten blocks, and fixed-design rank 125. The public inspection command also passes.
- All functional contrast/band and adjusted-mean tables reproduce saved outputs within 1e-12. Smoothing, propagated covariance, and standard errors recalculated from all 300 saved point fits match within 1e-10.
- All 4,999 bootstrap weight vectors per VI match exactly. Replaying the first three draws per VI matches saved bootstrap-t statistics within 1e-10. Joint and per-curve critical values reproduce the saved draw summaries.
- A synthetic repeated-measures fit gives identical original/cleaned coefficients, variance components, CR2 calculations, and diagnostics, excluding elapsed time.

Detailed checks and verification scripts are retained beside this numbered local package. This was bounded verification, not a full 30,300-model or 4,999-draw rerun.

## Interpretation retained

Cotton heritability is Vg/(Vg+Ve), conditional on block. No heritability confidence intervals or genotype-by-date variance components are added. Test labels identify field sections; the planned irrigation treatment was not implemented.

Functional bands are approximate, conditional on fitted covariance weights and smoothers. They account for dependence within plants, assume independence between plants after adjustment, and do not propagate variance-component or smoothing-parameter uncertainty. Multiplicity adjustment covers fifteen TM-1 contrasts over the quantile grid separately for each VI, not all VIs jointly. The additive model does not establish date-specific genotype differences in senescence.

The RData inputs contain the previously prepared image quantiles; this cleanup does not certify raw-image extraction. Large inputs and saved reference models still need public hosting and release citations. No GitHub push was performed.
