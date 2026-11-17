# Inputs

The [release guide](DOWNLOADS.md) lists the archives and explains where to extract them. The journal DOI is pending. The [G2F collections](G2F_COLLECTIONS.md) list source imagery separately.

The [configured input list](configured_inputs.csv) gives the expected paths and checksums, including the HAM images. [Published inputs](published_inputs.csv) covers the files used by the preparation validator. Relative paths in a configuration are resolved from that configuration file’s directory.

## Field assignments, BLUEs, and raw support

The manuscript uses the corrected V2 field map and the archived BLUEs and raw support table. The field join starts from existing image quantiles. For each genotype, date, and quantile, the support table records the unweighted mean and count of finite image measurements. The provenance file records the checksum of the script that built this table. You can move the files and update the paths while keeping the original checksums.

The reference yield file includes flowering traits. The yield fitting script regenerates yield only. To recover the archived BLUEs, use the `20250821_Elastic_Net_V9_Home.RData` workspace containing `blue.filt`. See the [preparation guide](../scripts/preprocess/README.md).

## Across-environment inputs

The [LOEO configuration](../config/loeo.json) records the historical WITHIN model archive (about 3.2 GB), yield data, and recovered BLUEs. It also records the checksum of the combined ExGR projection table (about 4.1 GB). Pass that table’s path directly to the partitioning command. The five archived parts are available if you want to check how the table was assembled. See the [LOEO guide](../scripts/maize/LOEO.md) and the [reference metrics](../results/reference/LOEO_primary_metrics.csv).

## Cotton inputs

The [cotton configuration](../config/cotton.json) lists the field map and the three original RData inputs for PSRI, CVI, and ExG2. To export saved results, use `Fit_Metrics_V6.csv` and the three `functional_model_V7.rds` files under `data/cotton_reference/`. Table 2 and the vegetation index summaries are in `results/reference/`. See the [cotton guide](../scripts/cotton/README.md).

## Maize heritability and HAM

The [maize heritability configuration](../config/maize_heritability.json) lists the saved components and flight times. The [HAM configuration](../config/ham.json) gives the checksums for the included `ham_*.csv` excerpts and the eleven TIFFs in the release. See the [heritability guide](../scripts/maize/HERITABILITY.md) and [heatmap guide](../scripts/figures/HAM.md).

## Remaining figure inputs

The [figure configuration](../config/figures.json) lists the saved tables and checksums under `data/figures/`. Table 1 and S1 metadata are in `metadata/`, with sources recorded in [figure and table provenance](../provenance/authored_artifacts.json).

A regenerated table can match the published numbers and still have a different checksum because of formatting. Check the values before updating a checksum for a new analysis. The release also keeps the optional phenotype files, BLUE workspace, projection parts, and saved V15 runs separate from the main figure inputs.
