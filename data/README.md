# Inputs

Large data are not bundled in this review package. `published_inputs.csv` lists the exact manuscript filenames and MD5 hashes; the configuration contains the same identities. Relative paths resolve from the configuration directory.

The frozen BLUEs and support table are the manuscript reference. The combined yield file includes flowering traits, but `fit_yield_blue.R` regenerates only yield. Archive recovery additionally requires the saved `20250821_Elastic_Net_V9_Home.RData` workspace containing `blue.filt`.

The field map is the corrected V2 map. Field joining starts from existing image-level quantiles and does not verify image extraction. Raw support is an unweighted mean of finite image-level quantiles, with finite counts, for each genotype/date/quantile.

The existing support provenance records its original builder hash. Relocating these files is supported; silently substituting a newer input or builder is not.

Before publication, supply a durable download location and source-release citation for these data. The current local files and checksums are verified; public data hosting is not yet established.

## Across-environment inputs

`config/loeo.json` pins the historical WITHIN model archive, the combined ExGR projection table, yield inputs, and recovered BLUEs. The model archive is about 3.2 GB; projections are about 4.1 GB. Place the archive and yield/BLUE files in `data/` or update the configured paths. The partition command takes the projection CSV path directly. Archived PART1-PART5 are needed only for optional assembly verification. These large files are not bundled and still need public download locations. Small verified primary metrics are in `results/reference/LOEO_primary_metrics.csv`.

## Cotton inputs

`config/cotton.json` records the field map and three original RData inputs for PSRI, CVI, and ExG2. Place them in `data/` or change their configured paths. The saved-result exporter additionally needs `Fit_Metrics_V6.csv` and the three `functional_model_V7.rds` files under `data/cotton_reference/` as configured. These files are not bundled and still need public hosting. The small Table 2 and VI summary reference tables are included in `results/reference/`.

## Maize heritability and HAM

The saved component CSV and flight-time table are configured in `config/maize_heritability.json`. HAM's small verified curve/date/image-membership excerpts are bundled as `ham_*.csv`; `config/ham.json` pins their hashes and the eleven external source TIFFs. See the two run guides for commands. Large component/image files still require public download locations.

## Remaining figure inputs

`config/figures.json` lists and hashes the saved plotting inputs. Place them in `data/figures/` or edit the paths. These external inputs are not bundled. Table 1 and S1 metadata are included in `metadata/`; their document provenance is recorded separately.

`configured_inputs.csv` covers all configured files, including HAM images. `published_inputs.csv` covers only the preparation validator. Optional command-line inputs (raw phenotype years, BLUE workspace archive, archived projection parts, and saved V15 run directories) are described in the corresponding guides and also need release locations. Hash-pinned plotting inputs are frozen artifacts: regenerated tables may be numerically equal without byte-identical hashes. Check values before changing hashes.

[Source imagery collections](G2F_COLLECTIONS.md) provides the relevant G2F DOIs. These links do not replace the processed-file downloads listed above.
