# Preparation

Run the commands from the repository root. Replace the uppercase placeholders with your file paths and quote any path that contains spaces. Choose new output paths.

## Check field assignments

Compare the existing export with the corrected field map.
```sh
Rscript scripts/preprocess/join_field_info.R --input RAW.csv --field-map MAP.csv --output output/field_check --mode verify --allow-unmatched true
```
Use `--mode join` to write the joined analysis table. `--allow-unmatched true` lets the check proceed with missing assignments and records the affected images. The historical export had 309 unmatched images. The default is to stop if any image lacks a match.

## Recover archived BLUEs

The command extracts the saved BLUEs and checks that the export matches the reference byte for byte.
```sh
Rscript scripts/preprocess/recover_blue_reference.R --archive WORKSPACE.RData --reference BLUE_REFERENCE.csv --output output/recovery
```

## Build raw support

Calculate the mean and count of finite image quantiles for each genotype, date, and quantile in the 11 manuscript environments.
```sh
Rscript scripts/preprocess/build_raw_support.R --raw RAW.csv --predictors BLUE_REFERENCE.csv --output output/raw_support.csv
```

## Fit yield BLUEs

Fit the historical model in all 19 environments and compare the results with the manuscript yield input.
```sh
Rscript scripts/preprocess/fit_yield_blue.R --phenotype-2020=PHENO_2020.csv --phenotype-2021=PHENO_2021.csv --output=output/yield --reference=YIELD_REFERENCE.csv
```
The yield model uses REML with Pedigree as a fixed effect and Range, Pass, and Replicate as random effects. It multiplies yield in bushels per acre by 0.0673 to convert to tonnes per hectare. The output includes diagnostics and counts of excluded rows.
