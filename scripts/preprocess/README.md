# Preparation

Run from the package root. Replace uppercase paths with your files; quote paths containing spaces. Output directories/files must be new.

**Field assignments** — verify the existing export against the corrected map:
```sh
Rscript scripts/preprocess/join_field_info.R --input RAW.csv --field-map MAP.csv --output output/field_check --mode verify --allow-unmatched true
```
Use `--mode join` to write the joined analysis table. The explicit unmatched option preserves the historical 309 unmatched images and records them; otherwise unmatched images fail validation.

**Archived BLUEs** — recover and require byte identity with the reference:
```sh
Rscript scripts/preprocess/recover_blue_reference.R --archive WORKSPACE.RData --reference BLUE_REFERENCE.csv --output output/recovery
```

**Raw support** — calculate finite image-quantile means and counts for the 11 manuscript environments:
```sh
Rscript scripts/preprocess/build_raw_support.R --raw RAW.csv --predictors BLUE_REFERENCE.csv --output output/raw_support.csv
```

**Yield BLUEs** — fit the historical 19-environment model and compare with the current yield input:
```sh
Rscript scripts/preprocess/fit_yield_blue.R --phenotype-2020=PHENO_2020.csv --phenotype-2021=PHENO_2021.csv --output=output/yield --reference=YIELD_REFERENCE.csv
```
The yield model remains fixed Pedigree with random Range, Pass, and Replicate, fitted by REML. The conversion remains 0.0673 from bu/acre to t/ha. Diagnostics and excluded-row counts are saved.
