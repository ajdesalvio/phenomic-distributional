# Optional ExGR refitting

This estimates new BLUEs; it is not required to use the frozen manuscript inputs. The formula remains fixed Pedigree with random Range, Pass, and Replicate, separately by environment/date/quantile, using REML.

```sh
Rscript scripts/preprocess/plan_exgr.R --input=RAW.csv --output=output/plan --batch-size=300 --exclusions=config/exgr_exclusions.csv
Rscript scripts/preprocess/fit_exgr_batch.R --plan=output/plan/ExGR_Plan.rds --batch=1 --output=output/chunks
```

Run the worker for **each batch listed in `output/plan/Batches.csv`**. The supplied exclusion preserves MOH1.2020 DAP 13, absent from the historical BLUE archive. Workers checkpoint every 25 keys and validate saved checkpoints when resumed.

```sh
Rscript scripts/preprocess/export_exgr.R --plan=output/plan/ExGR_Plan.rds --chunks=output/chunks --output=output/rebuilt
Rscript scripts/preprocess/compare_exgr.R --rebuilt=output/rebuilt/Quantile_BLUEs_G2F_Distributional_ExGR.csv --reference=BLUE_REFERENCE.csv --status=output/rebuilt/Model_Status.csv --output=output/comparison
```

A comparison failure preserves the rebuilt estimates and discrepancy reports. The 10^-7 tolerance is strictly a numerical diagnostic threshold. Do not automatically replace manuscript inputs after refitting. See [optimizer sensitivity](../../provenance/optimizer_sensitivity.md).
