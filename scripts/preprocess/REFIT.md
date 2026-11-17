# Optional ExGR refitting

Use this workflow to estimate new ExGR BLUEs. For manuscript reproduction, use the archived BLUEs in the data release. The refit uses REML and models Pedigree as a fixed effect and Range, Pass, and Replicate as random effects, separately for each environment, date, and quantile.

```sh
Rscript scripts/preprocess/plan_exgr.R --input=RAW.csv --output=output/plan --batch-size=300 --exclusions=config/exgr_exclusions.csv
Rscript scripts/preprocess/fit_exgr_batch.R --plan=output/plan/ExGR_Plan.rds --batch=1 --output=output/chunks
```

Run the worker for **each batch listed in `output/plan/Batches.csv`**. The supplied exclusion list omits DAP 13 in MOH1.2020 because that date is absent from the historical BLUE archive. Each worker saves a checkpoint every 25 model keys and checks it when resuming.

```sh
Rscript scripts/preprocess/export_exgr.R --plan=output/plan/ExGR_Plan.rds --chunks=output/chunks --output=output/rebuilt
Rscript scripts/preprocess/compare_exgr.R --rebuilt=output/rebuilt/Quantile_BLUEs_G2F_Distributional_ExGR.csv --reference=BLUE_REFERENCE.csv --status=output/rebuilt/Model_Status.csv --output=output/comparison
```

If the comparison fails, the rebuilt estimates and discrepancy reports are still available for inspection. The tolerance of 10^-7 is for numerical comparison. Keep refits separate from manuscript inputs. See the [optimizer sensitivity audit](../../provenance/optimizer_sensitivity.md).
