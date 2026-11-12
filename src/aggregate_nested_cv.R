# Validate and aggregate completed primary V15 runs. No fitting; preserves all older code/results.
# Rscript Aggregate_Nested_CV_V15_V1.R <Production_Primary_V15> [NEW_output_dir]
# Default output: <root>/Combined_Results_V1. Existing directories are rejected.
# Reconstructs fingerprints in the recorded C locale, restoring the
# caller's locale even on failure. R serialization v3 includes native encoding;
# UTF-8/CP1252 headers otherwise change the hash despite identical model inputs.
# This does not rewrite checkpoints, replace their fingerprints, or weaken guards.
suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)
v15a_envs <- c("IAH4.2021", "MIH1.2020", "MNH1.2020", "MNH1.2021", "TXH1.2020",
  "TXH2.2021", "WIH1.2021", "WIH2.2020", "WIH2.2021", "WIH3.2020", "WIH3.2021")
v15a_require <- function(d, cols, label) {
  if (!all(cols %in% names(d))) stop(label, ": missing columns")
  if (anyNA(d[, ..cols])) stop(label, ": missing values")
}
v15a_key <- function(d, cols) do.call(paste, c(as.list(d[, ..cols]), sep="|"))
v15a_exact <- function(d, expected, cols, label) {
  v15a_require(d, cols, label)
  if (anyDuplicated(d[, ..cols])) stop(label, ": duplicate keys")
  if (!setequal(v15a_key(d, cols), v15a_key(expected, cols))) stop(label, ": incomplete or unexpected keys")
}
v15a_validate <- function(p, of, inf, q, env) {
  p <- as.data.table(p); of <- as.data.table(of); inf <- as.data.table(inf); q <- as.data.table(q)
  v15a_require(q, c("Env", "N_Genotypes"), "Input QA")
  if (nrow(q) != 1L || q$Env != env || q$N_Genotypes < 5L || q$N_Genotypes != round(q$N_Genotypes))
    stop("Unexpected input QA")
  v15a_require(of, c("Env", "Outer_Repeat", "Pedigree", "Outer_Test_Fold"), "Outer manifest")
  v15a_require(p, c("Env", "Scenario", "Outer_Repeat", "Outer_Fold", "Pedigree", "Obs", "Pred", "Selected_Scenario", "N_Train"), "Predictions")
  v15a_require(inf, c("Env", "Outer_Repeat", "Outer_Fold", "Pedigree", "Inner_Validation_Fold"), "Inner manifest")
  if (any(of$Env != env) || any(p$Env != env) || any(inf$Env != env)) stop("Environment label mismatch")
  ids <- sort(unique(of$Pedigree))
  if (length(ids) != q$N_Genotypes || any(!nzchar(ids))) stop("Outer cohort disagrees with QA")
  expected_outer <- CJ(Env=env, Outer_Repeat=1:5, Pedigree=ids)
  v15a_exact(of, expected_outer, c("Env", "Outer_Repeat", "Pedigree"), "Outer manifest")
  fold_sets <- of[, .(Complete=setequal(Outer_Test_Fold, 1:5)), by=Outer_Repeat]
  if (!all(fold_sets$Complete)) stop("Missing or unexpected outer folds")
  expected_pred <- CJ(Env=env, Scenario=c("Median", "all"), Outer_Repeat=1:5, Pedigree=ids)
  v15a_exact(p, expected_pred, c("Env", "Scenario", "Outer_Repeat", "Pedigree"), "Predictions")
  if (any(!is.finite(p$Obs)) || any(!is.finite(p$Pred)) || any(p$Selected_Scenario != p$Scenario))
    stop("Invalid predictions or selected scenario")
  if (any(p[, .(N_Obs=uniqueN(Obs)), by=Pedigree]$N_Obs != 1L)) stop("Outcome changes across repetitions/scenarios")
  aligned <- merge(p, of, by=c("Env", "Outer_Repeat", "Pedigree"))
  if (nrow(aligned) != nrow(p) || any(aligned$Outer_Fold != aligned$Outer_Test_Fold)) stop("Prediction/test-fold mismatch")
  test_sizes <- of[, .(N_Test=.N), by=.(Env, Outer_Repeat, Outer_Test_Fold)]
  sized <- merge(p, test_sizes, by.x=c("Env", "Outer_Repeat", "Outer_Fold"),
    by.y=c("Env", "Outer_Repeat", "Outer_Test_Fold"))
  if (any(sized$N_Train != length(ids)-sized$N_Test)) stop("Outer training sample-size mismatch")
  expected_inner <- rbindlist(lapply(1:5, function(k) {
    d <- copy(of[Outer_Test_Fold != k, .(Env, Outer_Repeat, Pedigree)])
    d[, Outer_Fold := k]; d
  }))
  v15a_exact(inf, expected_inner, c("Env", "Outer_Repeat", "Outer_Fold", "Pedigree"), "Inner training complement")
  inner_sets <- inf[, .(Complete=setequal(Inner_Validation_Fold, 1:5)), by=.(Outer_Repeat, Outer_Fold)]
  if (nrow(inner_sets) != 25L || !all(inner_sets$Complete)) stop("Missing or unexpected inner folds")
  paired <- merge(p[Scenario=="Median"], p[Scenario=="all"],
    by=c("Env", "Outer_Repeat", "Outer_Fold", "Pedigree"), suffixes=c("_Median", "_All"))
  if (nrow(paired) != length(ids)*5L || !identical(paired$Obs_Median, paired$Obs_All)) stop("Unpaired scenario predictions")
  invisible(TRUE)
}
v15a_meta <- function(meta, env, model_script) {
  if (!identical(meta$session$locale, "C"))
    stop("Expected the recorded C locale of the frozen V15 production run")
  locale_categories <- c("LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_NUMERIC", "LC_TIME")
  caller_locale <- setNames(vapply(locale_categories, Sys.getlocale, character(1)), locale_categories)
  on.exit({
    for (category in locale_categories) {
      restored <- Sys.setlocale(category, caller_locale[[category]])
      if (!identical(Sys.getlocale(category), caller_locale[[category]]))
        stop("Could not restore caller locale category: ", category)
    }
  }, add=TRUE)
  Sys.setlocale("LC_ALL", "C")
  if (!identical(Sys.getlocale(), "C")) stop("Could not establish recorded C locale for fingerprint verification")
  cfg <- meta$config
  if (!identical(normalizePath(cfg$script_path, winslash="/", mustWork=TRUE), model_script))
    stop("Manifest points to a different model script")
  pkg <- c(meta$session$otherPkgs, meta$session$loadedOnly)$glmnet
  r_version <- meta$session$R.version$version.string
  if (is.null(pkg$Version) || is.null(r_version) || length(meta$input_md5)!=3L || anyNA(meta$input_md5))
    stop("Missing recorded version/input provenance")
  if (!identical(names(meta$input_md5),c(cfg$yield,cfg$predictors,cfg$raw_support)) ||
      !is.character(cfg$raw_support) || length(cfg$raw_support)!=1L || !nzchar(cfg$raw_support))
    stop("Manifest must identify yield, BLUE and raw-support input hashes in order")
  # sessionInfo retains DESCRIPTION spelling (4.1-10); packageVersion used by
  # the original fingerprint normalizes that spelling to 4.1.10.
  glmnet_version <- as.character(package_version(pkg$Version))
  identity <- cfg; identity$max_outer_splits <- NULL; identity$output <- NULL
  fingerprint <- v15_hash_object(list(config=identity, input_md5=meta$input_md5,
    source_md5=tools::md5sum(cfg$script_path), glmnet=glmnet_version, R=r_version))
  if (!identical(fingerprint, meta$fingerprint)) stop("Saved fingerprint disagrees with configuration/versions/current frozen V15 source")
  # Hash before numeric comparisons. Production sort() returns an ALTREP lambda
  # wrapper; all.equal() on that wrapper can alter its internal serialization
  # state without changing values. Compare detached vectors, keeping the saved
  # config intact for repeat validation and provenance export.
  stopifnot(identical(cfg$envs, env), identical(cfg$scenarios, c("Median", "all")),
    identical(cfg$adaptive, FALSE), cfg$mode=="full", cfg$version=="V15", cfg$target_vi=="ExGR",
    cfg$outer_folds==5L, cfg$outer_repeats==5L, cfg$inner_folds==5L,
    cfg$max_outer_splits >= 25L,
    isTRUE(all.equal(cfg$alphas + 0, seq(0,1,.1))),
    isTRUE(all.equal(cfg$lambdas + 0, exp(seq(log(10),log(1e-4),length.out=100)))))
  identity$envs <- NULL
  list(common_config=identity, input_md5=meta$input_md5, glmnet=glmnet_version, R=r_version)
}
v15a_support <- function(d, prep, env, inner=FALSE) {
  d <- as.data.table(d); prep <- as.data.table(prep)
  keys <- c("Env","Scenario","Outer_Repeat","Outer_Fold",if(inner)"Inner_Fold","Feature")
  v15a_require(prep,c(keys,"Retained"),"Preprocessing support coverage")
  if(!is.logical(prep$Retained)||anyDuplicated(prep[,..keys])||any(prep$Env!=env))
    stop("Invalid preprocessing support keys/retention")
  group_keys <- setdiff(keys,"Feature")
  expected_groups <- CJ(Env=env,Scenario=c("Median","all"),Outer_Repeat=1:5,Outer_Fold=1:5)
  if(inner) expected_groups <- expected_groups[,.(Inner_Fold=1:5),by=.(Env,Scenario,Outer_Repeat,Outer_Fold)]
  v15a_exact(unique(prep[,..group_keys]),expected_groups,group_keys,"Preprocessing fold support coverage")
  v15a_require(d,c(keys,"Max_Absolute_Z","N_Outside_Training_Range","N_Applied","Extreme_Support_Flag"),"Input support diagnostics")
  if(!is.logical(d$Extreme_Support_Flag))stop("Support flags must be logical")
  expected <- prep[Retained==TRUE,..keys]
  v15a_exact(d,expected,keys,"Retained-feature support diagnostics")
  if(any(!is.finite(d$Max_Absolute_Z))||any(d$Max_Absolute_Z<0)||
     any(!is.finite(d$N_Applied))||any(d$N_Applied<1)||any(d$N_Applied!=round(d$N_Applied))||
     any(!is.finite(d$N_Outside_Training_Range))||any(d$N_Outside_Training_Range<0)||
     any(d$N_Outside_Training_Range>d$N_Applied)||
     any(d$N_Outside_Training_Range!=round(d$N_Outside_Training_Range)))stop("Invalid input support diagnostic values")
  if(any(d$Extreme_Support_Flag != (d$Max_Absolute_Z>100)))stop("Support flags disagree with declared threshold")
  if(any(d$Extreme_Support_Flag))stop("Extreme predictor support flags require review; aggregation refused")
  invisible(TRUE)
}
v15a_run <- function(root, output, model_script) {
  if (dir.exists(output) || file.exists(output)) stop("Output already exists; choose a new numbered directory")
  model_script <- normalizePath(model_script, winslash="/", mustWork=TRUE)
  source(model_script[1L])
  preds <- outer <- inner <- qa <- configs <- sources <- outer_support <- inner_support <- list(); common <- NULL
  for (env in v15a_envs) {
    path <- file.path(root, env)
    files <- file.path(path, c("Run_Status_V15.txt", "Run_Manifest_V15.rds", "Outer_Predictions_V15.csv",
      "Outer_Fold_Manifest_V15.csv", "Inner_Fold_Manifest_V15.csv", "Input_QA_V15.csv",
      "Outer_Input_Support_V15.csv","Inner_Input_Support_V15.csv",
      "Outer_Preprocessing_V15.csv","Inner_Preprocessing_V15.csv"))
    if (any(!file.exists(files))) stop("Missing environment inputs: ", env)
    hashes <- tools::md5sum(files)
    status <- readLines(files[1])
    if (!all(c("Status: COMPLETE", "Expected outer splits: 25", "Completed outer splits: 25", "Completed scenario fits: 50") %in% status))
      stop("Environment incomplete: ", env)
    meta <- readRDS(files[2]); identity <- v15a_meta(meta, env, model_script)
    if (is.null(common)) common <- identity else if (!identical(identity, common)) stop("Environment configurations/input hashes/versions differ")
    if (!paste("Fingerprint:", meta$fingerprint) %in% status) stop("Run status fingerprint mismatch")
    p <- fread(files[3]); of <- fread(files[4]); inf <- fread(files[5]); q <- fread(files[6])
    v15a_validate(p, of, inf, q, env)
    os <- fread(files[7]); ins <- fread(files[8]); op <- fread(files[9]); ip <- fread(files[10])
    v15a_support(os,op,env,FALSE); v15a_support(ins,ip,env,TRUE)
    if (!identical(hashes, tools::md5sum(files))) stop("Environment files changed while reading: ", env)
    preds[[env]] <- p; outer[[env]] <- of; inner[[env]] <- inf; qa[[env]] <- q; configs[[env]] <- meta
    outer_support[[env]] <- os; inner_support[[env]] <- ins
    sources[[env]] <- data.table(Env=env, File=files, MD5=unname(hashes))
  }
  if (!dir.create(output, recursive=TRUE, showWarnings=FALSE)) stop("Cannot create new output directory")
  p <- rbindlist(preds)
  v15_summarize(p, output, complete=TRUE)
  fwrite(rbindlist(outer), file.path(output,"Outer_Fold_Manifest_V15.csv"))
  fwrite(rbindlist(inner), file.path(output,"Inner_Fold_Manifest_V15.csv"))
  fwrite(rbindlist(qa), file.path(output,"Input_QA_V15.csv"))
  fwrite(rbindlist(outer_support),file.path(output,"Outer_Input_Support_V15.csv"))
  fwrite(rbindlist(inner_support),file.path(output,"Inner_Input_Support_V15.csv"))
  fwrite(rbindlist(sources), file.path(output,"Aggregation_Input_Hashes_V1.csv"))
  saveRDS(configs, file.path(output,"Environment_Run_Manifests_V15.rds"))
  s <- fread(file.path(output,"Environment_Scenario_Summary_V15.csv"))
  if (nrow(s)!=22L || any(s$N_Repeats!=5L)) stop("Unexpected environment summary coverage")
  paired <- merge(s[Scenario=="all"], s[Scenario=="Median"], by="Env", suffixes=c("_All","_Median"))
  paired[, `:=`(Delta_Pearson=Mean_Pearson_All-Mean_Pearson_Median,
    Delta_RMSE=Mean_RMSE_All-Mean_RMSE_Median,
    Pct_RMSE_Reduction=ifelse(Mean_RMSE_Median>0, 100*(Mean_RMSE_Median-Mean_RMSE_All)/Mean_RMSE_Median, NA_real_))]
  fwrite(paired, file.path(output,"Primary_Environment_Comparisons_V15.csv"))
  overall <- paired[, .(N_Environments=.N,
    N_Environments_With_Pearson=sum(is.finite(Delta_Pearson)),
    Environments_Higher_Pearson=sum(Delta_Pearson>0, na.rm=TRUE),
    Environments_Lower_RMSE=sum(Delta_RMSE<0), Mean_Delta_Pearson=mean(Delta_Pearson),
    Mean_Pct_RMSE_Reduction=mean(Pct_RMSE_Reduction))]
  fwrite(overall, file.path(output,"Primary_Comparison_Summary_V15.csv"))
  writeLines(c("Status: COMPLETE", "11 environments; primary Median versus all quantiles only; adaptive selection disabled.",
    "550 scenario prediction groups and 275 paired outer splits passed complete-cohort/fold/isolation checks.",
    paste("Outer prediction rows:",nrow(p)),
    "Primary metrics: pooled outer predictions within each repeat, then arithmetic mean of five repeat metrics.",
    "RMSE percent reduction: ratio of environment mean-repeat RMSEs, not mean of repeat-specific percent reductions.",
    "Overall changes: equal-weight means across environments. Undefined correlations remain NA; counts identify finite comparisons.",
    "Repeated-CV SDs are descriptive, not independent-repeat confidence intervals.",
    "Both inner/outer retained-feature support diagnostics are complete and contain no declared extreme-support flags; this is not a guarantee against ordinary extrapolation.",
    "Coefficients, tuning grids and preprocessing remain in the original environment directories."),
    file.path(output,"Combined_Status_V15.txt"))
  print(overall)
  invisible(overall)
}
