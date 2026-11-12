# V15 analysis engine. Public entry point: scripts/maize/run_nested_cv.R.
# Statistical function bodies are preserved from Elastic_Net_V15.R.
v15_scenarios <- function() {
    list(Median = 50L, all = 1:100, drop_01_1 = 2:99, drop_01_02_99_1 = 3:98, drop_01_03_98_1 = 4:97, 
        drop_01_05_96_1 = 6:95, drop_01_10_91_1 = 11:90, drop_01_15_86_1 = 16:85, drop_01_20_81_1 = 21:80, 
        drop_01_25_76_1 = 26:75)
}

v15_hash_object <- function(x) {
    f <- tempfile("v15_fingerprint_", fileext = ".rds")
    on.exit(unlink(f), add = TRUE)
    saveRDS(x, f, compress = FALSE, version = 3)
    unname(tools::md5sum(f))
}

v15_atomic_rds <- function(x, path) {
    if (file.exists(path)) 
        stop("Refusing to overwrite completed checkpoint: ", path)
    tmp <- paste0(path, ".partial")
    saveRDS(x, tmp, compress = FALSE)
    if (!file.rename(tmp, path)) 
        stop("Could not finalize checkpoint: ", path)
}

v15_fold_ids <- function(ids, k, seed) {
    if (anyNA(ids) || anyDuplicated(ids)) 
        stop("Fold IDs must be unique nonmissing genotypes")
    if (length(ids) < k || k < 2L) 
        stop("Insufficient genotypes for requested folds")
    set.seed(seed)
    ord <- order(ids)
    out <- integer(length(ids))
    out[ord] <- sample(rep(seq_len(k), length.out = length(ids)))
    out
}

v15_numeric_constant <- function(z) {
    z <- z[is.finite(z)]
    if (length(z) < 2L) 
        return(TRUE)
    n <- length(z)
    mu <- mean(z)
    variance <- mean((z - mu)^2)
    eps <- .Machine$double.eps
    variance <= n * eps * variance + (n * mu * eps)^2
}

v15_support_stats <- function(z) {
    z <- sort(z[is.finite(z)])
    n <- length(z)
    if (!n) 
        return(c(N = 0, Resolution = NA, Unique = 0, UniquePct = 0, FrequencyRatio = Inf, DominantCount = 0, 
            Drop = 1))
    tolerance <- n * .Machine$double.eps * max(abs(z))
    u <- unique(z)
    counts <- tabulate(match(z, u), nbins = length(u))
    if (length(u) > 1L && any(diff(u) <= tolerance)) {
        groups <- integer(length(u))
        g <- 1L
        anchor <- u[1L]
        groups[1L] <- g
        for (i in seq.int(2L, length(u))) {
            if (u[i] - anchor > tolerance) {
                g <- g + 1L
                anchor <- u[i]
            }
            groups[i] <- g
        }
        counts <- as.numeric(rowsum(counts, groups, reorder = FALSE))
    }
    counts <- sort(counts, decreasing = TRUE)
    unique_pct <- 100 * length(counts)/n
    ratio <- if (length(counts) < 2L) 
        Inf
    else counts[1L]/counts[2L]
    drop <- v15_numeric_constant(z) || (ratio > 19 && unique_pct <= 10)
    c(N = n, Resolution = tolerance, Unique = length(counts), UniquePct = unique_pct, FrequencyRatio = ratio, 
        DominantCount = counts[1L], Drop = as.numeric(drop))
}

v15_preprocess_fit <- function(x, raw_support) {
    x <- as.matrix(x)
    raw_support <- as.matrix(raw_support)
    storage.mode(x) <- storage.mode(raw_support) <- "double"
    if (is.null(colnames(x)) || anyDuplicated(colnames(x))) 
        stop("Unique feature names required")
    if (!identical(dim(x), dim(raw_support)) || !identical(dimnames(x), dimnames(raw_support))) 
        stop("Raw support rows/features must exactly align with BLUE matrix")
    x[!is.finite(x)] <- NA_real_
    med <- apply(x, 2L, function(z) if (all(is.na(z))) 
        NA_real_
    else stats::median(z, na.rm = TRUE))
    observed <- colSums(!is.na(x))
    fill <- x
    for (j in seq_len(ncol(fill))) fill[is.na(fill[, j]), j] <- med[j]
    center <- colMeans(fill)
    scale <- apply(fill, 2L, stats::sd)
    numeric_constant <- vapply(seq_len(ncol(fill)), function(j) v15_numeric_constant(fill[, j]), logical(1))
    support <- t(vapply(seq_len(ncol(raw_support)), function(j) v15_support_stats(raw_support[, j]), 
        numeric(7)))
    keep <- is.finite(med) & is.finite(scale) & scale > 0 & !numeric_constant & !as.logical(support[, 
        "Drop"])
    reason <- ifelse(!is.finite(med), "No_finite_BLUE_training_values", ifelse(numeric_constant, "Numerically_constant_BLUE", 
        ifelse(support[, "N"] == 0, "No_finite_raw_training_values", ifelse(as.logical(support[, "Drop"]), 
            "Insufficient_raw_training_support", "Retained"))))
    list(features = colnames(x), median = med, center = center, scale = scale, keep = keep, n_finite_train = observed, 
        n_train = nrow(x), numeric_constant = numeric_constant, raw_support = support, drop_reason = reason, 
        train_min = vapply(seq_len(ncol(fill)), function(j) if (all(is.na(fill[, j]))) NA_real_ else min(fill[, 
            j]), numeric(1)), train_max = vapply(seq_len(ncol(fill)), function(j) if (all(is.na(fill[, 
            j]))) NA_real_ else max(fill[, j]), numeric(1)))
}

v15_prep_table <- function(p) {
    data.table::data.table(Feature = p$features, Median = p$median, Center = p$center, Scale = p$scale, 
        Retained = p$keep, N_Finite_Train = p$n_finite_train, N_Train = p$n_train, Numerical_Constant = p$numeric_constant, 
        Drop_Reason = p$drop_reason, N_Raw_Training = p$raw_support[, "N"], Raw_Resolution = p$raw_support[, 
            "Resolution"], Raw_Unique = p$raw_support[, "Unique"], Raw_Unique_Pct = p$raw_support[, "UniquePct"], 
        Raw_Frequency_Ratio = p$raw_support[, "FrequencyRatio"], Raw_Dominant_Count = p$raw_support[, 
            "DominantCount"], Train_Min = p$train_min, Train_Max = p$train_max)
}

v15_application_diagnostics <- function(x, prep) {
    scaled <- v15_preprocess_apply(x, prep)
    if (!ncol(scaled)) 
        return(data.table::data.table(Feature = character(), Max_Absolute_Z = numeric(), N_Outside_Training_Range = integer(), 
            N_Applied = integer(), Extreme_Support_Flag = logical()))
    raw <- x[, prep$keep, drop = FALSE]
    lo <- prep$train_min[prep$keep]
    hi <- prep$train_max[prep$keep]
    outside <- sweep(raw, 2L, lo, "<") | sweep(raw, 2L, hi, ">")
    mz <- apply(abs(scaled), 2L, max)
    data.table::data.table(Feature = colnames(scaled), Max_Absolute_Z = mz, N_Outside_Training_Range = colSums(outside, 
        na.rm = TRUE), N_Applied = nrow(x), Extreme_Support_Flag = mz > 100)
}

v15_preprocess_apply <- function(x, prep) {
    x <- as.matrix(x)
    if (!identical(colnames(x), prep$features)) 
        stop("Preprocessing feature names/order differ")
    storage.mode(x) <- "double"
    x[!is.finite(x)] <- NA_real_
    x <- x[, prep$keep, drop = FALSE]
    if (!ncol(x)) 
        return(x)
    for (j in seq_len(ncol(x))) x[is.na(x[, j]), j] <- prep$median[prep$keep][j]
    sweep(sweep(x, 2L, prep$center[prep$keep], "-"), 2L, prep$scale[prep$keep], "/")
}

v15_path <- function(x, y, alpha, lambdas) {
    stopifnot(all(is.finite(x)), all(is.finite(y)))
    if (!ncol(x) || stats::sd(y) == 0) {
        return(list(intercept_only = TRUE, mean_y = mean(y), lambdas = lambdas, feature_names = colnames(x)))
    }
    original_names <- colnames(x)
    if (ncol(x) == 1L) 
        x <- cbind(x, .V15_ZERO_DUMMY = 0)
    fit <- glmnet::glmnet(x = x, y = y, family = "gaussian", alpha = alpha, lambda = lambdas, standardize = FALSE, 
        intercept = TRUE, thresh = 1e-07, maxit = 100000L)
    if (fit$jerr != 0L || length(fit$lambda) != length(lambdas) || max(abs(fit$lambda - lambdas)) > 1e-10) {
        stop("Incomplete/nonconverged glmnet path; no silent interpolation or partial tuning")
    }
    list(intercept_only = FALSE, model = fit, feature_names = original_names)
}

v15_path_predict <- function(fit, x, lambdas) {
    if (fit$intercept_only) 
        return(matrix(fit$mean_y, nrow(x), length(lambdas)))
    if (length(fit$feature_names) == 1L) 
        x <- cbind(x, .V15_ZERO_DUMMY = 0)
    as.matrix(stats::predict(fit$model, newx = x, s = lambdas, type = "response"))
}

v15_path_coef <- function(fit, lambda, prep) {
    cf <- setNames(numeric(length(prep$features)), prep$features)
    intercept <- if (fit$intercept_only) 
        fit$mean_y
    else {
        raw <- as.matrix(stats::coef(fit$model, s = lambda))[, 1]
        cf[fit$feature_names] <- raw[fit$feature_names]
        unname(raw["(Intercept)"])
    }
    abs_cf <- abs(cf)
    data.table::data.table(Feature = c("(Intercept)", names(cf)), Coefficient_Standardized = c(intercept, 
        unname(cf)), Selected = c(NA, unname(abs_cf > 0)), Importance = c(NA, if (max(abs_cf) > 0) 100 * 
        abs_cf/max(abs_cf) else rep(0, length(cf))))
}

v15_inner_tune <- function(x, y, inner_fold, alphas, lambdas, raw_support) {
    stopifnot(nrow(x) == length(y), length(y) == length(inner_fold), all(is.finite(y)))
    folds <- sort(unique(inner_fold))
    if (length(folds) < 2L) 
        stop("At least two inner folds are required")
    sse <- matrix(0, length(alphas), length(lambdas))
    inner_rows <- list()
    prep_rows <- list()
    support_rows <- list()
    for (fold in folds) {
        train <- which(inner_fold != fold)
        valid <- which(inner_fold == fold)
        prep <- v15_preprocess_fit(x[train, , drop = FALSE], raw_support[train, , drop = FALSE])
        xtr <- v15_preprocess_apply(x[train, , drop = FALSE], prep)
        xva <- v15_preprocess_apply(x[valid, , drop = FALSE], prep)
        prep_rows[[length(prep_rows) + 1L]] <- v15_prep_table(prep)[, `:=`(Inner_Fold, fold)]
        support_rows[[length(support_rows) + 1L]] <- v15_application_diagnostics(x[valid, , drop = FALSE], 
            prep)[, `:=`(Inner_Fold, fold)]
        for (a in seq_along(alphas)) {
            fit <- v15_path(xtr, y[train], alphas[a], lambdas)
            p <- v15_path_predict(fit, xva, lambdas)
            if (any(!is.finite(p))) 
                stop("Nonfinite inner predictions")
            err <- colSums((p - y[valid])^2)
            sse[a, ] <- sse[a, ] + err
            inner_rows[[length(inner_rows) + 1L]] <- data.table::data.table(Inner_Fold = fold, Alpha = alphas[a], 
                Lambda = lambdas, N_Validation = length(valid), SSE = err, RMSE = sqrt(err/length(valid)))
        }
    }
    grid <- data.table::rbindlist(lapply(seq_along(alphas), function(a) {
        data.table::data.table(Alpha = alphas[a], Lambda = lambdas, Inner_RMSE = sqrt(sse[a, ]/length(y)))
    }))
    data.table::setorder(grid, Inner_RMSE, -Lambda, Alpha)
    list(best = grid[1L], grid = grid, fold_metrics = data.table::rbindlist(inner_rows), preprocessing = data.table::rbindlist(prep_rows), 
        application_diagnostics = data.table::rbindlist(support_rows))
}

v15_fit_outer_scenario <- function(x_train, y_train, x_test, inner_fold, alphas, lambdas, raw_support_train) {
    tune <- v15_inner_tune(x_train, y_train, inner_fold, alphas, lambdas, raw_support_train)
    prep <- v15_preprocess_fit(x_train, raw_support_train)
    xtr <- v15_preprocess_apply(x_train, prep)
    xte <- v15_preprocess_apply(x_test, prep)
    fit <- v15_path(xtr, y_train, tune$best$Alpha, lambdas)
    prediction <- as.numeric(v15_path_predict(fit, xte, tune$best$Lambda))
    if (any(!is.finite(prediction))) 
        stop("Nonfinite outer predictions")
    list(prediction = prediction, tuning = tune, preprocessing = prep, coefficients = v15_path_coef(fit, 
        tune$best$Lambda, prep), application_diagnostics = v15_application_diagnostics(x_test, prep), 
        solver = list(intercept_only = fit$intercept_only, jerr = if (fit$intercept_only) 0L else fit$model$jerr, 
            lambda_count = if (fit$intercept_only) length(lambdas) else length(fit$model$lambda)))
}

v15_adaptive_choice <- function(tuning_rows, scenario_order) {
    candidates <- tuning_rows[tuning_rows$Scenario != "Median", ]
    if (!nrow(candidates)) 
        return(NA_character_)
    candidates$Scenario_Order <- match(candidates$Scenario, scenario_order)
    candidates <- candidates[order(candidates$Inner_RMSE, candidates$Scenario_Order), ]
    candidates$Scenario[1L]
}

v15_read_inputs <- function(cfg) {
    for (p in c(cfg$yield, cfg$predictors, cfg$raw_support)) if (!file.exists(p)) 
        stop("Input is missing: ", p)
    yc <- c("Env", "Pedigree", "Yield.t.ha.BLUE")
    bc <- c("Env", "Pedigree", "DAP", "Probability", "Vegetation.Index", "Quantile.BLUE")
    y <- data.table::fread(cfg$yield, select = yc, showProgress = FALSE)
    b <- data.table::fread(cfg$predictors, select = bc, showProgress = FALSE)
    y[, `:=`(Yield_Source_Row, .I)]
    invalid_y <- is.na(y$Env) | is.na(y$Pedigree) | !nzchar(trimws(y$Env)) | !nzchar(trimws(y$Pedigree))
    yield_id_exclusions <- data.table::copy(y[invalid_y])
    yield_id_exclusions[, `:=`(Reason, "Missing environment/genotype join identifier")]
    y <- y[!invalid_y]
    y <- y[Env %in% cfg$envs]
    b <- b[Env %in% cfg$envs & Vegetation.Index == cfg$target_vi]
    if (!nrow(y) || !nrow(b)) 
        stop("No rows for requested environments/VI")
    if (anyNA(y[, .(Env, Pedigree)]) || anyNA(b[, .(Env, Pedigree, DAP, Probability)])) {
        stop("Missing join/feature keys in inputs")
    }
    if (anyDuplicated(y[, .(Env, Pedigree)])) 
        stop("Yield keys Env/Pedigree are not unique")
    if (!is.numeric(y$Yield.t.ha.BLUE) || !is.numeric(b$Quantile.BLUE)) 
        stop("Nonnumeric response/predictors")
    if (any(!is.finite(b$DAP)) || any(!is.finite(b$Probability))) 
        stop("Nonfinite date/quantile keys")
    if (any(b$Probability < 0.01 | b$Probability > 1 | abs(b$Probability * 100 - round(b$Probability * 
        100)) > 1e-07)) {
        stop("Quantiles must be integer percentages from 0.01 to 1.00")
    }
    b[, `:=`(ProbPct, as.integer(round(Probability * 100)))]
    if (anyDuplicated(b[, .(Env, Pedigree, DAP, ProbPct, Vegetation.Index)])) {
        stop("Duplicate genotype/date/quantile BLUE keys; refusing implicit aggregation")
    }
    b[!is.finite(Quantile.BLUE), `:=`(Quantile.BLUE, NA_real_)]
    sc <- c("Env", "Pedigree", "DAP", "Probability", "Vegetation.Index", "Quantile.Raw.Mean", "N.Raw.Finite")
    support <- data.table::fread(cfg$raw_support, select = sc, showProgress = FALSE)
    support <- support[Env %in% cfg$envs & Vegetation.Index == cfg$target_vi]
    keys <- c("Env", "Pedigree", "DAP", "Probability", "Vegetation.Index")
    if (anyNA(support[, ..keys]) || anyDuplicated(support[, ..keys])) 
        stop("Invalid raw support keys")
    if (!is.numeric(support$Quantile.Raw.Mean) || anyNA(support$N.Raw.Finite) || any(support$N.Raw.Finite < 
        0 | support$N.Raw.Finite != round(support$N.Raw.Finite))) 
        stop("Invalid raw support values/counts")
    if (any((support$N.Raw.Finite > 0) != is.finite(support$Quantile.Raw.Mean))) 
        stop("Raw support mean/count disagreement")
    support[, `:=`(Support_Row_Present, TRUE)]
    b <- merge(b, support, by = keys, all.x = TRUE, sort = FALSE)
    if (anyNA(b$Support_Row_Present)) 
        stop("Raw support is missing BLUE keys; refusing unverified support")
    if (any(is.finite(b$Quantile.BLUE) & b$N.Raw.Finite == 0)) 
        stop("Finite BLUE has no contributing raw measurements; source lineage is incomplete")
    list(yield = y, blue = b, yield_id_exclusions = yield_id_exclusions)
}

v15_environment_matrix <- function(inputs, env) {
    y <- inputs$yield[Env == env]
    b <- inputs$blue[Env == env]
    if (!nrow(y) || !nrow(b)) 
        stop("Requested environment lacks source rows: ", env)
    finite_y <- y[is.finite(Yield.t.ha.BLUE)]
    ids <- sort(intersect(finite_y$Pedigree, unique(b$Pedigree)))
    excluded <- data.table::data.table(Pedigree = sort(union(y$Pedigree, b$Pedigree)))
    excluded[, `:=`(Reason, ifelse(!Pedigree %in% y$Pedigree, "No yield row", ifelse(!Pedigree %in% finite_y$Pedigree, 
        "Nonfinite yield", ifelse(!Pedigree %in% b$Pedigree, "No VI BLUE row", "Included"))))]
    b <- data.table::copy(b[Pedigree %in% ids])
    b[, `:=`(Feature, paste0("DAP", DAP, "_Q", sprintf("%03d", ProbPct), "_", Vegetation.Index))]
    features <- unique(b[, .(Feature, DAP, ProbPct, Vegetation.Index)])
    data.table::setorder(features, DAP, ProbPct)
    wide <- data.table::dcast(b, Pedigree ~ Feature, value.var = "Quantile.BLUE", fill = NA_real_)
    wide <- wide[match(ids, Pedigree)]
    x <- as.matrix(wide[, features$Feature, with = FALSE])
    rownames(x) <- ids
    raw_wide <- data.table::dcast(b, Pedigree ~ Feature, value.var = "Quantile.Raw.Mean", fill = NA_real_)
    raw_wide <- raw_wide[match(ids, Pedigree)]
    raw_support <- as.matrix(raw_wide[, features$Feature, with = FALSE])
    rownames(raw_support) <- ids
    stopifnot(identical(dimnames(raw_support), dimnames(x)))
    yy <- finite_y$Yield.t.ha.BLUE[match(ids, finite_y$Pedigree)]
    stopifnot(identical(wide$Pedigree, ids), !anyNA(yy), !anyDuplicated(ids))
    list(Env = env, ids = ids, x = x, raw_support = raw_support, y = yy, features = features, exclusions = excluded, 
        n_yield_source = nrow(y), n_blue_source = nrow(b))
}

v15_metrics <- function(obs, pred) {
    stopifnot(length(obs) == length(pred), all(is.finite(obs)), all(is.finite(pred)))
    n <- length(obs)
    pearson <- if (n >= 3L && stats::sd(obs) > 0 && stats::sd(pred) > 0) 
        stats::cor(obs, pred)
    else NA_real_
    sse <- sum((obs - pred)^2)
    sst <- sum((obs - mean(obs))^2)
    list(N = n, Pearson = pearson, RMSE = sqrt(sse/n), MAE = mean(abs(obs - pred)), Bias = mean(pred - 
        obs), R2_SSE = if (sst > 0) 1 - sse/sst else NA_real_)
}

v15_summarize <- function(predictions, output, complete) {
    if (!nrow(predictions)) 
        return(invisible(NULL))
    if (anyDuplicated(predictions[, .(Env, Scenario, Outer_Repeat, Pedigree)])) {
        stop("Duplicate outer prediction per genotype/repeat/scenario")
    }
    data.table::fwrite(predictions, file.path(output, "Outer_Predictions_V15.csv"))
    folds <- predictions[, v15_metrics(Obs, Pred), by = .(Env, Scenario, Outer_Repeat, Outer_Fold)]
    repeats <- predictions[, v15_metrics(Obs, Pred), by = .(Env, Scenario, Outer_Repeat)]
    folds[, `:=`(Run_Complete, complete)]
    repeats[, `:=`(Run_Complete, complete)]
    data.table::fwrite(folds, file.path(output, "Outer_Fold_Metrics_V15.csv"))
    data.table::fwrite(repeats, file.path(output, "Outer_Repeat_Metrics_V15.csv"))
    if (!complete) 
        return(invisible(NULL))
    env_summary <- repeats[, .(N_Genotypes = unique(N), N_Repeats = .N, Mean_Pearson = mean(Pearson), 
        SD_Pearson = stats::sd(Pearson), Mean_RMSE = mean(RMSE), SD_RMSE = stats::sd(RMSE), Mean_MAE = mean(MAE), 
        Mean_Bias = mean(Bias), Mean_R2_SSE = mean(R2_SSE)), by = .(Env, Scenario)]
    median <- repeats[Scenario == "Median", .(Env, Outer_Repeat, Median_Pearson = Pearson, Median_RMSE = RMSE, 
        Median_MAE = MAE)]
    paired <- merge(repeats[Scenario != "Median"], median, by = c("Env", "Outer_Repeat"))
    paired[, `:=`(Delta_Pearson = Pearson - Median_Pearson, Delta_RMSE = RMSE - Median_RMSE, Pct_RMSE_Reduction = 100 * 
        (Median_RMSE - RMSE)/Median_RMSE, Comparison = ifelse(Scenario == "all", "Primary: all versus Median", 
        "Secondary/exploratory: versus Median"))]
    data.table::fwrite(env_summary, file.path(output, "Environment_Scenario_Summary_V15.csv"))
    data.table::fwrite(paired, file.path(output, "Paired_Repeat_Comparisons_V15.csv"))
    avg <- predictions[, .(Obs = Obs[1L], Pred = mean(Pred), N_Repeats = .N), by = .(Env, Scenario, Pedigree)]
    avg_metrics <- avg[, v15_metrics(Obs, Pred), by = .(Env, Scenario)]
    data.table::fwrite(avg_metrics, file.path(output, "Hybrid_Averaged_Ensemble_Metrics_Secondary_V15.csv"))
}

v15_run <- function(cfg) {
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
    if (!requireNamespace("glmnet", quietly = TRUE) || !requireNamespace("data.table", quietly = TRUE)) {
        stop("Installed glmnet and data.table are required; this script does not install packages")
    }
    data.table::setDTthreads(1L)
    stopifnot(cfg$outer_folds >= 2L, cfg$inner_folds >= 2L, cfg$outer_repeats >= 1L, all(c(cfg$outer_folds, 
        cfg$inner_folds, cfg$outer_repeats, cfg$seed) == round(c(cfg$outer_folds, cfg$inner_folds, cfg$outer_repeats, 
        cfg$seed))), cfg$max_outer_splits >= 1L, !anyNA(cfg$envs), !anyDuplicated(cfg$envs), all(cfg$scenarios %in% 
        names(v15_scenarios())), all(c("Median", "all") %in% cfg$scenarios), length(cfg$alphas) > 0, 
        all(cfg$alphas >= 0 & cfg$alphas <= 1), length(cfg$lambdas) >= 2L, all(cfg$lambdas > 0), !anyDuplicated(cfg$alphas), 
        !anyDuplicated(cfg$lambdas))
    cfg$lambdas <- sort(cfg$lambdas, decreasing = TRUE)
    inputs <- v15_read_inputs(cfg)
    dir.create(cfg$output, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(cfg$output, "checkpoints"), showWarnings = FALSE)
    cfg_identity <- cfg
    cfg_identity$max_outer_splits <- NULL
    cfg_identity$output <- NULL
    fingerprint <- v15_hash_object(list(config = cfg_identity, input_md5 = tools::md5sum(c(cfg$yield, 
        cfg$predictors, cfg$raw_support)), source_md5 = tools::md5sum(cfg$script_path), glmnet = as.character(utils::packageVersion("glmnet")), 
        R = R.version.string))
    meta_file <- file.path(cfg$output, "Run_Manifest_V15.rds")
    if (file.exists(meta_file)) {
        if (!identical(readRDS(meta_file)$fingerprint, fingerprint)) {
            stop("Output directory contains a different configuration/input/source; choose a new --out-dir")
        }
    }
    else {
        v15_atomic_rds(list(fingerprint = fingerprint, config = cfg, created = Sys.time(), input_md5 = tools::md5sum(c(cfg$yield, 
            cfg$predictors, cfg$raw_support)), session = utils::sessionInfo()), meta_file)
    }
    capture.output(str(cfg), file = file.path(cfg$output, "Configuration_V15.txt"))
    capture.output(utils::sessionInfo(), file = file.path(cfg$output, "SessionInfo_V15.txt"))
    data.table::fwrite(inputs$yield_id_exclusions, file.path(cfg$output, "Yield_Identifier_Exclusions_V15.csv"))
    writeLines(c("Status: RUNNING -- incomplete until all requested outer splits finish", paste("Mode:", 
        cfg$mode), paste("Fingerprint:", fingerprint)), file.path(cfg$output, "Run_Status_V15.txt"))
    writeLines(c("Scheme: nested within-environment whole-genotype CV.", "Primary comparison: all quantiles versus Median; scenarios fixed before evaluation.", 
        "Outer: balanced random genotype folds independent of outcomes; identical across scenarios.", 
        "Inner: training-only raw-support NZV screening, numerical-constant removal, BLUE median imputation, centering and scaling.", 
        "Tuning: pooled inner out-of-fold RMSE; deterministic ties prefer larger lambda, smaller alpha.", 
        "Adaptive distributional scenario uses inner RMSE only; no best outer-test scenario is selected.", 
        "Report mean repeat-wise outer metrics; repeated-fold SDs are descriptive, not independent CIs.", 
        "First-stage BLUE estimation/shared image-yield errors and historical ExGR selection remain outside this CV."), 
        file.path(cfg$output, "Validation_Scheme_V15.txt"))
    scenario_defs <- v15_scenarios()[cfg$scenarios]
    manifests <- list()
    inner_manifests <- list()
    all_results <- list()
    data_qa <- list()
    exclusions <- list()
    split_count <- 0L
    expected <- length(cfg$envs) * cfg$outer_repeats * cfg$outer_folds
    for (env_no in seq_along(cfg$envs)) {
        env <- cfg$envs[env_no]
        dat <- v15_environment_matrix(inputs, env)
        if (length(dat$ids) - ceiling(length(dat$ids)/cfg$outer_folds) < cfg$inner_folds * 2L) {
            stop("Too few outer-training genotypes for requested inner CV: ", env)
        }
        data_qa[[env]] <- data.table::data.table(Env = env, N_Genotypes = length(dat$ids), N_Source_Yield = dat$n_yield_source, 
            N_Source_BLUE = dat$n_blue_source, N_Features = ncol(dat$x), N_Missing_X = sum(!is.finite(dat$x)))
        exclusions[[env]] <- data.table::copy(dat$exclusions)[, `:=`(Env, env)]
        for (repeat_id in seq_len(cfg$outer_repeats)) {
            outer_seed <- cfg$seed + env_no * 100000L + repeat_id * 1000L
            fold <- v15_fold_ids(dat$ids, cfg$outer_folds, outer_seed)
            manifests[[length(manifests) + 1L]] <- data.table::data.table(Env = env, Outer_Repeat = repeat_id, 
                Pedigree = dat$ids, Outer_Test_Fold = fold, Seed = outer_seed)
            for (outer in seq_len(cfg$outer_folds)) {
                split_count <- split_count + 1L
                train <- which(fold != outer)
                test <- which(fold == outer)
                inner_seed <- outer_seed + outer
                inner <- v15_fold_ids(dat$ids[train], cfg$inner_folds, inner_seed)
                inner_manifests[[length(inner_manifests) + 1L]] <- data.table::data.table(Env = env, 
                  Outer_Repeat = repeat_id, Outer_Fold = outer, Pedigree = dat$ids[train], Inner_Validation_Fold = inner, 
                  Seed = inner_seed)
                if (split_count > cfg$max_outer_splits) 
                  next
                result_rows <- list()
                tune_rows <- list()
                for (scenario in cfg$scenarios) {
                  key <- sprintf("%s_R%02d_F%02d_%s", env, repeat_id, outer, scenario)
                  path <- file.path(cfg$output, "checkpoints", paste0(key, ".rds"))
                  if (file.exists(path)) {
                    res <- readRDS(path)
                    if (!identical(res$fingerprint, fingerprint)) 
                      stop("Checkpoint fingerprint mismatch: ", key)
                  }
                  else {
                    message("Nested CV ", key, " n_train=", length(train), " n_test=", length(test))
                    keep <- dat$features$ProbPct %in% scenario_defs[[scenario]]
                    if (!any(keep)) 
                      stop("Scenario has no available candidate features: ", key)
                    fit <- v15_fit_outer_scenario(dat$x[train, keep, drop = FALSE], dat$y[train], dat$x[test, 
                      keep, drop = FALSE], inner, cfg$alphas, cfg$lambdas, dat$raw_support[train, keep, 
                      drop = FALSE])
                    pred <- data.table::data.table(Env = env, Scenario = scenario, Outer_Repeat = repeat_id, 
                      Outer_Fold = outer, Pedigree = dat$ids[test], Obs = dat$y[test], Pred = fit$prediction, 
                      Alpha = fit$tuning$best$Alpha, Lambda = fit$tuning$best$Lambda, Inner_RMSE = fit$tuning$best$Inner_RMSE, 
                      Selected_Scenario = scenario, N_Train = length(train), N_Candidate_Features = sum(keep), 
                      N_Retained_Features = sum(fit$preprocessing$keep), Validation_Scheme = "Nested genotype CV V15")
                    res <- list(fingerprint = fingerprint, predictions = pred, fit = fit, train_ids = dat$ids[train], 
                      test_ids = dat$ids[test], inner_fold = inner, Env = env, Scenario = scenario, Outer_Repeat = repeat_id, 
                      Outer_Fold = outer)
                    v15_atomic_rds(res, path)
                  }
                  result_rows[[scenario]] <- res$predictions
                  tune_rows[[scenario]] <- data.table::data.table(Scenario = scenario, Inner_RMSE = res$fit$tuning$best$Inner_RMSE)
                  all_results[[key]] <- res
                }
                if (cfg$adaptive) {
                  winner <- v15_adaptive_choice(as.data.frame(data.table::rbindlist(tune_rows)), cfg$scenarios)
                  if (!is.na(winner)) {
                    adaptive <- data.table::copy(result_rows[[winner]])
                    adaptive[, `:=`(Scenario = "Adaptive_Inner_Distributional", Selected_Scenario = winner)]
                    all_results[[paste(env, repeat_id, outer, "adaptive", sep = "_")]] <- list(predictions = adaptive)
                  }
                }
            }
        }
    }
    data.table::fwrite(data.table::rbindlist(manifests), file.path(cfg$output, "Outer_Fold_Manifest_V15.csv"))
    data.table::fwrite(data.table::rbindlist(inner_manifests), file.path(cfg$output, "Inner_Fold_Manifest_V15.csv"))
    data.table::fwrite(data.table::rbindlist(data_qa), file.path(cfg$output, "Input_QA_V15.csv"))
    data.table::fwrite(data.table::rbindlist(exclusions), file.path(cfg$output, "Cohort_Inclusion_V15.csv"))
    predictions <- data.table::rbindlist(lapply(all_results, `[[`, "predictions"))
    complete <- cfg$max_outer_splits >= expected
    v15_summarize(predictions, cfg$output, complete)
    fits <- all_results[vapply(all_results, function(x) !is.null(x$fit), logical(1))]
    bind_fit_table <- function(getter) data.table::rbindlist(lapply(fits, function(x) {
        d <- data.table::as.data.table(getter(x$fit))
        d[, `:=`(Env = x$Env, Scenario = x$Scenario, Outer_Repeat = x$Outer_Repeat, Outer_Fold = x$Outer_Fold)]
        d
    }), fill = TRUE)
    data.table::fwrite(bind_fit_table(function(x) x$coefficients), file.path(cfg$output, "Outer_Coefficients_V15.csv"))
    data.table::fwrite(bind_fit_table(function(x) x$tuning$grid), file.path(cfg$output, "Inner_Tuning_Grid_V15.csv"))
    data.table::fwrite(bind_fit_table(function(x) x$tuning$fold_metrics), file.path(cfg$output, "Inner_Fold_Tuning_Metrics_V15.csv"))
    data.table::fwrite(bind_fit_table(function(x) x$tuning$preprocessing), file.path(cfg$output, "Inner_Preprocessing_V15.csv"))
    data.table::fwrite(bind_fit_table(function(x) {
        v15_prep_table(x$preprocessing)
    }), file.path(cfg$output, "Outer_Preprocessing_V15.csv"))
    support_outer <- bind_fit_table(function(x) x$application_diagnostics)
    support_inner <- bind_fit_table(function(x) x$tuning$application_diagnostics)
    data.table::fwrite(support_outer, file.path(cfg$output, "Outer_Input_Support_V15.csv"))
    data.table::fwrite(support_inner, file.path(cfg$output, "Inner_Input_Support_V15.csv"))
    flags <- sum(support_outer$Extreme_Support_Flag) + sum(support_inner$Extreme_Support_Flag)
    writeLines(c(paste("Extreme support flags (absolute standardized predictor >100):", flags), if (flags) "REVIEW REQUIRED: computational completion does not establish supported predictions." else "No extreme support flags under the declared diagnostic; ordinary extrapolation may remain."), 
        file.path(cfg$output, "Input_Support_Review_V15.txt"))
    writeLines(c(paste("Status:", if (complete) "COMPLETE" else "PARTIAL -- not a reportable CV result"), 
        paste("Mode:", cfg$mode), paste("Expected outer splits:", expected), paste("Completed outer splits:", 
            min(expected, cfg$max_outer_splits)), paste("Completed scenario fits:", length(fits)), paste("Fingerprint:", 
            fingerprint)), file.path(cfg$output, "Run_Status_V15.txt"))
    message("V15 ", if (complete) 
        "complete"
    else "partial", ": ", cfg$output)
    invisible(list(config = cfg, complete = complete, predictions = predictions, fingerprint = fingerprint))
}

