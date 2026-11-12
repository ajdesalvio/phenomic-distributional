parse_quantile <- function(x) {
    pattern <- "^FPC[0-9]+\\.[^.]+\\.(0(?:\\.[0-9]{1,2})?|1(?:\\.0{1,2})?)$"
    if (any(!grepl(pattern, x))) 
        stop("Malformed predictor: ", paste(x[!grepl(pattern, x)], collapse = ", "))
    q <- as.integer(round(as.numeric(sub(pattern, "\\1", x)) * 100))
    if (any(q < 1L | q > 100L)) 
        stop("Quantile outside 1:100")
    q
}

make_wide <- function(he) {
    meta <- manifest$partitions[manifest$partitions$heldout == he, ]
    p <- file.path(cache, meta$file)
    stopifnot(file.info(p)$size == meta$bytes)
    message("Reading saved projections for ", he)
    fpc <- fread(p)
    stopifnot(nrow(fpc) == meta$rows, all(fpc$Heldout.Env == he))
    fpc <- fpc[Pedigree.Env %in% yield$Pedigree.Env]
    fv <- grep("^FPC[0-9]+_FVE$", names(fpc), value = TRUE)
    fv <- fv[order(as.integer(sub("^FPC([0-9]+)_FVE$", "\\1", fv)))]
    vis <- sort(unique(fpc$Vegetation.Index))
    qs <- sort(unique(fpc$Quantile.Level))
    stopifnot(identical(as.integer(round(qs * 100)), 1:100))
    keys <- c("Pedigree.Env", "Pedigree", "Env", "Heldout.Env")
    pieces <- list()
    feature_info <- list()
    for (vi in vis) for (q in qs) {
        d <- fpc[Vegetation.Index == vi & Quantile.Level == q]
        setorderv(d, c("Env", "Pedigree"))
        stopifnot(nrow(d) > 0L, !anyDuplicated(d$Pedigree.Env))
        if (any(vapply(d[, ..fv], function(z) uniqueN(z, na.rm = FALSE) > 1L, logical(1)))) 
            stop("FVE metadata varies within heldout/VI/quantile")
        fve <- as.numeric(d[1, ..fv])
        fve[is.na(fve)] <- 0
        reached <- which(cumsum(fve) >= 0.95)
        k <- if (length(reached)) 
            reached[1]
        else length(fve)
        pc <- paste0("FPC", seq_len(k))
        stopifnot(all(pc %in% names(d)))
        keep <- c(keys, pc)
        piece <- d[, ..keep]
        new <- paste(pc, vi, format(q, trim = TRUE, scientific = FALSE), sep = ".")
        setnames(piece, pc, new)
        pieces[[length(pieces) + 1L]] <- piece
        feature_info[[length(feature_info) + 1L]] <- data.table(Heldout = he, Feature = new, VI = vi, 
            Quantile = parse_quantile(new), FPC = seq_len(k), FVE = fve[seq_len(k)], Cumulative_FVE = cumsum(fve)[seq_len(k)])
    }
    common <- Reduce(intersect, lapply(pieces, function(d) d$Pedigree.Env))
    yy <- yield[Pedigree.Env %in% common]
    pos <- match(yy$Pedigree.Env, pieces[[1]]$Pedigree.Env)
    wide <- cbind(yy, pieces[[1]][pos, setdiff(keys, "Pedigree.Env"), with = FALSE])
    for (piece in pieces) {
        ix <- match(yy$Pedigree.Env, piece$Pedigree.Env)
        stopifnot(identical(piece$Env[ix], wide$Env), identical(piece$Pedigree[ix], wide$Pedigree))
        wide <- cbind(wide, piece[ix, setdiff(names(piece), keys), with = FALSE])
    }
    fwrite(rbindlist(feature_info), file.path(output, "features", paste0(he, "_retained_FPCs.csv")))
    message("Wide predictor matrix: ", nrow(wide), " records; ", ncol(wide) - 6L, " predictors")
    wide
}

scale_env <- function(X, env) {
    for (e in unique(env)) {
        ix <- which(env == e)
        for (j in seq_along(X)) {
            z <- X[[j]][ix]
            s <- sd(z, na.rm = TRUE)
            m <- mean(z, na.rm = TRUE)
            X[[j]][ix] <- if (is.na(s) || s == 0) 
                0
            else (z - m)/s
        }
    }
    X
}

fit_one <- function(wide, he, scenario, scaling) {
    start <- Sys.time()
    predictors <- grep("^FPC[0-9]+\\.", names(wide), value = TRUE)
    predictors <- predictors[!parse_quantile(predictors) %in% filters[[scenario]]]
    stopifnot(length(unique(parse_quantile(predictors))) == 100L - length(filters[[scenario]]))
    tr <- wide[Env != he]
    te <- wide[Env == he]
    Xtr <- as.data.frame(tr[, ..predictors])
    Xte <- as.data.frame(te[, ..predictors])
    if (scaling == "WITHIN") {
        Xtr <- scale_env(Xtr, tr$Env)
        Xte <- scale_env(Xte, te$Env)
    }
    bad <- vapply(Xtr, function(z) any(!is.finite(z)), logical(1))
    Xtr <- Xtr[, !bad, drop = FALSE]
    nzv <- caret::nearZeroVar(Xtr)
    if (length(nzv)) 
        Xtr <- Xtr[, -nzv, drop = FALSE]
    keep <- complete.cases(Xtr) & is.finite(tr$Yield.t.ha.BLUE)
    Xtr <- Xtr[keep, , drop = FALSE]
    tr <- tr[keep]
    Xte <- Xte[, names(Xtr), drop = FALSE]
    if (any(!is.finite(as.matrix(Xte)))) 
        stop("Nonfinite heldout predictors")
    folds_out <- split(seq_len(nrow(tr)), tr$Env)
    folds_in <- lapply(folds_out, function(i) setdiff(seq_len(nrow(tr)), i))
    stopifnot(!he %in% names(folds_out))
    ctrl <- trainControl(method = "cv", index = folds_in, indexOut = folds_out, verboseIter = FALSE, 
        allowParallel = workers > 1L, returnResamp = "final", savePredictions = "final")
    pre <- if (scaling == "TRAIN") 
        c("zv", "center", "scale")
    else "zv"
    message("FIT ", he, " / ", scenario, " / ", scaling, " (", nrow(tr), " x ", ncol(Xtr), ")")
    warnings <- character()
    set.seed(123)
    fit <- withCallingHandlers(caret::train(x = Xtr, y = tr$Yield.t.ha.BLUE, method = "glmnet", trControl = ctrl, 
        preProcess = pre, tuneLength = 25), warning = function(w) {
        warnings <<- unique(c(warnings, conditionMessage(w)))
        invokeRestart("muffleWarning")
    })
    if (all(is.na(fit$results$RMSE))) 
        stop("All tuning RMSE values missing")
    chosen <- fit$bestTune$lambda
    path <- fit$finalModel$lambda
    jerr <- fit$finalModel$jerr
    if (length(chosen) != 1L || !is.finite(chosen) || chosen <= 0 || !length(path) || any(!is.finite(path)) || 
        any(path <= 0)) 
        stop("Invalid selected lambda/fitted path")
    if (length(jerr) != 1L || !is.finite(jerr) || jerr > 0) 
        stop("Fatal final glmnet path error")
    if (nrow(fit$resample) != length(folds_out) || any(!is.finite(fit$resample$RMSE))) 
        stop("Selected model lacks finite RMSE in an expected inner environment")
    if (nrow(fit$pred) != nrow(tr) || anyDuplicated(fit$pred$rowIndex) || !setequal(fit$pred$rowIndex, 
        seq_len(nrow(tr))) || any(!is.finite(fit$pred$pred))) 
        stop("Selected inner predictions do not cover every training record exactly once")
    clamp <- if (chosen < min(path)) 
        "below_path_min"
    else if (chosen > max(path)) 
        "above_path_max"
    else "within_path"
    effective <- max(min(path), min(max(path), chosen))
    pred <- as.numeric(predict(fit, newdata = Xte))
    obs <- te$Yield.t.ha.BLUE
    if (length(pred) != nrow(te) || any(!is.finite(pred))) 
        stop("Invalid heldout predictions")
    cf <- as.matrix(coef(fit$finalModel, s = fit$bestTune$lambda))[-1, 1]
    retained <- is.finite(pred) & is.finite(obs)
    pp <- pred[retained]
    oo <- obs[retained]
    cor_safe <- function(method) if (length(pp) >= 3 && sd(pp) > 0 && sd(oo) > 0) 
        cor(pp, oo, method = method)
    else NA_real_
    metrics <- data.table(Heldout = he, Scenario = scenario, Scaling = scaling, RustScore = "Yield.t.ha.BLUE", 
        n_test = sum(retained), n_train = nrow(tr), Pearson = cor_safe("pearson"), Spearman = cor_safe("spearman"), 
        RMSE = sqrt(mean((pp - oo)^2)), MAE = mean(abs(pp - oo)), SD_Pred = sd(pp), SD_Obs = sd(oo), 
        NZ_Coefs = sum(abs(cf) > 0), Alpha = fit$bestTune$alpha, Lambda = fit$bestTune$lambda, Path_Lambda_Min = min(path), 
        Path_Lambda_Max = max(path), Path_Lambda_Count = length(path), Lambda_Evaluation = clamp, Effective_Lambda = effective, 
        Final_Jerr = jerr, Final_NPasses = fit$finalModel$npasses, Distinct_Warning_Count = length(warnings), 
        Grid_Alpha_Count = uniqueN(fit$results$alpha), Grid_Lambda_Count = uniqueN(fit$results$lambda), 
        Grid_Candidates = nrow(fit$results), Missing_Tuning_RMSE = sum(!is.finite(fit$results$RMSE)), 
        Alpha_Grid_Boundary = fit$bestTune$alpha %in% range(fit$results$alpha), Lambda_Grid_Boundary = chosen %in% 
            range(fit$results$lambda), Predictors = ncol(Xtr), Quantiles = length(unique(parse_quantile(names(Xtr)))), 
        Elapsed_Seconds = as.numeric(difftime(Sys.time(), start, units = "secs")))
    predictions <- data.table(Heldout = he, Scenario = scenario, Scaling = scaling, Pedigree.Env = te$Pedigree.Env, 
        Pedigree = te$Pedigree, Env = te$Env, Obs = obs, Pred = pred, Residual = obs - pred)
    coefficients <- data.table(Heldout = he, Scenario = scenario, Scaling = scaling, Feature = names(cf), 
        Quantile = parse_quantile(names(cf)), Coefficient = cf, Selected = abs(cf) > 0, Importance = if (max(abs(cf)) > 
            0) 
            100 * abs(cf)/max(abs(cf))
        else 0)
    tuning <- cbind(data.table(Heldout = he, Scenario = scenario, Scaling = scaling), as.data.table(fit$results))
    folds <- data.table(Heldout = he, Scenario = scenario, Scaling = scaling, Pedigree.Env = tr$Pedigree.Env, 
        Inner_Validation_Env = tr$Env)
    fit$trainingData <- NULL
    list(config_hash = config_hash, metrics = metrics, predictions = predictions, coefficients = coefficients, 
        tuning = tuning, folds = folds, fit = fit, warnings = warnings, removed_nonfinite = predictors[bad], 
        removed_nearzero = if (length(nzv)) predictors[!bad][nzv] else character())
}

