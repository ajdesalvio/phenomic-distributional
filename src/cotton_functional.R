v7_formula <- function() response ~ Pedigree + Flight + Block + (1 | Plant)

v7_hash <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)

v7_signature <- function(provenance) {
    digest::digest(as.character(jsonlite::toJSON(provenance, auto_unbox = TRUE, digits = NA, null = "null")), 
        serialize = FALSE, algo = "sha256")
}

v7_save <- function(object, path) {
    if (file.exists(path)) 
        stop("Refusing to overwrite: ", path)
    tmp <- paste0(path, ".partial_", Sys.getpid())
    if (file.exists(tmp)) 
        stop("Inspect existing partial file: ", tmp)
    saveRDS(object, tmp)
    if (!file.rename(tmp, path)) 
        stop("Atomic checkpoint rename failed: ", tmp)
}

v7_csv <- function(x, path) {
    if (file.exists(path)) {
        character_columns <- names(x)[vapply(x, is.character, logical(1))]
        existing <- as.data.frame(data.table::fread(path, colClasses = list(character = character_columns)))
        if (!isTRUE(all.equal(existing, x, check.attributes = FALSE, tolerance = 1e-12))) 
            stop("Existing CSV differs: ", path)
    }
    else {
        tmp <- paste0(path, ".partial_", Sys.getpid())
        if (file.exists(tmp)) 
            stop("Partial CSV exists: ", tmp)
        data.table::fwrite(x, tmp)
        if (!file.rename(tmp, path)) 
            stop("CSV rename failed.")
    }
}

v7_prepare <- function(cfg) {
    if (grepl("[.]rdata$", cfg$input, ignore.case = TRUE)) {
        env <- new.env(parent = emptyenv())
        load(cfg$input, envir = env)
        if (!exists("data2", env, inherits = FALSE)) 
            stop("Workspace lacks data2.")
        src <- env$data2
    }
    else if (grepl("[.]rds$", cfg$input, ignore.case = TRUE)) 
        src <- readRDS(cfg$input)
    else stop("Input must be saved RData/data2 or an RDS data.frame, never an executable script.")
    required <- c("Test", "Tier", "Row", "Range", "Pedigree", "Pltg_ID", "DAT")
    if (!is.data.frame(src) || !all(required %in% names(src))) 
        stop("Field/date metadata are missing.")
    prefix <- paste0("^", cfg$vi, "_Quant_")
    cols <- grep(prefix, names(src), value = TRUE)
    prob <- suppressWarnings(as.numeric(sub(prefix, "", cols)))
    if (length(prob) != 100L || anyNA(prob) || anyDuplicated(prob) || max(abs(sort(prob) - seq(0.01, 
        1, 0.01))) > 1e-10) 
        stop("Invalid VI probability grid.")
    cols <- cols[order(prob)]
    prob <- sort(prob)
    d <- src[, required, drop = FALSE]
    Y <- as.matrix(src[, cols, drop = FALSE])
    if (anyNA(d) || !is.numeric(Y)) 
        stop("Incomplete metadata/non-numeric quantiles.")
    keep <- !d$DAT %in% cfg$excluded_dates & !d$Pltg_ID %in% cfg$excluded_plants
    d <- d[keep, , drop = FALSE]
    Y <- Y[keep, , drop = FALSE]
    if (any(!is.finite(Y)) || anyDuplicated(d[, c("Pltg_ID", "DAT")])) 
        stop("Nonfinite quantiles or duplicate plant/date; no automatic deletion.")
    if (any(Y[, -1, drop = FALSE] - Y[, -ncol(Y), drop = FALSE] < -1e-08)) 
        stop("Quantile ordering failed.")
    map <- as.data.frame(data.table::fread(cfg$fieldmap))
    if (!all(setdiff(required, "DAT") %in% names(map)) || anyDuplicated(map$Pltg_ID)) 
        stop("Invalid field map.")
    expected <- as.character(map$Pltg_ID[!map$Pltg_ID %in% cfg$excluded_plants])
    if (!setequal(unique(as.character(d$Pltg_ID)), expected)) 
        stop("Retained plant IDs differ from field map.")
    idx <- match(d$Pltg_ID, map$Pltg_ID)
    for (nm in setdiff(required, "DAT")) if (!identical(as.character(d[[nm]]), as.character(map[[nm]][idx]))) 
        stop("Field-map mismatch: ", nm)
    counts <- c(observations = nrow(d), plants = length(unique(d$Pltg_ID)), dates = length(unique(d$DAT)), 
        genotypes = length(unique(d$Pedigree)), blocks = length(unique(interaction(d$Test, d$Tier, drop = TRUE))))
    if (!identical(as.integer(counts), c(47807L, 479L, 101L, 16L, 10L))) 
        stop("The audited retained cohort changed.")
    cohort_hash <- digest::digest(list(metadata = d, probability = prob, Y = Y), algo = "sha256")
    metadata_hash <- digest::digest(d, algo = "sha256")
    full_dates <- sort(unique(d$DAT))
    if (cfg$mode == "smoke") {
        keep <- d$DAT %in% cfg$smoke_dates
        d <- d[keep, , drop = FALSE]
        Y <- Y[keep, cfg$smoke_indices, drop = FALSE]
        prob <- prob[cfg$smoke_indices]
    }
    d$Pedigree <- stats::relevel(factor(d$Pedigree), ref = "5016")
    d$Flight <- factor(d$DAT, levels = sort(unique(d$DAT)))
    d$Block <- interaction(d$Test, d$Tier, drop = TRUE, sep = "_")
    d$Plant <- factor(d$Pltg_ID)
    contrasts(d$Pedigree) <- stats::contr.treatment(levels(d$Pedigree), base = 1L)
    contrasts(d$Flight) <- stats::contr.sum(nlevels(d$Flight))
    contrasts(d$Block) <- stats::contr.sum(nlevels(d$Block))
    plant <- unique(d[, c("Plant", "Pedigree", "Block", "Row", "Range")])
    if (anyDuplicated(plant$Plant) || any(table(d$Plant) < 2L)) 
        stop("Plant metadata/repeated-observation check failed.")
    X <- stats::model.matrix(~Pedigree + Flight + Block, d)
    if (qr(X)$rank != ncol(X)) 
        stop("Fixed effects are rank deficient.")
    if (max(abs(diff(prob) - diff(prob)[1])) > 1e-10) 
        stop("Use an evenly spaced functional grid.")
    list(data = d, Y = Y, probability = prob, clusters = levels(d$Plant), summary = list(full_counts = counts, 
        full_dates = full_dates, full_missing_plant_dates = 572L, full_cohort_sha256 = cohort_hash, full_metadata_sha256 = metadata_hash, 
        observations = nrow(d), plants = nlevels(d$Plant), dates = nlevels(d$Flight), blocks = nlevels(d$Block), 
        fixed_columns = colnames(X), fixed_rank = qr(X)$rank, probability = prob, plant_ids = levels(d$Plant), 
        plants_per_block_genotype = table(plant$Block, plant$Pedigree)))
}

v7_components <- function(fit, cluster, targets, need_df = FALSE) {
    V <- clubSandwich::vcovCR(fit, cluster = cluster, type = "CR2")
    groups <- attr(V, "cluster")
    ids <- levels(groups)
    if (!identical(ids, levels(cluster))) 
        stop("CR2 cluster ordering changed.")
    X <- stats::model.matrix(fit)
    beta <- lme4::fixef(fit)
    e <- lme4::getME(fit, "y") - as.numeric(X %*% beta)
    rows <- split(seq_len(nrow(X)), groups)
    A <- attr(V, "adjustments")
    XW <- attr(V, "est_mats")
    B <- attr(V, "bread")/attr(V, "v_scale")
    selected <- match(targets, names(beta))
    if (anyNA(selected)) 
        stop("Genotype contrasts absent.")
    R <- length(targets)
    J <- length(ids)
    p <- ncol(X)
    Psi <- matrix(NA_real_, p, J, dimnames = list(names(beta), ids))
    C <- matrix(NA_real_, R, J, dimnames = list(targets, ids))
    K <- matrix(NA_real_, R * J, p, dimnames = list(NULL, names(beta)))
    adjustment_max <- numeric(J)
    for (i in seq_len(J)) {
        Xi <- X[rows[[i]], , drop = FALSE]
        ei <- e[rows[[i]]]
        fi <- as.numeric(A[[i]] %*% ei)
        Psi[, i] <- as.numeric(B %*% (XW[[i]] %*% fi))
        D <- B[selected, , drop = FALSE] %*% XW[[i]] %*% A[[i]]
        C[, i] <- as.numeric(D %*% fi)
        K[((i - 1L) * R + 1L):(i * R), ] <- D %*% Xi
        adjustment_max[i] <- max(abs(A[[i]]))
    }
    error <- max(abs(tcrossprod(Psi) - as.matrix(V)))
    scale <- max(abs(as.matrix(V)))
    if (!is.finite(error) || error > 1e-08 * max(scale, .Machine$double.eps)) 
        stop("CR2 contribution scaling check failed.")
    raw_df <- if (need_df) 
        as.data.frame(clubSandwich::coef_test(fit, vcov = V, coefs = targets, test = "Satterthwaite"))
    else NULL
    contribution_effective_count <- rowSums(Psi[selected, , drop = FALSE]^2)^2/rowSums(Psi[selected, 
        , drop = FALSE]^4)
    list(Psi = Psi, C = C, K = K, targets = targets, clusters = ids, coefficient_names = names(beta), 
        raw_cr2_se = sqrt(diag(V)), raw_contrast_df = raw_df, Gram_max_absolute_error = error, Gram_relative_error = error/max(scale, 
            .Machine$double.eps), contribution_effective_count = contribution_effective_count, maximum_absolute_adjustment_element = max(adjustment_max))
}

v7_fit <- function(d, y, probability, cfg) {
    d$response <- y
    warnings <- character()
    fit_timing <- system.time(fit <- withCallingHandlers(lme4::lmer(v7_formula(), data = d, REML = TRUE, 
        control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = cfg$maxfun), check.rankX = "stop.deficient")), 
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        }))
    conv <- fit@optinfo$conv
    singular <- lme4::isSingular(fit)
    bad_messages <- unlist(conv$lme4$messages)
    bad_messages <- bad_messages[!grepl("singular", bad_messages, ignore.case = TRUE)]
    bad_warnings <- warnings[!grepl("singular", warnings, ignore.case = TRUE)]
    ok <- (is.null(conv$opt) || all(conv$opt == 0)) && length(bad_messages) == 0L && length(bad_warnings) == 
        0L && all(is.finite(lme4::fixef(fit))) && is.finite(stats::sigma(fit)) && stats::sigma(fit) > 
        0 && stats::nobs(fit) == nrow(d)
    if (!ok) 
        return(list(ok = FALSE, probability = probability, warnings = warnings, convergence = conv, singular = singular, 
            fit_timing = fit_timing))
    targets <- paste0("Pedigree", 5001:5015)
    covariance_timing <- system.time(comp <- v7_components(fit, d$Plant, targets, need_df = any(abs(probability - 
        cfg$df_probabilities) < 1e-10)))
    list(ok = TRUE, probability = probability, beta = lme4::fixef(fit), components = comp, variance_components = as.data.frame(lme4::VarCorr(fit)), 
        warnings = warnings, convergence = conv, singular = singular, nobs = stats::nobs(fit), fit_timing = fit_timing, 
        covariance_timing = covariance_timing)
}

v7_point <- function(prepared, cfg, signature) {
    directory <- file.path(cfg$out, "point_checkpoints")
    dir.create(directory, showWarnings = FALSE)
    records <- vector("list", length(prepared$probability))
    for (j in seq_along(records)) {
        path <- file.path(directory, sprintf("quantile_%03d.rds", j))
        if (file.exists(path)) {
            r <- readRDS(path)
            if (!identical(r$signature, signature) || !identical(r$probability, prepared$probability[j])) 
                stop("Point checkpoint mismatch.")
        }
        else {
            cat("V7", cfg$vi, "point", j, "/", length(records), "q=", prepared$probability[j], "\n")
            flush.console()
            r <- v7_fit(prepared$data, prepared$Y[, j], prepared$probability[j], cfg)
            r$signature <- signature
            v7_save(r, path)
        }
        if (!isTRUE(r$ok)) 
            stop("Point diagnostic failure; preserve and inspect ", path)
        if (!identical(r$components$clusters, prepared$clusters)) 
            stop("Plant order differs across quantiles.")
        records[[j]] <- r
    }
    records
}

v7_smooth <- function(records, probability, cfg) {
    raw <- do.call(cbind, lapply(records, `[[`, "beta"))
    L <- ncol(raw)
    argvals <- seq_len(L)
    if (!all(vapply(records, function(r) identical(names(r$beta), rownames(raw)), logical(1)))) 
        stop("Coefficient order differs.")
    k <- min(cfg$max_basis, floor(L/2) + 1L)
    sm <- mgcv::smoothCon(mgcv::s(argvals, bs = "tp", k = k), data = data.frame(argvals = argvals), absorb.cons = TRUE)[[1]]
    B <- sm$X
    S <- sm$S[[1]]
    beta <- raw
    lambda <- edf <- numeric(nrow(raw))
    target_names <- paste0("Pedigree", 5001:5015)
    smoothers <- influences <- covariances <- vector("list", 15L)
    names(smoothers) <- names(influences) <- names(covariances) <- target_names
    smoothing_error <- numeric(nrow(raw))
    for (r in seq_len(nrow(raw))) {
        values <- raw[r, ]
        g <- mgcv::gam(values ~ s(argvals, bs = "tp", k = k), method = cfg$smooth_method)
        M <- B %*% solve(crossprod(B) + g$sp * S, t(B)) + matrix(1/L, L, L)
        beta[r, ] <- as.numeric(g$fitted.values)
        lambda[r] <- g$sp
        edf[r] <- sum(g$edf)
        smoothing_error[r] <- max(abs(as.numeric(M %*% values) - beta[r, ]))
        if (smoothing_error[r] > 1e-08 * max(1, max(abs(beta[r, ])))) 
            stop("Fitted smoother propagation mismatch.")
        if (rownames(raw)[r] %in% target_names) {
            nm <- rownames(raw)[r]
            E <- do.call(rbind, lapply(records, function(z) z$components$Psi[nm, ]))
            T <- M %*% E
            smoothers[[nm]] <- M
            influences[[nm]] <- T
            covariances[[nm]] <- tcrossprod(T)
        }
    }
    colnames(raw) <- colnames(beta) <- sprintf("%.2f", probability)
    se <- do.call(rbind, lapply(influences, function(x) sqrt(rowSums(x^2))))
    if (any(!is.finite(se) | se <= 0)) 
        stop("Degenerate smoothed contrast standard errors.")
    list(beta_raw = raw, beta_smooth = beta, probability = probability, targets = target_names, lambda = setNames(lambda, 
        rownames(raw)), edf = setNames(edf, rownames(raw)), smoothers = smoothers, influences = influences, 
        covariances = covariances, se = se, maximum_smoother_reconstruction_error = max(smoothing_error))
}

v7_wild_batch <- function(records, point, weights) {
    L <- length(records)
    R <- length(point$targets)
    J <- nrow(weights)
    B <- ncol(weights)
    raw_delta <- array(NA_real_, c(L, R, B))
    contributions <- array(NA_real_, c(L, R, J, B))
    for (q in seq_len(L)) {
        z <- records[[q]]$components
        delta <- z$Psi %*% weights
        raw_delta[q, , ] <- delta[point$targets, , drop = FALSE]
        c_weights <- matrix(as.vector(z$C), nrow = R * J, ncol = B) * weights[rep(seq_len(J), each = R), 
            , drop = FALSE]
        u <- c_weights - z$K %*% delta
        contributions[q, , , ] <- array(u, c(R, J, B))
    }
    t_abs <- array(NA_real_, c(R, L, B), dimnames = list(point$targets, NULL, NULL))
    min_se <- Inf
    for (r in seq_len(R)) {
        numerator <- point$smoothers[[r]] %*% matrix(raw_delta[, r, ], nrow = L, ncol = B)
        u <- point$smoothers[[r]] %*% matrix(contributions[, r, , ], nrow = L, ncol = J * B)
        variance <- t(rowsum(t(u^2), group = rep(seq_len(B), each = J), reorder = FALSE))
        if (any(!is.finite(variance) | variance <= 0)) 
            stop("Degenerate pseudo-dataset studentizing variance; no draw discarded.")
        se <- sqrt(variance)
        min_se <- min(min_se, se)
        t_abs[r, , ] <- abs(numerator/se)
    }
    if (any(!is.finite(t_abs))) 
        stop("Nonfinite bootstrap-t draw.")
    list(t_abs = t_abs, joint_max = apply(t_abs, 3, max), per_curve_max = apply(t_abs, c(1, 3), max), 
        minimum_studentizing_se = min_se)
}

v7_bootstrap <- function(records, point, prepared, cfg, signature) {
    directory <- file.path(cfg$out, "wild_checkpoints")
    dir.create(directory, showWarnings = FALSE)
    weights_path <- file.path(cfg$out, "wild_weights_V7.rds")
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(cfg$seed)
    weights <- matrix(sample(c(-1, 1), length(prepared$clusters) * cfg$draws, replace = TRUE), nrow = length(prepared$clusters), 
        dimnames = list(prepared$clusters, NULL))
    weight_hash <- digest::digest(weights, algo = "sha256")
    if (file.exists(weights_path)) {
        old <- readRDS(weights_path)
        if (!identical(old$signature, signature) || !identical(old$weights_sha256, weight_hash) || !identical(old$weights, 
            weights)) 
            stop("Bootstrap weight checkpoint mismatch.")
    }
    else v7_save(list(signature = signature, weights_sha256 = weight_hash, weights = weights, seed = cfg$seed, 
        RNG_kind = RNGkind()), weights_path)
    all_t <- array(NA_real_, c(15L, length(prepared$probability), cfg$draws))
    timings <- numeric()
    min_se <- Inf
    starts <- seq.int(1L, cfg$draws, by = cfg$batch_size)
    for (s in starts) {
        e <- min(cfg$draws, s + cfg$batch_size - 1L)
        path <- file.path(directory, sprintf("draws_%05d_%05d.rds", s, e))
        if (file.exists(path)) {
            batch <- readRDS(path)
            if (!identical(batch$signature, signature) || !identical(batch$weights_sha256, weight_hash) || 
                !identical(batch$range, c(s, e))) 
                stop("Bootstrap batch checkpoint mismatch.")
        }
        else {
            cat("V7", cfg$vi, "wild bootstrap-t draws", s, "-", e, "/", cfg$draws, "\n")
            flush.console()
            elapsed <- system.time(result <- v7_wild_batch(records, point, weights[, s:e, drop = FALSE]))
            batch <- list(signature = signature, weights_sha256 = weight_hash, range = c(s, e), result = result, 
                elapsed = elapsed)
            v7_save(batch, path)
        }
        all_t[, , s:e] <- batch$result$t_abs
        timings <- c(timings, batch$elapsed[["elapsed"]])
        min_se <- min(min_se, batch$result$minimum_studentizing_se)
    }
    joint <- apply(all_t, 3, max)
    curve <- apply(all_t, c(1, 3), max)
    critical_joint <- as.numeric(stats::quantile(joint, 0.95, type = 8, names = FALSE))
    critical_pointwise <- apply(all_t, c(1, 2), stats::quantile, probs = 0.95, type = 8, names = FALSE)
    critical_curve <- apply(curve, 1, stats::quantile, probs = 0.95, type = 8, names = FALSE)
    if (critical_joint + 1e-12 < max(critical_pointwise)) 
        stop("Joint critical smaller than a pointwise critical.")
    sorted <- sort(joint)
    ranks <- pmax(1L, pmin(cfg$draws, c(stats::qbinom(0.025, cfg$draws, 0.95), stats::qbinom(0.975, cfg$draws, 
        0.95) + 1L)))
    half <- floor(cfg$draws/2)
    list(critical_joint = critical_joint, critical_pointwise = critical_pointwise, critical_curve = critical_curve, 
        joint_max = joint, per_curve_max = curve, draws = cfg$draws, weights_sha256 = weight_hash, joint_critical_MC95_interval = sorted[ranks], 
        MC_interval_order_ranks = ranks, joint_critical_half_samples = c(stats::quantile(joint[seq_len(half)], 
            0.95, names = FALSE, type = 8), stats::quantile(joint[(half + 1L):cfg$draws], 0.95, names = FALSE, 
            type = 8)), expected_upper_tail_draws = cfg$draws * 0.05, batch_elapsed_seconds = timings, 
        minimum_studentizing_se = min_se, scope = "Joint family: fifteen TM-1 contrasts across this VI's observed quantile grid; ordinary pointwise95% also exported. Fixed working weights/smoothers; CR2 studentizer recomputed per draw.")
}

v7_tables <- function(point, boot, cfg) {
    curves <- do.call(rbind, lapply(seq_along(point$targets), function(r) {
        nm <- point$targets[r]
        estimate <- point$beta_smooth[nm, ]
        se <- point$se[r, ]
        data.frame(VI = cfg$vi, Coefficient = nm, Reference = "5016_TM1", Probability = point$probability, 
            Estimate = as.numeric(estimate), RawEstimate = as.numeric(point$beta_raw[nm, ]), CR2_SE = se, 
            LowerPointwise95 = estimate - boot$critical_pointwise[r, ] * se, UpperPointwise95 = estimate + 
                boot$critical_pointwise[r, ] * se, LowerCurve95 = estimate - boot$critical_curve[r] * 
                se, UpperCurve95 = estimate + boot$critical_curve[r] * se, LowerFamily15_95 = estimate - 
                boot$critical_joint * se, UpperFamily15_95 = estimate + boot$critical_joint * se, PointwiseCritical = boot$critical_pointwise[r, 
                ], Family15Critical = boot$critical_joint, Status = if (cfg$mode == "smoke") 
                "SMOKE_ONLY_NOT_MANUSCRIPT_RESULTS"
            else "APPROXIMATE_CR2_PLANT_WILD_BOOTSTRAP_T_FIXED_WEIGHTS_AND_SMOOTHERS")
    }))
    means <- do.call(rbind, lapply(c("5016", as.character(5001:5015)), function(g) {
        estimate <- point$beta_smooth["(Intercept)", ]
        if (g != "5016") 
            estimate <- estimate + point$beta_smooth[paste0("Pedigree", g), ]
        data.frame(VI = cfg$vi, Genotype = g, Probability = point$probability, AdjustedMean = as.numeric(estimate), 
            Status = "POINT_ESTIMATES_ONLY_EQUAL_DATE_EQUAL_BLOCK_REFERENCE_NO_MEAN_CONFIDENCE_BANDS")
    }))
    rownames(curves) <- rownames(means) <- NULL
    list(curves = curves, means = means)
}

v7_validate <- function(out, expected_vi = NULL, source_path = NULL) {
    manifest_path <- file.path(out, "manifest.rds")
    completion_path <- file.path(out, "completion_status.json")
    if (!file.exists(manifest_path)) 
        stop("ABSENT_OR_INVALID: no V7 manifest.")
    manifest <- readRDS(manifest_path)
    cfg <- manifest$provenance$configuration
    if (!is.null(expected_vi) && !identical(cfg$vi, expected_vi)) 
        stop("VI mismatch.")
    if (!is.null(source_path)) 
        cfg$script_path <- normalizePath(source_path, winslash = "/", mustWork = TRUE)
    cfg <- v7_setup(cfg)
    if (!identical(v7_signature(v7_provenance(cfg)), manifest$signature)) 
        stop("INVALID: source/input/map/package/config identity changed.")
    if (!file.exists(completion_path)) 
        stop("INCOMPLETE: no committed completion record.")
    done <- jsonlite::read_json(completion_path, simplifyVector = TRUE)
    if (!identical(done$status, "COMPLETE") || !identical(done$signature, manifest$signature)) 
        stop("Invalid completion record.")
    if (cfg$mode != "full" || cfg$draws < 1999L) 
        stop("NOT_PRODUCTION: smoke/inspection is not a full inference run.")
    for (nm in names(done$output_sha256)) {
        p <- file.path(out, nm)
        if (!file.exists(p) || !identical(tolower(v7_hash(p)), tolower(done$output_sha256[[nm]]))) 
            stop("Output hash mismatch: ", nm)
    }
    model <- readRDS(file.path(out, "functional_model_V7.rds"))
    tab <- as.data.frame(data.table::fread(file.path(out, "genotype_contrasts_V7.csv")))
    if (length(model$point$probability) != 100L || nrow(tab) != 1500L || max(abs(as.numeric(model$point$probability) - 
        seq(0.01, 1, 0.01))) > 1e-10 || !setequal(unique(tab$Coefficient), paste0("Pedigree", 5001:5015)) || 
        anyDuplicated(tab[, c("Coefficient", "Probability")]) || any(!is.finite(as.matrix(tab[, vapply(tab, 
        is.numeric, logical(1)), drop = FALSE])))) 
        stop("Incomplete/invalid family grid.")
    if (any(tab$LowerFamily15_95 > tab$LowerPointwise95 + 1e-10) || any(tab$UpperFamily15_95 < tab$UpperPointwise95 - 
        1e-10)) 
        stop("Band nesting check failed.")
    if (!identical(as.integer(model$input_design$full_counts), c(47807L, 479L, 101L, 16L, 10L))) 
        stop("Cohort counts changed.")
    point_files <- file.path(out, "point_checkpoints", sprintf("quantile_%03d.rds", 1:100))
    if (!all(file.exists(point_files))) 
        stop("Missing point checkpoints.")
    for (j in 1:100) {
        z <- readRDS(point_files[j])
        if (!isTRUE(z$ok) || !identical(z$signature, manifest$signature) || abs(z$probability - j/100) > 
            1e-10 || !identical(z$components$clusters, model$input_design$plant_ids)) 
            stop("Invalid point checkpoint.")
    }
    starts <- seq.int(1L, cfg$draws, by = cfg$batch_size)
    for (s in starts) {
        e <- min(cfg$draws, s + cfg$batch_size - 1L)
        p <- file.path(out, "wild_checkpoints", sprintf("draws_%05d_%05d.rds", s, e))
        if (!file.exists(p)) 
            stop("Missing wild bootstrap batch.")
        z <- readRDS(p)
        if (!identical(z$signature, manifest$signature) || !identical(z$weights_sha256, model$bootstrap$weights_sha256) || 
            !identical(z$range, c(s, e)) || any(!is.finite(z$result$t_abs))) 
            stop("Invalid wild bootstrap batch.")
    }
    list(status = "COMPLETE_VALIDATED", vi = cfg$vi, quantiles = 100L, contrasts = 15L, draws = cfg$draws, 
        signature = manifest$signature)
}

v7_run <- function(cfg) {
    cfg <- v7_setup(cfg)
    if (cfg$mode == "validate") 
        return(v7_validate(cfg$out, cfg$vi, cfg$script_path))
    if (cfg$mode == "inspect") {
        prepared <- v7_prepare(cfg)
        print(prepared$summary[c("full_counts", "observations", "dates", "fixed_rank")])
        return(invisible(prepared$summary))
    }
    if (dir.exists(cfg$out) && !isTRUE(cfg$resume)) 
        stop("Output exists; choose new output or explicit identical resume.")
    if (!dir.exists(cfg$out)) 
        dir.create(cfg$out, recursive = TRUE)
    lock <- file.path(cfg$out, ".active_owner")
    if (!dir.create(lock, showWarnings = FALSE)) 
        stop("Owner marker exists; verify no active writer before manual recovery.")
    owner <- file.path(lock, "owner.rds")
    saveRDS(list(pid = Sys.getpid(), started = Sys.time(), script = cfg$script_path), owner)
    on.exit({
        if (identical(normalizePath(dirname(lock), winslash = "/"), normalizePath(cfg$out, winslash = "/")) && 
            identical(list.files(lock, all.files = TRUE, no.. = TRUE), "owner.rds")) {
            file.remove(owner)
            unlink(lock, recursive = TRUE)
        }
    }, add = TRUE)
    provenance <- v7_provenance(cfg)
    signature <- v7_signature(provenance)
    path <- file.path(cfg$out, "manifest.rds")
    if (file.exists(path)) {
        if (!identical(readRDS(path)$signature, signature)) 
            stop("Resume identity mismatch.")
    }
    else {
        if (isTRUE(cfg$resume) || length(setdiff(list.files(cfg$out, all.files = TRUE, no.. = TRUE), 
            ".active_owner"))) 
            stop("Refusing to adopt unmanifested output.")
        v7_save(list(signature = signature, provenance = provenance, session = utils::sessionInfo()), 
            path)
    }
    prepared <- v7_prepare(cfg)
    design_path <- file.path(cfg$out, "input_design.rds")
    if (!file.exists(design_path)) 
        v7_save(prepared$summary, design_path)
    if (!identical(readRDS(design_path), prepared$summary)) 
        stop("Input/design signature changed.")
    print(prepared$summary[c("full_counts", "observations", "dates", "fixed_rank")])
    if (cfg$mode == "inspect") 
        return(invisible(prepared$summary))
    final <- file.path(cfg$out, "functional_model_V7.rds")
    if (file.exists(file.path(cfg$out, "completion_status.json"))) {
        if (cfg$mode == "full") 
            return(v7_validate(cfg$out, cfg$vi, cfg$script_path))
        return(invisible(readRDS(final)))
    }
    records <- v7_point(prepared, cfg, signature)
    point <- v7_smooth(records, prepared$probability, cfg)
    point_path <- file.path(cfg$out, "point_estimates_V7.rds")
    if (!file.exists(point_path)) 
        v7_save(list(signature = signature, point = point), point_path)
    boot <- v7_bootstrap(records, point, prepared, cfg, signature)
    tables <- v7_tables(point, boot, cfg)
    v7_csv(tables$curves, file.path(cfg$out, "genotype_contrasts_V7.csv"))
    v7_csv(tables$means, file.path(cfg$out, "genotype_adjusted_means_point_only_V7.csv"))
    diagnostics <- lapply(records, function(z) list(probability = z$probability, singular = z$singular, 
        warnings = z$warnings, convergence = z$convergence, fit_timing = z$fit_timing, covariance_timing = z$covariance_timing, 
        raw_contrast_df = z$components$raw_contrast_df, Gram_relative_error = z$components$Gram_relative_error, 
        contribution_effective_count = z$components$contribution_effective_count, maximum_absolute_adjustment_element = z$components$maximum_absolute_adjustment_element))
    result <- list(signature = signature, formula = v7_formula(), point = point, bootstrap = boot, point_diagnostics = diagnostics, 
        input_design = prepared$summary, provenance = provenance, scope = "Fifteen TM-1 contrasts jointly over the observed quantile grid, separately for each VI; fixed observed block/date adjustment.", 
        limitations = "Approximate fixed-weight/fixed-smoother CR2 Plant wild bootstrap-t. Assumes independent plants after adjustment; permits arbitrary within-plant date/quantile covariance. Does not propagate variance-component or smoothing-parameter estimation, remove smoothing bias, validate spatial independence, or establish genotype-by-date interactions.")
    if (!file.exists(final)) 
        v7_save(result, final)
    outputs <- c("functional_model_V7.rds", "genotype_contrasts_V7.csv", "genotype_adjusted_means_point_only_V7.csv", 
        "point_estimates_V7.rds", "input_design.rds", "wild_weights_V7.rds")
    done <- list(status = if (cfg$mode == "full") "COMPLETE" else "SMOKE_COMPLETE", signature = signature, 
        VI = cfg$vi, quantiles = length(prepared$probability), contrasts = 15L, draws = cfg$draws, completed_UTC = format(Sys.time(), 
            tz = "UTC", usetz = TRUE), output_sha256 = as.list(setNames(vapply(file.path(cfg$out, outputs), 
            v7_hash, character(1)), outputs)))
    done_path <- file.path(cfg$out, "completion_status.json")
    tmp <- paste0(done_path, ".partial_", Sys.getpid())
    if (file.exists(done_path) || file.exists(tmp)) 
        stop("Completion record already exists.")
    jsonlite::write_json(done, tmp, auto_unbox = TRUE, pretty = TRUE)
    if (!file.rename(tmp, done_path)) 
        stop("Completion commit failed.")
    cat("V7", cfg$vi, done$status, "joint family critical", boot$critical_joint, "output", cfg$out, "\n")
    invisible(result)
}


v7_setup <- function (cfg) 
{
    if (!cfg$vi %in% c("ExG2", "PSRI", "CVI")) 
        stop("Only the three unique VIs ExG2/PSRI/CVI are supported.")
    if (!cfg$mode %in% c("inspect", "smoke", "full", "validate") || cfg$cores != 1L || !identical(as.numeric(cfg$alpha), 
        0.05) || !is.finite(cfg$draws) || cfg$draws%%1 != 0 || cfg$draws < 39L || !is.finite(cfg$batch_size) || 
        cfg$batch_size%%1 != 0 || cfg$batch_size < 1L || !is.finite(cfg$seed) || cfg$seed < 1L || cfg$seed > 
        .Machine$integer.max - 100L) 
        stop("Invalid V7 configuration.")
    if (cfg$mode == "full" && cfg$draws < 1999L) 
        stop("Full inference requires at least 1999 draws; default4999.")
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1", 
        NUMEXPR_NUM_THREADS = "1")
    for (p in c("lme4", "mgcv", "Matrix", "clubSandwich", "data.table", "digest", "jsonlite")) if (!requireNamespace(p, 
        quietly = TRUE)) 
        stop("Missing package: ", p, "; install the version listed in environment/R_versions.csv.")
    if (as.character(utils::packageVersion("clubSandwich")) != "0.7.0") 
        stop("Revalidate a new script version before changing clubSandwich0.7.0.")
    data.table::setDTthreads(1L)
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) 
        RhpcBLASctl::blas_set_num_threads(1L)
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    cfg$input <- normalizePath(cfg$input, winslash = "/", mustWork = TRUE)
    cfg$fieldmap <- normalizePath(cfg$fieldmap, winslash = "/", mustWork = TRUE)
    cfg$out <- normalizePath(cfg$out, winslash = "/", mustWork = FALSE)
    cfg
}
v7_provenance <- function (cfg) 
{
    cfg$resume <- NULL
    packages <- c("lme4", "mgcv", "Matrix", "clubSandwich", "data.table", "digest", "jsonlite")
    internals <- c("vcov_CR", "CR2", "vcovCR.lmerMod", "targetVariance.lmerMod", "weightMatrix.lmerMod")
    list(configuration = cfg, script_sha256 = v7_hash(cfg$script_path), input_sha256 = v7_hash(cfg$input), 
        fieldmap_sha256 = v7_hash(cfg$fieldmap), R_version = R.version.string, RNG_kind = RNGkind(), 
        packages = setNames(vapply(packages, function(p) as.character(utils::packageVersion(p)), character(1)), 
            packages), CR2_internals_sha256 = setNames(vapply(internals, function(nm) digest::digest(deparse(get(nm, 
            asNamespace("clubSandwich"))), algo = "sha256"), character(1)), internals), engine_sha256 = v7_hash(cfg$engine_path), io_sha256 = v7_hash(cfg$io_path))
}
