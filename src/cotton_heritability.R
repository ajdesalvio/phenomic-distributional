block_only_formula_v5 <- function() "response ~ 1 + (1 | Pedigree) + (1 | Block)"

quantile_id_v5 <- function(p) {
    if (!is.numeric(p) || any(!is.finite(p))) 
        stop("Probability must be finite numeric")
    id <- as.integer(round(100 * p))
    if (any(id < 1L | id > 100L) || any(abs(p - id/100) > 1e-10)) 
        stop("Probability must be on the hundredth grid")
    id
}

cohort_hash_v5 <- function(x, response) {
    o <- order(as.character(x$Plant))
    fields <- lapply(x[, c("Plant", "Pedigree", "Block"), drop = FALSE], as.character)
    if (any(vapply(fields, function(z) any(grepl("[|\n\r]", z)), logical(1)))) 
        stop("Invalid metadata delimiter")
    rows <- do.call(paste, c(lapply(fields, function(z) z[o]), list(sprintf("%a", response[o] + 0), sep = "|")))
    digest::digest(paste(rows, collapse = "\n"), algo = "sha256", serialize = FALSE)
}

v6_hash <- function(p) digest::digest(file = p, algo = "sha256", serialize = FALSE)

v6_save <- function(x, p) {
    if (file.exists(p)) 
        stop("Refuse overwrite: ", p)
    t <- paste0(p, ".partial_", Sys.getpid())
    if (file.exists(t)) 
        stop("Partial exists: ", t)
    saveRDS(x, t)
    if (!file.rename(t, p)) 
        stop("Could not finalize ", p)
}

v6_input <- function(path, vi, field) {
    e <- new.env(parent = emptyenv())
    load(path, envir = e)
    stopifnot(all(c("data", "data2") %in% ls(e)))
    d <- data.table::as.data.table(e$data2)
    raw <- data.table::as.data.table(e$data)
    meta <- c("Test", "Tier", "Row", "Range", "Pedigree", "Pltg_ID", "DAT", "Replicate")
    stopifnot(all(meta %in% names(d)), all(c(meta, "vi_name", "probability", "quantile") %in% names(raw)))
    cols <- grep(paste0("^", vi, "_Quant_"), names(d), value = TRUE)
    ids <- quantile_id_v5(as.numeric(sub(paste0("^", vi, "_Quant_"), "", cols)))
    stopifnot(length(ids) == 100L, !anyDuplicated(ids), setequal(ids, 1:100))
    cols <- cols[match(1:100, ids)]
    dates <- sort(unique(d$DAT))
    stopifnot(length(dates) == 101L, !any(dates %in% c(46, 71, 72, 99, 100)), all(is.finite(dates)))
    take <- d$Pltg_ID != "5011_018"
    x <- as.data.frame(d[take, ..meta])
    y <- as.matrix(d[take, ..cols])
    stopifnot(nrow(x) == 47807L, !anyNA(x), all(is.finite(y)), !anyDuplicated(x[, c("DAT", "Pltg_ID")]))
    fi <- match(x$Pltg_ID, field$Pltg_ID)
    for (n in c("Test", "Tier", "Row", "Range", "Pedigree", "Pltg_ID")) if (anyNA(fi) || !identical(as.character(x[[n]]), 
        as.character(field[[n]][fi]))) 
        stop("Field identity failed: ", n)
    stopifnot(all(as.integer(as.character(x$Replicate)) == as.integer(sub(".*_", "", x$Pltg_ID))))
    r <- raw[DAT %in% dates & Pltg_ID != "5011_018", c(meta, "vi_name", "probability", "quantile"), with = FALSE]
    r[, `:=`(QuantileID, quantile_id_v5(probability))]
    stopifnot(!anyNA(r), all(is.finite(r$quantile)), all(r$vi_name == vi), nrow(r) == nrow(x) * 100L, 
        !anyDuplicated(r[, .(DAT, Pltg_ID, QuantileID)]))
    r[, `:=`(WideRow, match(paste(DAT, Pltg_ID, sep = "|"), paste(x$DAT, x$Pltg_ID, sep = "|")))]
    stopifnot(!anyNA(r$WideRow))
    data.table::setorder(r, QuantileID, WideRow)
    stopifnot(identical(as.numeric(r$quantile), as.numeric(y)))
    for (n in meta) stopifnot(identical(as.character(r[[n]]), rep(as.character(x[[n]]), 100L)))
    x$Plant <- factor(x$Pltg_ID)
    x$Pedigree <- factor(x$Pedigree)
    x$Block <- interaction(x$Test, x$Tier, drop = TRUE, sep = "_")
    stopifnot(nlevels(x$Plant) == 479L, nlevels(x$Pedigree) == 16L, nlevels(x$Block) == 10L)
    counts <- data.table::as.data.table(x)[, .(N = .N, Genotypes = data.table::uniqueN(Pedigree), Blocks = data.table::uniqueN(Block)), 
        by = DAT]
    stopifnot(all(counts$Genotypes == 16L), all(counts$Blocks == 10L))
    oldcounts <- c(`57` = 476L, `86` = 478L, `96` = 479L, `113` = 479L, `130` = 477L)
    stopifnot(identical(as.integer(table(x$DAT)[names(oldcounts)]), unname(oldcounts)))
    list(data = x, Y = y, dates = dates, counts = counts, identity = list(input = path, VI = vi, rows = nrow(x), 
        plants = 479L, dates = dates, excluded_dates = c(46, 71, 72, 99, 100), excluded_plant = "5011_018", 
        no_imputation = TRUE, all_long_wide_responses_and_metadata_identical = TRUE, field_identity = TRUE, 
        source_quantile_rows = nrow(r)))
}

v6_fit_date <- function(task) {
    data.table::setDTthreads(1L)
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) 
        RhpcBLASctl::blas_set_num_threads(1L)
    started <- Sys.time()
    p <- readRDS(task$cache)
    take <- p$data$DAT == task$DAT
    x <- droplevels(p$data[take, , drop = FALSE])
    stopifnot(!anyDuplicated(x$Plant), nlevels(x$Pedigree) == 16L, nlevels(x$Block) == 10L)
    records <- metrics <- components <- vector("list", length(task$qids))
    for (k in seq_along(task$qids)) {
        q <- task$qids[k]
        x$response <- p$Y[take, q]
        warnings <- messages <- character()
        err <- NULL
        h <- cohort_hash_v5(x, x$response)
        fit <- tryCatch(withCallingHandlers(lme4::lmer(stats::as.formula(block_only_formula_v5(), env = baseenv()), 
            data = x, REML = TRUE, control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+05))), 
            warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            }, message = function(m) {
                messages <<- c(messages, conditionMessage(m))
                invokeRestart("muffleMessage")
            }), error = function(e) {
            err <<- conditionMessage(e)
            NULL
        })
        if (is.null(fit)) 
            stop("Fit failed ", task$VI, " DAT", task$DAT, " Q", q, ": ", err)
        vc <- as.data.frame(lme4::VarCorr(fit))
        stopifnot(nrow(vc) == 3L, setequal(vc$grp, c("Pedigree", "Block", "Residual")))
        vg <- vc$vcov[vc$grp == "Pedigree"]
        ve <- vc$vcov[vc$grp == "Residual"]
        vb <- vc$vcov[vc$grp == "Block"]
        stopifnot(all(is.finite(vc$vcov)), all(vc$vcov >= 0), vg + ve > 0, stats::nobs(fit) == nrow(x), 
            identical(names(lme4::fixef(fit)), "(Intercept)"))
        conv <- fit@optinfo$conv
        bad <- c(warnings, messages, unlist(conv$lme4$messages))
        bad <- bad[!grepl("singular", bad, ignore.case = TRUE)]
        ok <- !length(bad) && (is.null(conv$opt) || all(conv$opt == 0))
        singular <- lme4::isSingular(fit)
        records[[k]] <- list(VI = task$VI, DAT = task$DAT, QuantileID = q, Probability = q/100, formula = block_only_formula_v5(), 
            cohort_sha256 = h, nobs = nrow(x), variance = vc, fixed_intercept = unname(lme4::fixef(fit)), 
            converged = ok, singular = singular, convergence = conv, warnings = warnings, messages = messages, 
            optimizer = "bobyqa", maxfun = 200000L, REML = TRUE)
        metrics[[k]] <- data.frame(Vegetation.Index = task$VI, DAT = task$DAT, QuantileID = q, Probability = q/100, 
            Vg = vg, Vb = vb, Ve = ve, H2Individual = vg/(vg + ve), Model = "Actual_Block_Only_Focal_FullSeason_V6", 
            Formula = block_only_formula_v5(), N_Observations = nrow(x), N_Genotypes = 16L, N_Blocks = 10L, 
            Converged = ok, Singular = singular, ExactZeroComponents = paste(vc$grp[vc$vcov == 0], collapse = ";"), 
            Warnings = paste(warnings, collapse = " | "), ConvergenceMessages = paste(unlist(conv$lme4$messages), 
                collapse = " | "), CohortSHA256 = h)
        vc$Percent <- 100 * vc$vcov/sum(vc$vcov)
        vc$Vegetation.Index <- task$VI
        vc$DAT <- task$DAT
        vc$QuantileID <- q
        vc$Probability <- q/100
        vc$H2Individual <- vg/(vg + ve)
        components[[k]] <- vc
    }
    result <- list(task = task, records = records, metrics = data.table::rbindlist(metrics), components = data.table::rbindlist(components), 
        elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
    v6_save(result, task$output)
    if (!all(result$metrics$Converged)) 
        stop("Unresolved convergence warning in ", task$output)
    list(path = task$output, fits = nrow(result$metrics), seconds = result$elapsed_seconds)
}

