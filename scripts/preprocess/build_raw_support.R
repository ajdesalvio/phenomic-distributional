# Unweighted finite image-quantile means and counts; no fitting or yield input.
# Rscript build_raw_support.R --raw CSV --predictors CSV --output NEW_CSV
support_envs <- function() c("IAH4.2021", "MIH1.2020", "MNH1.2020", "MNH1.2021", "TXH1.2020", "TXH2.2021", 
    "WIH1.2021", "WIH2.2020", "WIH2.2021", "WIH3.2020", "WIH3.2021")

support_build <- function(raw, blue, envs = support_envs()) {
    requireNamespace("data.table", quietly = TRUE)
    raw <- data.table::copy(data.table::as.data.table(raw))
    blue <- data.table::copy(data.table::as.data.table(blue))
    rc <- c("file_name", "vi_name", "probability", "quantile", "Env", "Pedigree", "DAP")
    keys <- c("Env", "Pedigree", "DAP", "Probability", "Vegetation.Index")
    if (!all(rc %in% names(raw)) || !all(keys %in% names(blue))) 
        stop("Required raw or BLUE columns missing")
    raw <- raw[Env %in% envs & vi_name == "ExGR"]
    blue <- blue[Env %in% envs & Vegetation.Index == "ExGR", ..keys]
    if (!nrow(raw) || !nrow(blue)) 
        stop("Empty requested input")
    if (anyNA(raw[, ..rc][, c("file_name", "Env", "Pedigree", "DAP", "probability"), with = FALSE]) || 
        any(!nzchar(raw$file_name) | !nzchar(raw$Env) | !nzchar(raw$Pedigree))) 
        stop("Missing raw identity/feature key")
    if (anyNA(blue) || any(!nzchar(blue$Env) | !nzchar(blue$Pedigree))) 
        stop("Missing BLUE identity/feature key")
    if (!is.numeric(raw$quantile) || !is.numeric(raw$probability) || !is.numeric(blue$Probability) || 
        !is.numeric(raw$DAP) || !is.numeric(blue$DAP)) 
        stop("Nonnumeric response/quantile/date")
    valid_prob <- function(z) all(is.finite(z) & z >= 0.01 & z <= 1 & abs(100 * z - round(100 * z)) < 
        1e-07)
    if (!valid_prob(raw$probability) || !valid_prob(blue$Probability) || any(!is.finite(raw$DAP)) || 
        any(!is.finite(blue$DAP))) 
        stop("Invalid quantile/date key")
    raw[, `:=`(probability, round(probability, 2))]
    blue[, `:=`(Probability, round(Probability, 2))]
    if (anyDuplicated(raw[, .(file_name, vi_name, probability)])) 
        stop("Duplicate image/VI/quantile rows")
    image_map <- unique(raw[, .(file_name, Env, Pedigree, DAP)])
    if (anyDuplicated(image_map$file_name)) 
        stop("One image maps to multiple genotype/environment/date keys")
    if (anyDuplicated(blue)) 
        stop("Duplicate BLUE keys")
    support <- raw[, .(Quantile.Raw.Mean = if (any(is.finite(quantile))) 
        mean(quantile[is.finite(quantile)])
    else NA_real_, N.Raw.Finite = sum(is.finite(quantile))), by = .(Env, Pedigree, DAP, Probability = probability, 
        Vegetation.Index = vi_name)]
    missing <- blue[!support, on = keys]
    if (nrow(missing)) 
        stop("Current raw input lacks ", nrow(missing), " BLUE cells; no support file written. First: ", 
            paste(unlist(missing[1]), collapse = " | "))
    result <- support[blue, on = keys]
    data.table::setorderv(result, keys)
    stopifnot(nrow(result) == nrow(blue), !anyDuplicated(result[, ..keys]), all((result$N.Raw.Finite > 
        0) == is.finite(result$Quantile.Raw.Mean)))
    result
}

support_main <- function(args = commandArgs(trailingOnly = TRUE)) {
    suppressPackageStartupMessages(library(data.table))
    setDTthreads(1)
    own <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
    cfg <- list(raw = NULL, predictors = NULL, output = NULL)
    if (!length(args) || length(args)%%2) 
        stop("Arguments require --name value pairs, including --raw <current raw CSV>")
    seen <- character()
    for (i in seq.int(1, length(args), by = 2)) {
        k <- sub("^--", "", args[i])
        if (!k %in% names(cfg) || k %in% seen) 
            stop("Unknown option: ", args[i])
        cfg[[k]] <- args[i + 1]
        seen <- c(seen, k)
    }
    if (any(vapply(cfg, is.null, logical(1))) || !file.exists(cfg$raw) || !file.exists(cfg$predictors)) 
        stop("Supply existing current raw CSV with --raw; predictor CSV must exist")
    cfg <- lapply(cfg, function(p) normalizePath(p, winslash = "/", mustWork = FALSE))
    manifest <- paste0(cfg$output, ".provenance.rds")
    summary <- paste0(cfg$output, ".provenance.txt")
    paths <- c(cfg$output, manifest, summary)
    if (any(file.exists(c(paths, paste0(paths, ".partial"))))) 
        stop("Refusing to overwrite output/provenance/partial")
    hashes <- tools::md5sum(c(cfg$raw, cfg$predictors, own))
    raw <- fread(cfg$raw, select = c("file_name", "vi_name", "probability", "quantile", "Env", "Pedigree", 
        "DAP"), showProgress = FALSE)
    blue <- fread(cfg$predictors, select = c("Env", "Pedigree", "DAP", "Probability", "Vegetation.Index"), 
        showProgress = FALSE)
    result <- support_build(raw, blue)
    hashes_after <- tools::md5sum(c(cfg$raw, cfg$predictors, own))
    if (!identical(hashes, hashes_after)) 
        stop("Source input/script changed while building; no output written")
    dir.create(dirname(cfg$output), recursive = TRUE, showWarnings = FALSE)
    fwrite(result, paste0(cfg$output, ".partial"))
    prov <- list(version = "V1", created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), configuration = cfg, 
        source_md5 = hashes, source_raw_rows = nrow(raw), source_BLUE_rows = nrow(blue), exported_rows = nrow(result), 
        environments = support_envs(), session = capture.output(sessionInfo()), method = "Unweighted mean of finite image-level quantiles per Env/Pedigree/DAP/Probability/VI; no BLUEs, yields, folds or outcomes used in support values.", 
        raw_provenance_qualification = "Caller must provide the current, provenance-confirmed plot quantile export; key coverage alone does not prove image processing provenance.")
    saveRDS(prov, paste0(manifest, ".partial"))
    writeLines(c(prov$method, prov$raw_provenance_qualification, paste("Exported cells:", nrow(result)), 
        paste(names(hashes), unname(hashes), sep = " MD5=")), paste0(summary, ".partial"))
    for (p in paths) if (!file.rename(paste0(p, ".partial"), p)) 
        stop("Could not finalize ", p)
    cat("Built", nrow(result), "support cells:", cfg$output, "\n")
}

if (sys.nframe() == 0L) support_main()

