suppressPackageStartupMessages(library(data.table))

setDTthreads(2)

model_keys <- c("Env", "DAP", "Probability", "Vegetation.Index")

blue_columns <- c("Env", "Pedigree", "Probability", "Quantile.BLUE", "DAP", "Vegetation.Index")

blue_keys <- setdiff(blue_columns, "Quantile.BLUE")

options_cli <- function(required, optional = character()) {
    a <- commandArgs(TRUE)
    if (any(!grepl("^--[^=]+=.+$", a))) 
        stop("Use --name=value arguments")
    k <- sub("^--([^=]+)=.*$", "\\1", a)
    if (anyDuplicated(k) || any(!k %in% c(required, optional))) 
        stop("Unknown/duplicate argument")
    x <- setNames(as.list(sub("^--[^=]+=", "", a)), k)
    if (any(!required %in% k)) 
        stop("Required: ", paste(required, collapse = ", "))
    x
}

new_directory <- function(p) {
    if (file.exists(p)) 
        stop("Output exists: ", p)
    if (!dir.create(p, recursive = TRUE)) 
        stop("Cannot create: ", p)
}

hash <- function(p) unname(tools::md5sum(p))

valid_keys <- function(x, keys) {
    if (!all(keys %in% names(x)) || anyNA(x[, ..keys]) || anyDuplicated(x[, ..keys])) 
        stop("Missing/duplicate keys")
}

load_plan <- function(path) {
    p <- readRDS(path)
    valid_keys(p$keys, model_keys)
    if (any(!p$keys$Vegetation.Index %in% "ExGR")) 
        stop("Non-ExGR plan")
    p
}

compare_blue <- function(x, reference, output, tolerance = 1e-07) {
    y <- fread(reference, select = blue_columns)
    valid_keys(x, blue_keys)
    valid_keys(y, blue_keys)
    setorderv(x, blue_keys)
    setorderv(y, blue_keys)
    same_keys <- identical(x[, ..blue_keys], y[, ..blue_keys])
    if (!same_keys) 
        stop("Reference key mismatch")
    if (any(!is.finite(x$Quantile.BLUE)) || any(!is.finite(y$Quantile.BLUE))) 
        stop("Nonfinite BLUE")
    delta <- abs(x$Quantile.BLUE - y$Quantile.BLUE)
    qa <- data.table(Rows = nrow(x), MaxAbsoluteDifference = max(delta), Tolerance = tolerance, Pass = all(delta <= 
        tolerance))
    fwrite(qa, file.path(output, "Reference_Comparison.csv"))
    if (!qa$Pass) 
        stop("Reference values differ beyond tolerance")
    qa
}

