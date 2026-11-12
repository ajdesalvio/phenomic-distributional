# Recover the archived BLUE table without fitting. All paths are explicit.
# Rscript this.R --archive FILE --reference CSV --output NEW_DIR
recover_blue_reference <- function(archive, reference, output, script) {
    suppressPackageStartupMessages(library(data.table))
    setDTthreads(2)
    if (file.exists(output) || dir.exists(output)) 
        stop("Output exists; choose a new numbered directory")
    sources <- normalizePath(c(archive, reference, script), winslash = "/", mustWork = TRUE)
    before <- tools::md5sum(sources)
    e <- new.env(parent = baseenv())
    load(sources[1], envir = e)
    if (!exists("blue.filt", envir = e, inherits = FALSE)) 
        stop("Archive lacks blue.filt")
    cols <- c("Env", "Pedigree", "Probability", "Quantile.BLUE", "DAP", "Vegetation.Index")
    x <- as.data.table(e$blue.filt)
    if (!all(cols %in% names(x))) 
        stop("Archive BLUE schema differs")
    x <- copy(x[, ..cols])
    rm(e)
    gc()
    keys <- setdiff(cols, "Quantile.BLUE")
    if (anyNA(x[, ..keys]) || anyDuplicated(x[, ..keys]) || any(x$Vegetation.Index != "ExGR")) 
        stop("Invalid BLUE identity")
    dir.create(output, recursive = TRUE)
    candidate <- file.path(output, "Quantile_BLUEs_G2F_Distributional_ExGR_Reference.csv.partial")
    fwrite(x, candidate)
    candidate_md5 <- unname(tools::md5sum(candidate))
    equal <- identical(candidate_md5, unname(before[2]))
    stable <- identical(before, tools::md5sum(sources))
    if (!equal || !stable) 
        stop("Recovery identity failed; only a .partial candidate remains")
    final <- sub("[.]partial$", "", candidate)
    if (!file.rename(candidate, final)) 
        stop("Cannot finalize recovered CSV")
    jsonlite::write_json(list(status = "PASS", object = "blue.filt", rows = nrow(x), columns = cols, 
        input_md5 = as.list(before), output = list(path = normalizePath(final, winslash = "/"), md5 = candidate_md5), 
        reference_byte_identical = equal, inputs_unchanged = stable, session = capture.output(sessionInfo()), 
        scope = "Recovery of a saved BLUE table; no first-stage model refit or claim of original fit-batch provenance."), 
        file.path(output, "BLUE_Recovery_Manifest.json"), pretty = TRUE, auto_unbox = TRUE)
    message("PASS: byte-identical recovery of ", nrow(x), " BLUE rows")
}

if (sys.nframe() == 0L) {
    args <- commandArgs(TRUE)
    cfg <- list(archive = NULL, reference = NULL, output = NULL)
    seen <- character()
    if (length(args) != 6L) 
        stop("Supply --archive FILE --reference CSV --output NEW_DIR")
    for (i in seq(1, 6, 2)) {
        n <- sub("^--", "", args[i])
        if (!n %in% names(cfg) || n %in% seen) 
            stop("Unknown/repeated argument")
        cfg[[n]] <- args[i + 1]
        seen <- c(seen, n)
    }
    own <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
    recover_blue_reference(cfg$archive, cfg$reference, cfg$output, own)
}

