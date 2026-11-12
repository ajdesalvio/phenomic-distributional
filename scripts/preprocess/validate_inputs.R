# Rscript validate_inputs.R CONFIG.json NEW_RECEIPT.json
# Relative input paths resolve from the configuration directory.
validate_inputs <- function(config, output) {
    if (file.exists(output)) 
        stop("Receipt already exists")
    cfg <- jsonlite::fromJSON(config, simplifyVector = FALSE)
    roles <- c("raw", "blue", "support", "yield", "field_map", "support_provenance")
    if (!setequal(names(cfg$inputs), roles)) 
        stop("Unexpected input roles")
    paths <- vapply(cfg$inputs, function(x) x$path, character(1))
    relative <- !grepl("^([A-Za-z]:[/\\]|/|\\\\)", paths)
    paths[relative] <- file.path(dirname(normalizePath(config)), paths[relative])
    paths <- vapply(paths, normalizePath, character(1), winslash = "/", mustWork = TRUE)
    expected <- vapply(cfg$inputs, function(x) x$md5, character(1))
    actual <- tools::md5sum(paths)
    if (!identical(tolower(unname(actual)), tolower(unname(expected)))) 
        stop("Input hash mismatch")
    support <- readRDS(paths[["support_provenance"]])
    required <- c(expected[c("raw", "blue")], cfg$support_builder_md5)
    if (length(cfg$support_builder_md5) != 1L || !all(tolower(required) %in% tolower(unname(support$source_md5)))) 
        stop("Support provenance mismatch")
    if (!setequal(unlist(cfg$environments), support$environments)) 
        stop("Support environment mismatch")
    receipt <- list(status = "PASS", inputs = as.list(actual), environments = cfg$environments, support_provenance_matches = TRUE, 
        scope = "Frozen manuscript inputs; no model refitting.")
    dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(receipt, output, pretty = TRUE, auto_unbox = TRUE)
}

if (sys.nframe() == 0L) {
    a <- commandArgs(TRUE)
    if (length(a) != 2L) 
        stop("Supply CONFIG.json NEW_RECEIPT.json")
    validate_inputs(a[1], a[2])
}

