# Join image quantiles to a unique field map, or verify an existing joined export.
# Rscript this.R --input CSV --field-map CSV --output NEW_DIR
#   [--mode join|verify] [--allow-unmatched true|false]
join_field_info <- function(raw, field) {
    raw <- data.table::copy(data.table::as.data.table(raw))
    field <- data.table::copy(data.table::as.data.table(field))
    measures <- c("file_name", "vi_name", "probability", "quantile")
    key <- "File.Name.TIF.Unified.Masked"
    fields <- c("Env", "Pedigree", "DAP", "Range", "Pass", "Replicate")
    if (!all(measures %in% names(raw)) || !all(c(key, fields) %in% names(field))) 
        stop("Missing measurement or field-map columns")
    for (n in c("file_name", "vi_name")) if (anyNA(raw[[n]]) || any(!nzchar(trimws(raw[[n]])))) 
        stop("Empty raw identifier: ", n)
    if (anyNA(field[[key]]) || any(!nzchar(trimws(field[[key]]))) || anyDuplicated(field[[key]])) 
        stop("Field-map image keys must be nonempty and unique")
    if (!is.numeric(raw$probability) || any(!is.finite(raw$probability)) || any(raw$probability < 0.01 | 
        raw$probability > 1 | abs(raw$probability * 100 - round(raw$probability * 100)) > 1e-07)) 
        stop("Invalid quantile levels")
    if (!is.numeric(raw$quantile)) 
        stop("Nonnumeric quantile values")
    if (anyDuplicated(raw[, .(file_name, vi_name, probability)])) 
        stop("Duplicate image/VI/quantile rows")
    if (any(grepl("[.]tiff?$", raw$file_name, ignore.case = TRUE))) 
        stop("Expected extensionless file_name")
    images <- unique(raw[, .(file_name)])
    images[, `:=`((key), paste0(file_name, ".tif"))]
    ix <- match(images[[key]], field[[key]])
    unmatched <- images[is.na(ix)]
    assignments <- cbind(images, field[ix, ..fields])
    ri <- match(raw$file_name, assignments$file_name)
    joined <- raw[, ..measures]
    joined[, `:=`((key), assignments[[key]][ri])]
    for (n in fields) joined[, `:=`((n), assignments[[n]][ri])]
    list(data = joined, assignments = assignments, unmatched = unmatched)
}

main <- function(args = commandArgs(TRUE)) {
    suppressPackageStartupMessages(library(data.table))
    setDTthreads(2)
    cfg <- list(input = NULL, `field-map` = NULL, output = NULL, mode = "join", `allow-unmatched` = "false")
    if (!length(args) || length(args)%%2) 
        stop("Supply --name value pairs")
    seen <- character()
    for (i in seq(1, length(args), 2)) {
        n <- sub("^--", "", args[i])
        if (!n %in% names(cfg) || n %in% seen) 
            stop("Unknown/repeated option: ", n)
        cfg[[n]] <- args[i + 1]
        seen <- c(seen, n)
    }
    if (any(vapply(cfg[c("input", "field-map", "output")], is.null, logical(1)))) 
        stop("Input, field-map and output are required")
    if (!cfg$mode %in% c("join", "verify") || !cfg[["allow-unmatched"]] %in% c("true", "false")) 
        stop("Invalid mode/options")
    if (file.exists(cfg$output) || dir.exists(cfg$output)) 
        stop("Output exists; choose a new numbered directory")
    own <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
    inputs <- normalizePath(c(cfg$input, cfg[["field-map"]], own), winslash = "/", mustWork = TRUE)
    before <- tools::md5sum(inputs)
    fields <- c("Env", "Pedigree", "DAP", "Range", "Pass", "Replicate")
    cols <- c("file_name", "vi_name", "probability", "quantile", if (cfg$mode == "verify") fields)
    message("Reading measurements and rebuilding field assignments")
    raw <- fread(inputs[1], select = cols, showProgress = FALSE)
    field <- fread(inputs[2], select = c("File.Name.TIF.Unified.Masked", fields), showProgress = FALSE)
    ans <- join_field_info(raw, field)
    dir.create(cfg$output, recursive = TRUE)
    fwrite(ans$unmatched, file.path(cfg$output, "Unmatched_Images.csv"))
    fwrite(ans$assignments, file.path(cfg$output, "Image_Field_Assignments.csv"))
    checks <- NULL
    if (cfg$mode == "verify") {
        checks <- rbindlist(lapply(fields, function(n) {
            x <- as.character(raw[[n]])
            y <- as.character(ans$data[[n]])
            x[!is.na(x) & x == ""] <- NA_character_
            y[!is.na(y) & y == ""] <- NA_character_
            data.table(Field = n, MismatchedRows = sum(is.na(x) != is.na(y) | (!is.na(x) & !is.na(y) & 
                x != y)))
        }))
        fwrite(checks, file.path(cfg$output, "Field_Comparison.csv"))
    }
    matched <- nrow(ans$unmatched) == 0 || cfg[["allow-unmatched"]] == "true"
    equal <- is.null(checks) || all(checks$MismatchedRows == 0)
    stable <- identical(before, tools::md5sum(inputs))
    passed <- matched && equal && stable
    if (passed && cfg$mode == "join") 
        fwrite(ans$data, file.path(cfg$output, "Quantiles_FieldInfo.csv"))
    outputs <- list.files(cfg$output, full.names = TRUE)
    receipt <- list(status = if (passed) "PASS" else "FAIL", mode = cfg$mode, rows = nrow(raw), images = nrow(ans$assignments), 
        unmatched_images = nrow(ans$unmatched), allow_unmatched = cfg[["allow-unmatched"]] == "true", 
        comparison_rule = "Exact populated fields; empty strings and NA both denote missing fields", 
        field_comparison = checks, input_md5 = as.list(before), inputs_unchanged = stable, output_md5 = as.list(tools::md5sum(outputs)), 
        session = capture.output(sessionInfo()), scope = "Core analysis field assignments only. Verify mode reuses saved image quantiles; it does not validate image extraction or masking.")
    jsonlite::write_json(receipt, file.path(cfg$output, "Field_Join_Manifest.json"), pretty = TRUE, auto_unbox = TRUE)
    if (!passed) 
        stop("Field join failed; see diagnostics. No joined quantile CSV published.")
    message("PASS: ", nrow(raw), " rows; ", nrow(ans$assignments), " images")
}

if (sys.nframe() == 0L) main()

