# --rebuilt=CSV --reference=CSV --output=NEW_DIR --status=MODEL_STATUS_CSV
own <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])

source(file.path(dirname(own), "..", "..", "src", "exgr.R"))

o <- options_cli(c("rebuilt", "reference", "output", "status"))

new_directory(o$output)

x <- fread(o$rebuilt, select = blue_columns)

y <- fread(o$reference, select = blue_columns)

valid_keys(x, blue_keys)

valid_keys(y, blue_keys)

missing <- y[!x, on = blue_keys]

extra <- x[!y, on = blue_keys]

fwrite(missing, file.path(o$output, "Missing_BLUE_Keys.csv"))

fwrite(extra, file.path(o$output, "Extra_BLUE_Keys.csv"))

if (nrow(missing) || nrow(extra)) {
    writeLines("FAIL: different genotype/model keys", file.path(o$output, "RESULT.txt"))
    stop("BLUE key mismatch")
}

setorderv(x, blue_keys)

setorderv(y, blue_keys)

x[, `:=`(Reference = y$Quantile.BLUE, Difference = Quantile.BLUE - y$Quantile.BLUE)]

tolerance <- 1e-07

bad <- !is.finite(x$Difference) | abs(x$Difference) > tolerance

fwrite(x[bad], file.path(o$output, "BLUE_Discrepancies.csv"))

qa <- x[, .(Rows = .N, Nonfinite = sum(!is.finite(Difference)), AboveTolerance = sum(!is.finite(Difference) | 
    abs(Difference) > tolerance), MaxAbsoluteDifference = max(abs(Difference))), by = model_keys]

fwrite(qa, file.path(o$output, "Comparison_By_Model.csv"))

summary <- data.table(Rows = nrow(x), Models = nrow(qa), Tolerance = tolerance, MaxAbsoluteDifference = max(abs(x$Difference)), 
    RowsAboveTolerance = sum(bad), ModelsAboveTolerance = sum(qa$AboveTolerance > 0), Pass = !any(bad))

fwrite(summary, file.path(o$output, "Comparison_Summary.csv"))

status <- fread(o$status)

fwrite(status[, .(Models = .N), by = .(Status, Singular, OptimizerCode)], file.path(o$output, "Fit_Diagnostic_Summary.csv"))

fwrite(data.table(Path = c(o$rebuilt, o$reference, o$status), MD5 = vapply(c(o$rebuilt, o$reference, 
    o$status), hash, character(1))), file.path(o$output, "Compared_File_Hashes.csv"))

writeLines(if (summary$Pass) "PASS: exact keys; BLUE values within 1e-7 VI units" else "FAIL: discrepancies saved; rebuilt outputs preserved", 
    file.path(o$output, "RESULT.txt"))

print(summary)

if (!summary$Pass) stop("Numerical comparison failed; inspect saved discrepancies")

