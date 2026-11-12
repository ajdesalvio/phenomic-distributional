# Rscript this.R --input=RAW_EXGR_CSV --output=NEW_DIR [--batch-size=1200] [--exclusions=CSV]
own <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])

source(file.path(dirname(own), "..", "..", "src", "exgr.R"))

o <- options_cli(c("input", "output"), c("batch-size", "exclusions"))

n <- if (is.null(o[["batch-size"]])) 1200L else as.integer(o[["batch-size"]])

if (is.na(n) || n < 1) stop("Invalid batch size")

input <- normalizePath(o$input, winslash = "/", mustWork = TRUE)

new_directory(o$output)

message("Reading actual ExGR model keys")

d <- fread(input, select = c("Env", "DAP", "probability"))

setnames(d, "probability", "Probability")

ok <- !is.na(d$Env) & nzchar(d$Env) & is.finite(d$DAP) & is.finite(d$Probability)

fwrite(data.table(InputRows = nrow(d), InvalidKeyRows = sum(!ok)), file.path(o$output, "Row_Accounting.csv"))

k <- unique(d[ok])

rm(d)

gc()

if (any(k$Probability < 0 | k$Probability > 1)) stop("Invalid probabilities")

k[, `:=`(Vegetation.Index, "ExGR")]

setorderv(k, model_keys)

k[, `:=`(ModelID = .I, Batch = (.I - 1L)%/%n + 1L, Exclusion = "")]

if (!is.null(o$exclusions)) {
    ex <- fread(o$exclusions)
    valid_keys(ex, c("Env", "DAP"))
    if (!"Reason" %in% names(ex) || anyNA(ex$Reason) || any(!nzchar(ex$Reason))) 
        stop("Exclusions need reasons")
    if (nrow(ex[!unique(k[, .(Env, DAP)]), on = .(Env, DAP)])) 
        stop("Unknown exclusion")
    k[ex, on = .(Env, DAP), `:=`(Exclusion, i.Reason)]
    file.copy(o$exclusions, file.path(o$output, "Exclusions.csv"))
}

batches <- unique(k[, .(Batch)])

batches[, `:=`(File, sprintf("ExGR_Batch_%05d.rds", Batch))]

message("Hashing source input")

p <- list(input = input, input_md5 = hash(input), keys = k, batches = batches, scope = "full raw-data key plan", 
    created = as.character(Sys.time()))

saveRDS(p, file.path(o$output, "ExGR_Plan.rds"))

fwrite(k, file.path(o$output, "Model_Keys.csv"))

fwrite(batches, file.path(o$output, "Batches.csv"))

writeLines(capture.output(sessionInfo()), file.path(o$output, "Session_Info.txt"))

writeLines("COMPLETE", file.path(o$output, "COMPLETE.txt"))

message("Planned ", nrow(k), " keys in ", nrow(batches), " batches")

