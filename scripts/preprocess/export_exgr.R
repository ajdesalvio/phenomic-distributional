# Rscript this.R --plan=PLAN_RDS --chunks=DIR --output=NEW_DIR [--reference=CSV]
own <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])

source(file.path(dirname(own), "..", "..", "src", "exgr.R"))

o <- options_cli(c("plan", "chunks", "output"), "reference")

p <- load_plan(o$plan)

files <- file.path(o$chunks, p$batches$File)

if (any(!file.exists(files))) stop("Missing declared chunks: ", paste(basename(files[!file.exists(files)]), 
    collapse = ", "))

new_directory(o$output)

g <- v <- logs <- list()

provenance <- character()

for (i in seq_along(files)) {
    x <- readRDS(files[i])
    if (!identical(x$source, "raw_model_refit")) 
        stop("Full verification requires raw model refits")
    expected <- p$keys[Batch == p$batches$Batch[i]]
    if (!identical(x$plan_md5, hash(o$plan)) || x$batch != p$batches$Batch[i]) 
        stop("Wrong plan/batch in chunk")
    valid_keys(x$status, model_keys)
    if (nrow(x$status[!expected, on = model_keys]) || nrow(expected[!x$status, on = model_keys])) 
        stop("Chunk status coverage mismatch")
    if (any(!x$status$Status %in% c("fitted", "archived", "excluded"))) 
        stop("Failed/unknown status")
    exclusions <- expected[nzchar(Exclusion)]
    skipped <- x$status[Status == "excluded"]
    if (nrow(skipped[!exclusions, on = model_keys]) || nrow(exclusions[!skipped, on = model_keys])) 
        stop("Unapproved/missing exclusion")
    if (nrow(x$blue)) {
        valid_keys(x$blue, blue_keys)
        counts <- x$blue[, .(Actual = .N), by = model_keys]
    }
    else counts <- expected[0, c(model_keys), with = FALSE][, `:=`(Actual, integer())]
    successful <- x$status[Status != "excluded"]
    if (nrow(counts[!successful, on = model_keys]) || nrow(successful[!counts, on = model_keys])) 
        stop("BLUE coverage mismatch")
    if (any(counts[successful, on = model_keys]$Actual != successful$Genotypes)) 
        stop("Genotype count mismatch")
    if (nrow(successful)) {
        valid_keys(x$components, c(model_keys, "grp"))
        vc_keys <- unique(x$components[, ..model_keys])
        if (nrow(vc_keys[!successful, on = model_keys]) || nrow(successful[!vc_keys, on = model_keys])) 
            stop("Variance-component coverage mismatch")
        groups <- x$components[, .(Valid = setequal(grp, c("Pass", "Range", "Replicate", "Residual"))), 
            by = model_keys]
        if (any(!groups$Valid)) 
            stop("Unexpected variance-component groups")
    }
    else if (nrow(x$components)) 
        stop("Components for excluded batch")
    if (nrow(x$blue) && any(!is.finite(x$blue$Quantile.BLUE))) 
        stop("Nonfinite BLUE")
    g[[i]] <- x$blue
    v[[i]] <- x$components
    logs[[i]] <- x$status
    provenance <- c(provenance, x$source)
}

y <- rbindlist(g)[, ..blue_columns]

valid_keys(y, blue_keys)

setorderv(y, blue_keys)

fwrite(y, file.path(o$output, "Quantile_BLUEs_G2F_Distributional_ExGR.csv"))

fwrite(rbindlist(v, fill = TRUE), file.path(o$output, "ExGR_Variance_Components.csv"))

fwrite(rbindlist(logs), file.path(o$output, "Model_Status.csv"))

fwrite(data.table(File = files, MD5 = vapply(files, hash, character(1)), Source = provenance), file.path(o$output, 
    "Chunk_Provenance.csv"))

writeLines(c(paste("Plan MD5:", hash(o$plan)), paste("Scope:", p$scope), paste("Sources:", paste(unique(provenance), 
    collapse = ", ")), capture.output(sessionInfo())), file.path(o$output, "Provenance.txt"))

if (!is.null(o$reference)) compare_blue(y, o$reference, o$output)

writeLines("COMPLETE", file.path(o$output, "COMPLETE.txt"))

