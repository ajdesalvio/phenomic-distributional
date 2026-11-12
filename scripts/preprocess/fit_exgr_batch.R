# Resumable worker: --plan=FILE --batch=N --output=DIR. Checkpoints every 25 models.
own <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])

source(file.path(dirname(own), "..", "..", "src", "exgr.R"))

suppressPackageStartupMessages(library(lme4))

setDTthreads(1)

o <- options_cli(c("plan", "batch", "output"))

p <- load_plan(o$plan)

id <- as.integer(o$batch)

plan_hash <- hash(o$plan)

k <- copy(p$keys[Batch == id])

if (!nrow(k)) stop("Unknown batch")

k[, `:=`(Piece, (seq_len(.N) - 1L)%/%25L + 1L)]

dir.create(o$output, recursive = TRUE, showWarnings = FALSE)

dest <- file.path(o$output, p$batches[Batch == id, File])

check <- function(x, expected) {
    if (!identical(x$plan_md5, plan_hash) || x$batch != id || x$source != "raw_model_refit") 
        stop("Wrong checkpoint provenance")
    valid_keys(x$status, model_keys)
    if (nrow(x$status[!expected, on = model_keys]) || nrow(expected[!x$status, on = model_keys])) 
        stop("Checkpoint key mismatch")
    if (any(x$status$Status == "failed")) 
        stop("Checkpoint contains failed models; inspect diagnostics")
}

if (file.exists(dest)) {
    check(readRDS(dest), k)
    message("Already complete: batch ", id)
    quit(status = 0L)
}

journal <- file.path(o$output, sprintf("Checkpoints_%05d", id))

dir.create(journal, showWarnings = FALSE)

piece_ids <- unique(k$Piece)

piece_files <- file.path(journal, sprintf("Piece_%04d.rds", piece_ids))

for (j in which(file.exists(piece_files))) check(readRDS(piece_files[j]), k[Piece == piece_ids[j]])

pending <- piece_ids[!file.exists(piece_files)]

if (length(pending)) {
    if (!identical(hash(p$input), p$input_md5)) 
        stop("Raw input hash changed")
    d <- fread(p$input, select = c("Env", "DAP", "probability", "quantile", "Pedigree", "Range", "Pass", 
        "Replicate"))
    setnames(d, "probability", "Probability")
    d <- d[k[Piece %in% pending], on = .(Env, DAP, Probability), nomatch = 0L]
    for (piece in pending) {
        selected <- k[Piece == piece]
        g <- vc <- ledger <- list()
        for (i in seq_len(nrow(selected))) {
            z <- selected[i]
            started <- proc.time()[["elapsed"]]
            status <- copy(z[, ..model_keys])
            status[, `:=`(Status = "excluded", Detail = z$Exclusion, NFit = 0L, Genotypes = 0L, Singular = NA, 
                OptimizerCode = "", ElapsedSeconds = 0)]
            if (!nzchar(z$Exclusion)) {
                s <- copy(d[ModelID == z$ModelID])
                for (col in c("Pedigree", "Range", "Pass", "Replicate")) set(s, j = col, value = factor(s[[col]]))
                s[!is.finite(quantile), `:=`(quantile, NA_real_)]
                warnings <- character()
                message("Fitting ", z$Env, " / ", z$DAP, " / ", z$Probability)
                fit <- tryCatch(withCallingHandlers(lmer(quantile ~ Pedigree + (1 | Range) + (1 | Pass) + 
                  (1 | Replicate), data = s, REML = TRUE, na.action = na.omit), warning = function(w) {
                  warnings <<- c(warnings, conditionMessage(w))
                  invokeRestart("muffleWarning")
                }), error = identity)
                if (inherits(fit, "error")) {
                  status[, `:=`(Status = "failed", Detail = conditionMessage(fit))]
                }
                else {
                  ids <- levels(model.frame(fit)$Pedigree)
                  means <- predict(fit, newdata = data.frame(Pedigree = factor(ids, levels = ids)), re.form = NA)
                  g[[i]] <- cbind(z[rep(1L, length(ids)), ..model_keys], data.table(Pedigree = ids, Quantile.BLUE = as.numeric(means)))
                  v <- as.data.table(as.data.frame(VarCorr(fit)))
                  vc[[i]] <- cbind(z[rep(1L, nrow(v)), ..model_keys], v)
                  status[, `:=`(Status = "fitted", NFit = nobs(fit), Genotypes = length(ids), Singular = isSingular(fit), 
                    OptimizerCode = paste(fit@optinfo$conv$opt, collapse = ";"), Detail = paste(c(warnings, 
                      fit@optinfo$conv$lme4$messages), collapse = ";"))]
                }
            }
            status[, `:=`(ElapsedSeconds, proc.time()[["elapsed"]] - started)]
            ledger[[i]] <- status
        }
        x <- list(plan_md5 = plan_hash, batch = id, source = "raw_model_refit", blue = rbindlist(g), 
            components = rbindlist(vc), status = rbindlist(ledger), session = capture.output(sessionInfo()), 
            worker_md5 = hash(own), helper_md5 = hash(file.path(dirname(own), "..", "..", "src", "exgr.R")))
        target <- piece_files[match(piece, piece_ids)]
        tmp <- paste0(target, ".", Sys.getpid(), ".partial")
        saveRDS(x, tmp)
        if (file.exists(target) || !file.rename(tmp, target)) 
            stop("Cannot finalize checkpoint")
        fwrite(x$status, sub("[.]rds$", "_Status.csv", target))
        if (any(x$status$Status == "failed")) 
            stop("Failed models recorded; resolve before resume")
    }
    rm(d)
    gc()
}

parts <- lapply(piece_files, readRDS)

for (j in seq_along(parts)) check(parts[[j]], k[Piece == piece_ids[j]])

x <- list(plan_md5 = plan_hash, batch = id, source = "raw_model_refit", blue = rbindlist(lapply(parts, 
    `[[`, "blue"), fill = TRUE), components = rbindlist(lapply(parts, `[[`, "components"), fill = TRUE), 
    status = rbindlist(lapply(parts, `[[`, "status")), session = capture.output(sessionInfo()), worker_md5 = hash(own))

tmp <- paste0(dest, ".", Sys.getpid(), ".partial")

saveRDS(x, tmp)

if (file.exists(dest) || !file.rename(tmp, dest)) stop("Cannot finalize batch")

message("Completed batch ", id)

