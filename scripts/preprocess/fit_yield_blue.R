# Rscript fit_yield_blue.R --phenotype-2020=FILE --phenotype-2021=FILE
#   --output=NEW_DIRECTORY [--reference=FILE] [--environments=ENV1,ENV2]
# Yield-only export; the optional reference is used for verification, never fitting.
suppressPackageStartupMessages({
    library(data.table)
    library(lme4)
})

setDTthreads(1)

args <- commandArgs(trailingOnly = TRUE)

if (any(!grepl("^--[^=]+=.+$", args))) stop("Use --name=value arguments.")

keys <- sub("^--([^=]+)=.*$", "\\1", args)

values <- sub("^--[^=]+=", "", args)

allowed <- c("phenotype-2020", "phenotype-2021", "output", "reference", "environments")

if (anyDuplicated(keys) || any(!keys %in% allowed)) stop("Duplicate or unknown argument.")

opt <- setNames(as.list(values), keys)

for (key in c("phenotype-2020", "phenotype-2021", "output")) {
    if (is.null(opt[[key]])) 
        stop("Missing --", key)
}

inputs <- unlist(opt[c("phenotype-2020", "phenotype-2021")])

if (any(!file.exists(inputs))) stop("Phenotype input not found.")

if (!is.null(opt$reference) && !file.exists(opt$reference)) stop("Reference not found.")

out <- opt$output

if (dir.exists(out) || file.exists(out)) stop("Output must be a new directory: ", out)

dir.create(out, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(out)) stop("Cannot create output directory.")

writeLines(capture.output(sessionInfo()), file.path(out, "Session_Info.txt"))

writeLines(args, file.path(out, "Arguments.txt"))

sources <- c(inputs, opt$reference)

fwrite(data.table(Path = normalizePath(sources, winslash = "/"), MD5 = unname(tools::md5sum(sources))), 
    file.path(out, "Input_Hashes.csv"))

envs <- c("DEH1.2020", "IAH4.2021", "MIH1.2020", "MNH1.2020", "MNH1.2021", "MOH1.2020", "NEH1.2021", 
    "TXH1.2020", "TXH1.2021", "TXH2.2020", "TXH2.2021", "TXH3.2020", "TXH3.2021", "WIH1.2020", "WIH1.2021", 
    "WIH2.2020", "WIH2.2021", "WIH3.2020", "WIH3.2021")

if (!is.null(opt$environments)) envs <- strsplit(opt$environments, ",", fixed = TRUE)[[1]]

if (anyDuplicated(envs) || any(!nzchar(envs))) stop("Invalid environment list.")

cols <- c("Field-Location", "Year", "Pedigree", "Grain Yield (bu/A)", "Range", "Pass", "Replicate")

d <- rbindlist(lapply(inputs, function(p) fread(p, select = cols)))

d[, `:=`(Env, paste(`Field-Location`, Year, sep = "."))]

d <- d[Env %in% envs]

if (any(!envs %in% d$Env)) stop("Requested environment absent from phenotype inputs.")

numeric_yield <- suppressWarnings(as.numeric(d[["Grain Yield (bu/A)"]]))

bad <- !is.na(d[["Grain Yield (bu/A)"]]) & nzchar(d[["Grain Yield (bu/A)"]]) & is.na(numeric_yield)

if (any(bad)) stop("Non-numeric yield values; resolve explicitly before fitting.")

d[, `:=`(response, numeric_yield * 0.0673)]

factors <- c("Pedigree", "Range", "Pass", "Replicate")

valid <- is.finite(d$response)

for (col in factors) valid <- valid & !is.na(d[[col]]) & nzchar(d[[col]])

d[, `:=`(Usable, valid)]

fwrite(d[, .(InputRows = .N, UsableRows = sum(Usable), ExcludedRows = sum(!Usable)), by = Env], file.path(out, 
    "Row_Accounting.csv"))

results <- components <- diagnostics <- list()

for (env in envs) {
    message("Fitting yield: ", env)
    s <- copy(d[Env == env & Usable])
    for (col in factors) set(s, j = col, value = factor(s[[col]]))
    warnings <- character()
    fit <- tryCatch(withCallingHandlers(lmer(response ~ Pedigree + (1 | Range) + (1 | Pass) + (1 | Replicate), 
        data = s, REML = TRUE, na.action = na.fail), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    }), error = function(e) e)
    if (inherits(fit, "error")) {
        writeLines(paste(env, conditionMessage(fit)), file.path(out, "FAILED.txt"))
        stop("Model failed: ", env, ": ", conditionMessage(fit))
    }
    ids <- levels(model.frame(fit)$Pedigree)
    newdata <- data.frame(Pedigree = factor(ids, levels = ids))
    estimates <- as.numeric(predict(fit, newdata = newdata, re.form = NA))
    results[[env]] <- data.table(Env = env, Pedigree = ids, Yield.t.ha.BLUE = estimates, Pedigree.Env = paste(ids, 
        env, sep = "."))
    vc <- as.data.table(as.data.frame(VarCorr(fit)))
    vc[, `:=`(Env, env)]
    components[[env]] <- vc
    diagnostics[[env]] <- data.table(Env = env, NFit = nobs(fit), Genotypes = length(ids), Singular = isSingular(fit), 
        OptimizerCode = paste(fit@optinfo$conv$opt, collapse = ";"), ConvergenceMessage = paste(fit@optinfo$conv$lme4$messages, 
            collapse = ";"), Warnings = paste(unique(warnings), collapse = ";"))
    fwrite(rbindlist(diagnostics), file.path(out, "Model_Diagnostics.csv"))
}

y <- rbindlist(results)

if (anyDuplicated(y[, .(Env, Pedigree)]) || any(!is.finite(y$Yield.t.ha.BLUE))) stop("Invalid output keys/values.")

setorder(y, Env, Pedigree)

fwrite(y, file.path(out, "Single_Env_Yield_BLUEs_G2F_2020_2021.csv"))

fwrite(rbindlist(components), file.path(out, "Yield_Variance_Components.csv"))

if (!is.null(opt$reference)) {
    ref <- fread(opt$reference, select = c("Env", "Pedigree", "Yield.t.ha.BLUE"))
    ref <- ref[Env %in% envs]
    if (anyDuplicated(ref[, .(Env, Pedigree)])) 
        stop("Duplicate reference keys.")
    fwrite(ref[!is.finite(Yield.t.ha.BLUE)], file.path(out, "Reference_Rows_Without_Yield.csv"))
    ref <- ref[is.finite(Yield.t.ha.BLUE)]
    cmp <- merge(y[, .(Env, Pedigree, Refit = Yield.t.ha.BLUE)], ref, by = c("Env", "Pedigree"), all = TRUE)
    cmp[, `:=`(Difference, Refit - Yield.t.ha.BLUE)]
    qa <- cmp[, .(Keys = .N, MissingComparisons = sum(!is.finite(Refit) | !is.finite(Yield.t.ha.BLUE)), 
        MaxAbsoluteDifference = if (any(is.finite(Difference))) 
            max(abs(Difference), na.rm = TRUE)
        else NA_real_), by = Env]
    qa[, `:=`(Pass, MissingComparisons == 0 & is.finite(MaxAbsoluteDifference) & MaxAbsoluteDifference <= 
        1e-08)]
    fwrite(qa, file.path(out, "Reference_Comparison.csv"))
    fwrite(cmp[!is.finite(Difference) | abs(Difference) > 1e-08], file.path(out, "Reference_Discrepancies.csv"))
    if (any(!qa$Pass) || !setequal(qa$Env, envs)) 
        stop("Reference comparison failed; inspect output.")
}

writeLines(capture.output(sessionInfo()), file.path(out, "Session_Info.txt"))

writeLines(if (is.null(opt$reference)) "COMPLETE (reference not supplied)" else "COMPLETE: reference comparison passed (absolute tolerance 1e-8)", 
    file.path(out, "COMPLETE.txt"))

message("Completed: ", out)

