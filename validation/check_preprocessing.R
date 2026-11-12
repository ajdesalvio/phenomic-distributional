# Rscript validation/check_preprocessing.R NEW_TEST_DIRECTORY
suppressPackageStartupMessages({
    library(data.table)
    library(jsonlite)
})

own <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])

root <- normalizePath(file.path(dirname(own), ".."), winslash = "/")

out <- commandArgs(TRUE)

if (length(out) != 1L || file.exists(out)) stop("Supply a new test directory")

dir.create(out, recursive = TRUE)

out <- normalizePath(out, winslash = "/")

scripts <- file.path(root, "scripts/preprocess")

checks <- list()

check <- function(name, value) {
    checks[[name]] <<- data.table(Check = name, Pass = isTRUE(value))
    if (!isTRUE(value)) 
        stop(name)
}

fails <- function(expr) inherits(tryCatch({
    force(expr)
    NULL
}, error = identity), "error")

for (p in list.files(root, pattern = "[.]R$", recursive = TRUE, full.names = TRUE)) parse(p)

check("All R scripts parse", TRUE)

source(file.path(scripts, "join_field_info.R"))

source(file.path(scripts, "build_raw_support.R"))

source(file.path(scripts, "validate_inputs.R"))

raw <- data.table(file_name = c("a", "b", "c"), vi_name = "ExGR", probability = 0.5, quantile = c(1, 
    3, NA_real_))

map <- data.table(File.Name.TIF.Unified.Masked = c("a.tif", "b.tif", "c.tif"), Env = "E", Pedigree = c("G", 
    "G", "H"), DAP = 1L, Range = 1L, Pass = 1L, Replicate = 1L)

j <- join_field_info(raw, map)

check("Field join preserves measurement rows", nrow(j$data) == 3L && identical(j$data$quantile, raw$quantile))

check("Duplicate field keys rejected", fails(join_field_info(raw, rbind(map, map[1]))))

check("Unmatched image recorded", nrow(join_field_info(raw, map[-1])$unmatched) == 1L)

blue <- unique(j$data[, .(Env, Pedigree, DAP, Probability = probability, Vegetation.Index = vi_name)])

support <- support_build(j$data, blue, envs = "E")

check("Finite mean and counts", support[Pedigree == "G", Quantile.Raw.Mean] == 2 && support[Pedigree == 
    "G", N.Raw.Finite] == 2)

check("All missing is retained with count zero", is.na(support[Pedigree == "H", Quantile.Raw.Mean]) && 
    support[Pedigree == "H", N.Raw.Finite] == 0)

check("Duplicate measurement rejected", fails(support_build(rbind(j$data, j$data[1]), blue, envs = "E")))

check("Missing raw support rejected", fails(support_build(j$data[Pedigree == "G"], blue, envs = "E")))

run <- function(script, args, tag, expect = 0L) {
    code <- system2(file.path(R.home("bin"), if (.Platform$OS.type == "windows") 
        "Rscript.exe"
    else "Rscript"), c("--vanilla", shQuote(file.path(scripts, script)), shQuote(args)), stdout = file.path(out, 
        paste0(tag, ".log")), stderr = file.path(out, paste0(tag, ".log")))
    check(tag, if (expect == 0L) 
        code == 0L
    else code != 0L)
}

e <- new.env()

e$blue.filt <- data.table(Env = "E", Pedigree = c("G", "H"), Probability = 0.5, Quantile.BLUE = c(1, 
    2), DAP = 1L, Vegetation.Index = "ExGR")

save(list = "blue.filt", envir = e, file = file.path(out, "archive.RData"))

fwrite(e$blue.filt, file.path(out, "reference.csv"))

run("recover_blue_reference.R", c("--archive", file.path(out, "archive.RData"), "--reference", file.path(out, 
    "reference.csv"), "--output", file.path(out, "recovery")), "Archive recovery")

check("Recovered bytes identical", identical(unname(tools::md5sum(file.path(out, "reference.csv"))), 
    unname(tools::md5sum(file.path(out, "recovery/Quantile_BLUEs_G2F_Distributional_ExGR_Reference.csv")))))

run("recover_blue_reference.R", c("--archive", file.path(out, "archive.RData"), "--reference", file.path(out, 
    "reference.csv"), "--output", file.path(out, "recovery")), "Overwrite rejected", 1L)

roles <- c("raw", "blue", "support", "yield", "field_map")

inputs <- setNames(lapply(roles, function(role) {
    p <- paste0(role, ".txt")
    writeLines(role, file.path(out, p))
    list(path = p, md5 = unname(tools::md5sum(file.path(out, p))))
}), roles)

builder_hash <- "11111111111111111111111111111111"

prov <- list(source_md5 = c(old_raw = inputs$raw$md5, old_blue = inputs$blue$md5, old_builder = builder_hash), 
    environments = "E")

saveRDS(prov, file.path(out, "support.rds"))

inputs$support_provenance <- list(path = "support.rds", md5 = unname(tools::md5sum(file.path(out, "support.rds"))))

cfg <- list(inputs = inputs, support_builder_md5 = builder_hash, environments = "E")

write_json(cfg, file.path(out, "config.json"), auto_unbox = TRUE)

validate_inputs(file.path(out, "config.json"), file.path(out, "receipt.json"))

check("Relative paths and relocated provenance", TRUE)

writeLines("changed", file.path(out, "raw.txt"))

check("Changed input rejected", fails(validate_inputs(file.path(out, "config.json"), file.path(out, "invalid.json"))))

rawcsv <- file.path(out, "raw.csv")

fwrite(j$data, rawcsv)

run("plan_exgr.R", c(paste0("--input=", rawcsv), paste0("--output=", file.path(out, "plan")), "--batch-size=1"), 
    "Explicit batch planning")

plan <- readRDS(file.path(out, "plan/ExGR_Plan.rds"))

check("Plan uses actual model keys", nrow(plan$keys) == 1L)

fwrite(rbindlist(checks), file.path(out, "Checks.csv"))

versions <- data.table(Package = c("R", "data.table", "jsonlite", "lme4"), Version = c(as.character(getRversion()), 
    vapply(c("data.table", "jsonlite", "lme4"), function(x) as.character(packageVersion(x)), character(1))))

fwrite(versions, file.path(out, "R_versions.csv"))

print(rbindlist(checks))

