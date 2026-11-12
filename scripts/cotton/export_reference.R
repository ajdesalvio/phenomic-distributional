# Rscript scripts/cotton/export_reference.R CONFIG.json NEW_OUTPUT
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(TRUE);if(length(args)!=2L)stop('Usage: CONFIG.json NEW_OUTPUT')
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))
src <- file.path(dirname(script),'../../src')
for(n in c('cotton_io.R','cotton_summary.R','cotton_functional.R')) source(file.path(src,n))
cfg <- jsonlite::read_json(args[1],simplifyVector=TRUE)
paths <- vapply(cfg$reference,cotton_input,character(1),config=args[1])
out <- args[2];if(file.exists(out))stop('Choose a new output directory')
d <- fread(paths[['heritability']])
stopifnot(nrow(d)==30300L,all(d$Converged),uniqueN(d$DAT)==101L,
  setequal(d$Vegetation.Index,c('PSRI','CVI','ExG2')),
  !anyDuplicated(d[,.(Vegetation.Index,DAT,QuantileID)]),
  max(abs(d$H2Individual-d$Vg/(d$Vg+d$Ve)))<1e-12)
tables <- cotton_summary(d)
for(vi in c('PSRI','CVI','ExG2')) {
  model <- readRDS(paths[[vi]])
  stopifnot(model$provenance$configuration$vi==vi,model$bootstrap$draws==4999L)
  z <- v7_tables(model$point,model$bootstrap,model$provenance$configuration)
  stopifnot(nrow(z$curves)==1500L,nrow(z$means)==1600L)
  tables[[paste0(vi,'_genotype_contrasts_V7')]] <- z$curves
  tables[[paste0(vi,'_genotype_adjusted_means_point_only_V7')]] <- z$means
}
dir.create(out,recursive=TRUE)
for(n in names(tables)) fwrite(tables[[n]],file.path(out,paste0(n,'.csv')))
jsonlite::write_json(list(status='EXPORTED_SAVED_RESULTS_NO_REFIT',inputs=cfg$reference),
  file.path(out,'provenance.json'),pretty=TRUE,auto_unbox=TRUE)
