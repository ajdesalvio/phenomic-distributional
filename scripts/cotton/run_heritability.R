# Rscript scripts/cotton/run_heritability.R CONFIG.json NEW_OUTPUT smoke|full
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(TRUE)
if(length(args)!=3L||!args[3]%in%c('smoke','full')) stop('Usage: CONFIG.json NEW_OUTPUT smoke|full')
script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)))
src <- file.path(dirname(script),'../../src')
source(file.path(src,'cotton_io.R')); source(file.path(src,'cotton_heritability.R'))
source(file.path(src,'cotton_summary.R'))
cfg <- jsonlite::read_json(args[1],simplifyVector=TRUE)
stopifnot(identical(names(cfg$inputs),c('PSRI','CVI','ExG2')))
paths <- vapply(cfg$inputs,cotton_input,character(1),config=args[1])
field <- fread(cotton_input(cfg$fieldmap,args[1]),data.table=FALSE)
stopifnot(nrow(field)==480L,!anyDuplicated(field$Pltg_ID))
out <- args[2]; if(file.exists(out)) stop('Choose a new output directory')
dir.create(out,recursive=TRUE);dir.create(file.path(out,'input_cache'));dir.create(file.path(out,'checkpoints'))
setDTthreads(1L);Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1)
smoke <- args[3]=='smoke'; metrics <- components <- list(); common <- NULL
for(vi in names(paths)) {
  message('Validating ',vi)
  prepared <- v6_input(paths[[vi]],vi,field)
  keys <- prepared$data[,c('DAT','Pltg_ID','Pedigree','Block')]
  if(is.null(common)) common <- keys else stopifnot(identical(common,keys))
  cache <- file.path(out,'input_cache',paste0(vi,'_Unsmoothed_Quantiles_V6.rds'))
  v6_save(prepared,cache)
  for(day in if(smoke)57L else prepared$dates) {
    task <- list(VI=vi,DAT=day,qids=if(smoke)c(1L,50L,100L)else 1:100,
      cache=cache,output=file.path(out,'checkpoints',sprintf('%s_DAT%03d_V6.rds',vi,day)))
    v6_fit_date(task); z <- readRDS(task$output)
    metrics[[length(metrics)+1L]] <- z$metrics;components[[length(components)+1L]] <- z$components
  }
}
m <- rbindlist(metrics);v <- rbindlist(components)
stopifnot(nrow(m)==if(smoke)9L else 30300L,all(m$Converged),
  !anyDuplicated(m[,.(Vegetation.Index,DAT,QuantileID)]),all(is.finite(m$H2Individual)),
  max(abs(m$H2Individual-m$Vg/(m$Vg+m$Ve)))<1e-14)
curves <- m[,.(Quantiles=.N,Complete100=identical(sort(as.integer(QuantileID)),1:100),
  AllConverged=all(Converged),SingularQuantiles=sum(Singular),MinH2=min(H2Individual),
  MedianH2=median(H2Individual),MaxH2=max(H2Individual)),by=.(Vegetation.Index,DAT)]
if(!smoke) stopifnot(nrow(curves)==303L,all(curves$Complete100))
fwrite(m,file.path(out,'Fit_Metrics_V6.csv'));fwrite(v,file.path(out,'Full_Components_V6.csv'))
fwrite(curves,file.path(out,'Curve_Index_V6.csv'));fwrite(common,file.path(out,'Shared_Cohort_Keys_V6.csv'))
if(!smoke) {
  fwrite(m,file.path(out,'MapReady_Complete_Curves_V6.csv'))
  tables <- cotton_summary(m)
  for(n in names(tables)) fwrite(tables[[n]],file.path(out,paste0(n,'.csv')))
}
sources <- c(script,file.path(src,c('cotton_io.R','cotton_heritability.R','cotton_summary.R')))
jsonlite::write_json(list(status=if(smoke)'SMOKE_ONLY'else'COMPLETE',fits=nrow(m),
  config=cfg,formula=block_only_formula_v5(),source_md5=as.list(tools::md5sum(sources))),
  file.path(out,'run_manifest.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
