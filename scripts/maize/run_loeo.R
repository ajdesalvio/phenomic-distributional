# Rscript scripts/maize/run_loeo.R CONFIG.json CACHE NEW_OUTPUT [HELDOUT]
suppressPackageStartupMessages({library(data.table); library(caret); library(glmnet)})
args <- commandArgs(TRUE); stopifnot(length(args)%in%c(3L,4L))
script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)))
engine <- file.path(dirname(script),'../../src/loeo.R'); source(engine)
setDTthreads(1L); Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1)
Sys.setlocale('LC_ALL','C'); foreach::registerDoSEQ()
cfg <- jsonlite::read_json(args[1],simplifyVector=TRUE)
cache <- args[2]; output <- args[3]; workers <- 1L
manifest <- jsonlite::read_json(file.path(cache,'partition_manifest.json'),simplifyVector=TRUE)
stopifnot(identical(manifest$source_sha256,cfg$projection_sha256),
          setequal(manifest$partitions$heldout,cfg$environments),!anyDuplicated(manifest$partitions$heldout))
if (!'sha256'%in%names(manifest$partitions)) stop('Create a new cache with partition_projections.py')
paths <- file.path(cache,manifest$partitions$file)
stopifnot(all(vapply(paths,function(p)digest::digest(file=p,algo='sha256'),character(1))==manifest$partitions$sha256))
yield_path <- cfg$yield$path
if (!grepl('^([A-Za-z]:[/\\\\]|/)',yield_path)) yield_path <- file.path(dirname(args[1]),yield_path)
stopifnot(unname(tools::md5sum(yield_path))==cfg$yield$md5)
yield <- fread(yield_path)[!is.na(Pedigree.Env)&nzchar(Pedigree.Env),.(Pedigree.Env,Yield.t.ha.BLUE)]
stopifnot(!anyDuplicated(yield$Pedigree.Env)); yield[,source_order:=.I]
filters <- list(Median=setdiff(1:100,50L),all=integer())
envs <- cfg$environments
if(length(args)==4L) {stopifnot(args[4]%in%envs);envs<-args[4]}
config <- list(inputs=cfg,partitions=manifest$partitions,
  sources=tools::md5sum(c(script,engine)),R=R.version.string,
  caret=as.character(packageVersion('caret')),glmnet=as.character(packageVersion('glmnet')))
config_hash <- digest::digest(config,algo='sha256')
record <- file.path(output,'run_configuration.rds')
if(file.exists(record)) stopifnot(identical(readRDS(record)$hash,config_hash)) else {
  if(dir.exists(output)&&length(list.files(output,all.files=TRUE,no..=TRUE))) stop('Choose an empty output directory')
  dir.create(output,recursive=TRUE); saveRDS(list(hash=config_hash,configuration=config),record)
}
dir.create(file.path(output,'features'),showWarnings=FALSE)
dir.create(file.path(output,'checkpoints'),showWarnings=FALSE)
for(he in envs) {
  targets <- file.path(output,'checkpoints',paste0(he,'__',names(filters),'__WITHIN.rds'))
  pending <- which(!file.exists(targets)); if(!length(pending)) next
  wide <- make_wide(he)
  for(i in pending) {
    result <- fit_one(wide,he,names(filters)[i],'WITHIN')
    saveRDS(result,paste0(targets[i],'.partial'),compress=FALSE)
    stopifnot(file.rename(paste0(targets[i],'.partial'),targets[i]))
  }
}
results <- lapply(list.files(file.path(output,'checkpoints'),pattern='[.]rds$',full.names=TRUE),readRDS)
stopifnot(all(vapply(results,function(z)identical(z$config_hash,config_hash),logical(1))))
for(n in c('metrics','predictions','coefficients','tuning','folds'))
  fwrite(rbindlist(lapply(results,`[[`,n)),file.path(output,paste0('LOEO_',n,'.csv')))
metrics <- rbindlist(lapply(results,`[[`,'metrics'))
stopifnot(!anyDuplicated(metrics[,.(Heldout,Scenario,Scaling)]))
summary <- metrics[,.(Completed_Heldouts=.N,Mean_Pearson=mean(Pearson),
  Mean_Spearman=mean(Spearman),Mean_RMSE=mean(RMSE),Mean_MAE=mean(MAE)),by=.(Scaling,Scenario)]
fwrite(summary,file.path(output,'LOEO_summary.csv'))
fwrite(data.table(Completed=length(results),Expected=22L),file.path(output,'completion_status.csv'))
writeLines(capture.output(sessionInfo()),file.path(output,'sessionInfo.txt'))
