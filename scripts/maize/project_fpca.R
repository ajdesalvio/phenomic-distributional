# Rscript scripts/maize/project_fpca.R CONFIG.json NEW.csv [HELDOUT] [QUANTILE]
suppressPackageStartupMessages({library(data.table);library(dplyr);library(tidyr);library(fdapace)})
args <- commandArgs(TRUE); stopifnot(length(args)%in%2:4)
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))
source(file.path(dirname(script),'../../src/fpca_projection.R'))
cfg <- jsonlite::read_json(args[1],simplifyVector=TRUE)
path <- cfg$blue$path
if(!grepl('^([A-Za-z]:[/\\\\]|/)',path)) path <- file.path(dirname(args[1]),path)
stopifnot(unname(tools::md5sum(path))==cfg$blue$md5)
if(any(file.exists(paste0(args[2],c('','.diagnostics.pdf','.sessionInfo.txt'))))) stop('Choose a new output file')
heldouts <- if(length(args)>=3L)args[3] else cfg$environments
quantiles <- if(length(args)==4L)as.numeric(args[4]) else (1:100)/100
stopifnot(all(heldouts%in%cfg$environments),all(quantiles%in%((1:100)/100)))
VI <- as.data.frame(fread(path))
stopifnot(!anyDuplicated(VI[,c('Pedigree','Env','DAP','Vegetation.Index','Probability')]))
dir.create(dirname(args[2]),recursive=TRUE,showWarnings=FALSE)
pdf(paste0(args[2],'.diagnostics.pdf'))
result <- tryCatch(project_fpca(VI,cfg$environments,heldouts,quantiles),finally=dev.off())
stopifnot(uniqueN(result[,.(Heldout.Env,Vegetation.Index,Quantile.Level)])==length(heldouts)*length(quantiles)*length(unique(VI$Vegetation.Index)))
fwrite(result,args[2])
writeLines(capture.output(sessionInfo()),paste0(args[2],'.sessionInfo.txt'))
