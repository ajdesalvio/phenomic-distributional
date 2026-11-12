# Rscript scripts/maize/recalculate_heritability.R CONFIG.json NEW_OUTPUT
suppressPackageStartupMessages(library(data.table));setDTthreads(1L)
args <- commandArgs(TRUE);if(length(args)!=2L)stop('Usage: CONFIG.json NEW_OUTPUT')
script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)),winslash='/')
src <- file.path(dirname(script),'../../src')
source(file.path(src,'maize_heritability.R'));source(file.path(src,'maize_heritability_summary.R'))
cfg <- jsonlite::read_json(args[1],simplifyVector=TRUE)
path <- cfg$components$path
if(!grepl('^([A-Za-z]:[/\\\\]|/)',path))path<-file.path(dirname(args[1]),path)
hash <- function(p)digest::digest(file=p,algo='sha256')
stopifnot(hash(path)==cfg$components$sha256)
out <- args[2];if(file.exists(out))stop('Choose a new output directory')
message('Reading saved components; no model fitting')
d <- fread(path)
stopifnot(is.character(d$vcov),all(c('Env','DAP','Probability')%in%names(d)))
stopifnot(!anyNA(d$vcov),all(grepl('^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$',d$vcov)))
raw_vcov <- d$vcov;d[,vcov:=as.numeric(vcov)]
stopifnot(all(is.finite(d$vcov)),all(d$vcov>=0))
d <- recalculate_within_date_h2(d)
d[,`:=`(vcov_Numeric_For_H2=vcov,H2Individual=H2_Individual_Adjusted)]
d[,vcov:=raw_vcov]
rm(raw_vcov)
keys <- c('Env','DAP','Vegetation.Index','Probability')
d[,PercentUnrounded:=100*vcov_Numeric_For_H2/sum(vcov_Numeric_For_H2),by=keys]
valid <- d[Env%in%cfg$environments & !DAP%in%cfg$excluded_dates]
ped <- valid[grp=='Pedigree']
stopifnot(nrow(d)==5748865L,nrow(ped)==629983L,all(ped$H2_Recalculation_Status=='computed_from_saved_components'))
coverage <- ped[,.(Quantiles=.N,Complete100=identical(sort(as.integer(round(Probability*100))),1:100)),by=.(Env,DAP,Vegetation.Index)]
map <- ped[coverage[Complete100==TRUE,.(Env,DAP,Vegetation.Index)],on=.(Env,DAP,Vegetation.Index),nomatch=0L]
setorderv(map,keys)
tables <- summarize_maize(ped,cfg$environments)
tables$Maize_Selected_11_Environments_All_Valid_Fits_V1 <- ped
tables$Maize_Valid_Curve_Coverage_V1 <- coverage
tables$Maize_MapReady <- map
tables$Maize_Figure_Component_Data_V1 <- valid[Vegetation.Index=='ExGR'|
  (Env=='TXH1.2020'&Vegetation.Index%in%c('GminusR','COM1','EXR')&DAP%in%c(47,68,88,96,111))]
timepath <- cfg$flight_times$path
if(!grepl('^([A-Za-z]:[/\\\\]|/)',timepath))timepath<-file.path(dirname(args[1]),timepath)
stopifnot(hash(timepath)==cfg$flight_times$sha256)
times <- fread(timepath);stopifnot(!anyDuplicated(times[,.(Env,DAP)]))
ftime <- merge(ped[Env=='TXH1.2020'&Vegetation.Index%in%c('ExGR','COM1','GminusR','CIVE')&abs(Probability-.5)<1e-10],
  times[,.(Env,DAP,Flight_Time)],by=c('Env','DAP'),all.x=TRUE)
ftime <- ftime[!is.na(Flight_Time)&nzchar(Flight_Time)]
stopifnot(all(grepl('^[0-9]{1,2}:[0-9]{2}$',ftime$Flight_Time)))
ftime[,TimeHours:=as.numeric(sub(':.*','',Flight_Time))+as.numeric(sub('.*:','',Flight_Time))/60]
tables$Maize_Flight_Time_Data_V1 <- ftime
tables$Maize_Flight_Time_Correlation_V1 <- ftime[,{t<-cor.test(TimeHours,H2Individual);
  .(N=.N,Pearson=unname(t$estimate),P=t$p.value)},by=Vegetation.Index]
stopifnot(hash(path)==cfg$components$sha256)
dir.create(out,recursive=TRUE)
for(n in names(tables))fwrite(tables[[n]],file.path(out,paste0(n,'.csv')))
jsonlite::write_json(list(status='COMPLETE_SAVED_COMPONENT_RECALCULATION',config=cfg,
  formula='Vg/(Vg+Ve); no residual replication divisor',model_fitting=FALSE,
  historical_convergence_not_recertified=TRUE,canonical_indices=31L,
  excluded_redundant_labels=c('GCC','NDI'),source_md5=as.list(tools::md5sum(c(script,
  file.path(src,c('maize_heritability.R','maize_heritability_summary.R')))))),
  file.path(out,'provenance.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
