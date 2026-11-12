# Rscript scripts/figures/render_results.R CONFIG.json NEW_OUTPUT [all|3|4|6|7|8]
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(cowplot)})
setDTthreads(1L)
args<-commandArgs(TRUE)
if(!length(args)%in%2:3)stop('Usage: CONFIG.json NEW_OUTPUT [all|3|4|6|7|8]')
which<-if(length(args)==3L)args[3]else'all'
stopifnot(which%in%c('all','3','4','6','7','8'))
if(file.exists(args[2]))stop('Use a new output directory')
script<-normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)),winslash='/')
src<-file.path(dirname(script),'../../src')
cfg<-jsonlite::read_json(args[1],simplifyVector=TRUE);envs<-cfg$environments
requirements<-list(`3`='contrasts',`4`=c('cot_metrics','cot_components','maize_components','maize_flight','cotton_flight'),
  `6`=c('nested_summary','nested_repeats','loeo_metrics'),`7`='importance',`8`='correlations')
selected<-if(which=='all')names(requirements)else which
inputs<-vapply(unique(unlist(requirements[selected])),function(n){
  entry<-cfg$inputs[[n]];p<-entry$path
  if(!grepl('^([A-Za-z]:[/\\\\]|/)',p))p<-file.path(dirname(args[1]),p)
  p<-normalizePath(p,winslash='/',mustWork=TRUE)
  stopifnot(digest::digest(file=p,algo='sha256')==entry$sha256);p
},character(1))
if(.Platform$OS.type=='windows')grDevices::windowsFonts(Arial=grDevices::windowsFont('Arial'))
dir.create(args[2],recursive=TRUE)
modules<-c(`3`='figure3.R',`4`='figure4_support.R',`6`='figure6.R',`7`='figure7_support.R')
for(number in selected){
  cowplot::set_null_device(if(number%in%c('6','7'))cowplot::png_null_device else cowplot::pdf_null_device)
  out<-file.path(args[2],paste0('Figure',number));message('Rendering ',number)
  if(number=='8') {
    status<-system2(file.path(R.home('bin'),if(.Platform$OS.type=='windows')'Rscript.exe'else'Rscript'),c('--vanilla',
      shQuote(file.path(dirname(script),'render_correlations.R')),shQuote(inputs[['correlations']]),shQuote(out)))
    stopifnot(status==0L)
  } else sys.source(file.path(src,modules[[number]]),envir=environment())
}
jsonlite::write_json(list(status='COMPLETE',inputs=cfg$inputs,numerical_refitting=FALSE,
  numbering='Figures 3,4,6,7,8; supporting figures S1-S25',dpi=600),
  file.path(args[2],'render_manifest.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(args[2],'sessionInfo.txt'))
