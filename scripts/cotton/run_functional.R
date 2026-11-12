# Rscript scripts/cotton/run_functional.R CONFIG.json VI OUTPUT inspect|smoke|full|validate [resume]
args <- commandArgs(TRUE)
if(!length(args)%in%4:5||!args[4]%in%c('inspect','smoke','full','validate'))
  stop('Usage: CONFIG.json VI OUTPUT inspect|smoke|full|validate [resume]')
if(length(args)==5L&&args[5]!='resume') stop('The last argument must be resume')
script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)),winslash='/')
src <- file.path(dirname(script),'../../src')
source(file.path(src,'cotton_io.R'));source(file.path(src,'cotton_functional.R'))
config <- jsonlite::read_json(args[1],simplifyVector=TRUE)
stopifnot(args[2]%in%c('PSRI','CVI','ExG2'))
cfg <- config$functional
cfg$vi <- args[2];cfg$mode <- args[4];cfg$out <- args[3];cfg$resume <- length(args)==5L
cfg$script_path <- script;cfg$engine_path <- normalizePath(file.path(src,'cotton_functional.R'),winslash='/')
cfg$io_path <- normalizePath(file.path(src,'cotton_io.R'),winslash='/')
cfg$input <- cotton_input(config$inputs[[cfg$vi]],args[1])
cfg$fieldmap <- cotton_input(config$fieldmap,args[1])
if(cfg$mode=='smoke') cfg$draws <- 199L
result <- v7_run(cfg)
if(cfg$mode=='validate') print(result)
