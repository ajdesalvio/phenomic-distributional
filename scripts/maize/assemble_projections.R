# Rscript scripts/maize/assemble_projections.R NEW.csv PART1.csv ... PART5.csv
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(TRUE); stopifnot(length(args)>=2L)
if(file.exists(args[1])) stop('Choose a new output file')
d <- rbindlist(lapply(args[-1],fread),use.names=TRUE,fill=TRUE)
# Environment occupies the final two dotted fields; genotype may contain dots.
d[,Env:=sub('^.*[.]([^.]+[.][^.]+)$','\\1',Pedigree.Env)]
d[,Pedigree:=substr(Pedigree.Env,1L,nchar(Pedigree.Env)-nchar(Env)-1L)]
d[,Env.VI:=paste(Env,Vegetation.Index,sep='.')]
stopifnot(!anyDuplicated(d[,.(Pedigree.Env,Heldout.Env,Vegetation.Index,Quantile.Level)]))
dir.create(dirname(args[1]),recursive=TRUE,showWarnings=FALSE)
fwrite(d,args[1])
