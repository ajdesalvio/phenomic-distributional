# Rscript scripts/maize/export_loeo_reference.R CONFIG.json NEW_OUTPUT
suppressPackageStartupMessages({library(data.table); library(glmnet)})
args <- commandArgs(TRUE)
stopifnot(length(args)==2L)
cfg <- jsonlite::read_json(args[1],simplifyVector=TRUE)
path <- cfg$archive$path
if (!grepl('^([A-Za-z]:[/\\\\]|/)',path)) path <- file.path(dirname(args[1]),path)
stopifnot(digest::digest(file=path,algo="sha256")==cfg$archive$sha256)
if (file.exists(args[2])) stop("Choose a new output directory")
saved <- readRDS(path)
keys <- as.vector(outer(cfg$environments,c("Median","all"),
                       function(e,s)paste(e,"Yield.t.ha.BLUE",s,sep="_")))
stopifnot(all(keys %in% names(saved$Final_Models)),all(keys %in% names(saved$LOEO_Preds)))
metrics <- predictions <- coefficients <- list()
for (key in keys) {
  d <- as.data.table(saved$LOEO_Preds[[key]])
  model <- saved$Final_Models[[key]]
  stopifnot(length(model$lambdaOpt)==1L,is.finite(model$lambdaOpt))
  cf <- as.matrix(coef(model,s=model$lambdaOpt))[,1]
  cf <- cf[names(cf)!="(Intercept)"]
  z <- d[is.finite(Pred)&is.finite(Obs)]
  stopifnot(nrow(z)>=3L,!anyDuplicated(d$Pedigree.Env),uniqueN(d$Heldout)==1L,uniqueN(d$Scenario)==1L)
  metrics[[key]] <- data.table(Heldout=d$Heldout[1],Scenario=d$Scenario[1],Scaling="WITHIN",
    n_test=nrow(z),Pearson=cor(z$Pred,z$Obs),Spearman=cor(z$Pred,z$Obs,method="spearman"),
    RMSE=sqrt(mean((z$Pred-z$Obs)^2)),MAE=mean(abs(z$Pred-z$Obs)),
    SD_Pred=sd(z$Pred),SD_Obs=sd(z$Obs),NZ_Coefs=sum(abs(cf)>0),
    SourceStatus="historical_not_refit")
  predictions[[key]] <- d[,.(Heldout,Scenario,Pedigree.Env,Obs,Pred,Residual)]
  coefficients[[key]] <- data.table(Heldout=d$Heldout[1],Scenario=d$Scenario[1],
    Feature=names(cf),Coefficient=unname(cf),Selected=abs(cf)>0)
}
metrics <- rbindlist(metrics)
stopifnot(nrow(metrics)==22L,!anyDuplicated(metrics[,.(Heldout,Scenario)]))
summary <- metrics[,.(Heldouts=.N,Mean_Pearson=mean(Pearson),Mean_Spearman=mean(Spearman),
  Mean_RMSE=mean(RMSE),Mean_MAE=mean(MAE)),by=.(Scaling,Scenario)]
dir.create(args[2],recursive=TRUE)
for (name in c("metrics","predictions","coefficients","summary")) {
  x <- get(name); if(is.list(x)&&!is.data.frame(x)) x <- rbindlist(x)
  fwrite(x,file.path(args[2],paste0("LOEO_",name,".csv")))
}
writeLines(capture.output(sessionInfo()),file.path(args[2],"sessionInfo.txt"))
jsonlite::write_json(list(archive=normalizePath(path,winslash='/'),sha256=cfg$archive$sha256,
  models=22L,status='historical_not_refit',metrics='Recalculated from archived predictions'),
  file.path(args[2],'provenance.json'),pretty=TRUE,auto_unbox=TRUE)
