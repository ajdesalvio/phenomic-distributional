# Rscript scripts/maize/run_nested_cv.R config/nested_cv.json NEW_OUTPUT [ENV]
# Each environment runs separately, preserving the manuscript's fold seeds.
nested_config <- function(inputs, env, output, sources) {
  list(version="V15", mode="full", seed=123,
       envs=env, target_vi="ExGR", outer_folds=5, outer_repeats=5, inner_folds=5,
       alphas=seq(0,1,by=0.1), lambdas=exp(seq(log(10),log(1e-4),length.out=100)),
       scenarios=c("Median","all"), adaptive=FALSE,
       predictors=inputs[['blue']], raw_support=inputs[['support']],
       preprocessing_policy="Raw-support NZV 19/10; numerical-constant roundoff bound V1",
       yield=inputs[['yield']], output=output, max_outer_splits=Inf,
       script_path=normalizePath(sources,winslash="/",mustWork=TRUE))
}

nested_main <- function(args=commandArgs(TRUE)) {
  if (!length(args) %in% 2:3) stop("Usage: CONFIG.json NEW_OUTPUT [ENV]")
  if (!requireNamespace("jsonlite",quietly=TRUE)) stop("Install jsonlite first")
  own <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
  engine <- file.path(dirname(own),"..","..","src","nested_cv.R")
  source(engine)
  if (Sys.setlocale("LC_ALL","C") != "C") stop("C locale is required for the recorded fold ordering")
  config <- jsonlite::fromJSON(args[1],simplifyVector=FALSE)
  roles <- c("blue","support","yield")
  if (!setequal(names(config$inputs),roles)) stop("Expected blue, support, and yield inputs")
  paths <- vapply(config$inputs,function(x)x$path,character(1))
  relative <- !grepl('^([A-Za-z]:[/\\]|/|\\\\)',paths)
  paths[relative] <- file.path(dirname(normalizePath(args[1])),paths[relative])
  paths <- vapply(paths,normalizePath,character(1),winslash="/",mustWork=TRUE)
  expected <- vapply(config$inputs,function(x)x$md5,character(1))
  if (!identical(tolower(unname(tools::md5sum(paths))),tolower(unname(expected)))) stop("Frozen input hash mismatch")
  envs <- unlist(config$environments,use.names=FALSE)
  allowed <- c("IAH4.2021","MIH1.2020","MNH1.2020","MNH1.2021","TXH1.2020","TXH2.2021",
               "WIH1.2021","WIH2.2020","WIH2.2021","WIH3.2020","WIH3.2021")
  if (!length(envs)||anyNA(envs)||anyDuplicated(envs)||any(!envs %in% allowed)) stop("Invalid environment scope")
  if (length(args)==3L) {
    if (!args[3] %in% envs) stop("Environment absent from configuration")
    envs <- args[3]
  }
  for (env in envs) {
    output <- file.path(args[2],env)
    if (dir.exists(output) && length(list.files(output,all.files=TRUE,no..=TRUE)) &&
        !file.exists(file.path(output,"Run_Manifest_V15.rds"))) stop("Existing output is not a resumable V15 run: ",output)
    cfg <- nested_config(paths,env,output,c(engine,own))
    v15_run(cfg)
  }
}
if (sys.nframe()==0L) nested_main()
