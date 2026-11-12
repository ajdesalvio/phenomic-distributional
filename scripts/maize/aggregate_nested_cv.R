# Rscript scripts/maize/aggregate_nested_cv.R RUN_ROOT NEW_SUMMARY_DIRECTORY
args <- commandArgs(TRUE)
if (length(args)!=2L) stop("Supply RUN_ROOT NEW_SUMMARY_DIRECTORY")
own <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
root <- file.path(dirname(own),"..","..")
source(file.path(root,"src","aggregate_nested_cv.R"))
v15a_run(args[1],args[2],c(file.path(root,"src","nested_cv.R"),
                          file.path(dirname(own),"run_nested_cv.R")))
