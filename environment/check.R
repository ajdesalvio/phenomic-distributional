# Run from the package root: Rscript environment/check.R
expected <- read.csv("environment/R_versions.csv", stringsAsFactors=FALSE)
actual <- vapply(expected$Package, function(p) {
  if (p == "R") return(as.character(getRversion()))
  if (requireNamespace(p, quietly=TRUE)) as.character(packageVersion(p)) else NA_character_
}, character(1))
report <- data.frame(Package=expected$Package, Tested=expected$Version,
                     Installed=unname(actual), Matches=unname(actual)==expected$Version)
print(report, row.names=FALSE)
if (anyNA(actual)) stop("Install missing packages before running the workflows")
if (actual[["clubSandwich"]] != "0.7.0") stop("Functional analysis requires clubSandwich 0.7.0")
if (any(!report$Matches)) warning("Versions differ from the tested environment; this check does not certify numerical parity")
