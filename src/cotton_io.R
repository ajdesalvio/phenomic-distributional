cotton_input <- function(entry, config) {
  p <- entry$path
  if(!grepl('^([A-Za-z]:[/\\\\]|/)',p)) p <- file.path(dirname(config),p)
  p <- normalizePath(p,winslash='/',mustWork=TRUE)
  actual <- if(!is.null(entry$md5)) unname(tools::md5sum(p)) else digest::digest(file=p,algo='sha256')
  expected <- if(!is.null(entry$md5)) entry$md5 else entry$sha256
  if(!identical(tolower(actual),tolower(expected))) stop('Input hash mismatch: ',p)
  p
}
