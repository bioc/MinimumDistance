THISPKG <- "MinimumDistance"
.mdEnv <- new.env(parent=emptyenv())

#' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  version <- packageDescription("MinimumDistance", fields="Version")
  packageStartupMessage(paste("Welcome to MinimumDistance version ", version))
}
