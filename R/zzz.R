THISPKG <- "MinimumDistance"
.mdEnv <- new.env(parent=emptyenv())

.onAttach <- function(libname, pkgname) {
  version <- packageDescription("MinimumDistance", fields="Version")
  packageStartupMessage(paste("Welcome to MinimumDistance version ", version))
}
