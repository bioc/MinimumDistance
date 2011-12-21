loader <- function(theFile, envir=.mdEnv, pkgname="MinimumDistance"){
	theFile <- file.path(system.file(package=pkgname),
			     "extdata", theFile)
	if (!file.exists(theFile))
		stop("File ", theFile, " does not exist in ", pkgname)
	load(theFile, envir=envir)
}
