loader <- function(theFile, envir=.mdEnv, pkgname="MinimumDistance"){
	theFile <- file.path(system.file(package=pkgname),
			     "extdata", theFile)
	if (!file.exists(theFile))
		stop("File ", theFile, " does not exist in ", pkgname)
	load(theFile, envir=envir)
}

isLoaded <- function(dataset, environ=.mdEnv)
	exists(dataset, envir=environ)

getVarInEnv <- function(dataset, environ=.mdEnv){
	if (!isLoaded(dataset, environ=environ))
		stop("Variable ", dataset, " not found in supplied environment")
	environ[[dataset]]
}
