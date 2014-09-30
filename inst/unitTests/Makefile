TOP=../..
PKG=${shell cd ${TOP};pwd}
SUITE=doRUnit.R

all: inst test

inst: # Install package
	cd ${TOP}/..;\
	R CMD INSTALL ${PKG}

test: # Run unit tests
	export RCMDCHECK=FALSE;\
	export RUNITFILEPATTERN="$(file)";\
	cd ${TOP}/tests;\
	R --vanilla --slave < ${SUITE}

release:
	R_LIBS_USER=~/Library/R/3.1-bioc-release/library R-3.1 CMD INSTALL --no-lock ${PKG}
	export RCMDCHECK=FALSE;\
	export RUNITFILEPATTERN="$(file)";\
	cd ${TOP}/tests;\
	R_LIBS_USER=~/Library/R/3.0-bioc-release/library R-3.1 --vanilla --slave < ${SUITE}

devel:
	R_LIBS_USER=~/Library/R/3.1-bioc-devel/library R-3.1 CMD INSTALL --no-lock ${PKG}
	export RCMDCHECK=FALSE;\
	export RUNITFILEPATTERN="$(file)";\
	cd ${TOP}/tests;\
	R_LIBS_USER=~/Library/R/3.1-bioc-devel/library R-3.1 --vanilla --slave < ${SUITE}

installRelease:
	R_LIBS_USER=~/Library/R/3.0-bioc-release/library R-3.0 CMD INSTALL --no-lock ${PKG}
