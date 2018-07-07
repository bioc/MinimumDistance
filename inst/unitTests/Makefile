## Note that the tests are run on the scripts in the installed directly.
## Need to run 'make' to test any changes to the files in inst/unitTests
TOP=../..
PKG=${shell cd ${TOP};pwd}
SUITE=doRUnit.R
R=${HOME}/bin/R

all: inst test

inst: # Install package
	cd ${TOP}/..;\
	R CMD INSTALL ${PKG}

test: # Run unit tests
	export RCMDCHECK=FALSE;\
	export RUNITFILEPATTERN="$(file)";\
	cd ${TOP}/tests;\
	R_LIBS_USER=~/Library/R/3.8-bioc-devel R --vanilla --slave < ${SUITE}

release: # Run unit tests
	cd ${TOP}/..;\
	R_LIBS_USER=~/Library/R/3.7-bioc-release R CMD INSTALL ${PKG}
	export RCMDCHECK=FALSE;\
	export RUNITFILEPATTERN="$(file)";\
	cd ${TOP}/tests;\
	R_LIBS_USER=~/Library/R/3.7-bioc-release R --vanilla --slave < ${SUITE}

devel:
	R_LIBS_USER=~/Library/R/3.8-bioc-devel $R CMD INSTALL --no-lock ${PKG}
	export RCMDCHECK=FALSE;\
	export RUNITFILEPATTERN="$(file)";\
	cd ${TOP}/tests;\
	R_LIBS_USER=~/Library/R/3.8-bioc-devel $R --vanilla --slave < ${SUITE}

install:
	R_LIBS_USER=~/Library/R/3.8-bioc-devel R CMD INSTALL --no-lock ${PKG}
