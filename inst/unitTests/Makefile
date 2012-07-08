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
	R_LIBS_USER=~/Library/R/2.15-bioc-release/library R-2.15 CMD INSTALL --no-lock ${PKG}
	export RCMDCHECK=FALSE;\
	export RUNITFILEPATTERN="$(file)";\
	cd ${TOP}/tests;\
	R_LIBS_USER=~/Library/R/2.15-bioc-release/library R-2.15 --vanilla --slave < ${SUITE}

devel:
	R_LIBS_USER=~/Library/R/2.15-bioc-devel/library R-2.15 CMD INSTALL --no-lock ${PKG}
	export RCMDCHECK=FALSE;\
	export RUNITFILEPATTERN="$(file)";\
	cd ${TOP}/tests;\
	R_LIBS_USER=~/Library/R/2.15-bioc-devel/library R-2.15 --vanilla --slave < ${SUITE}

