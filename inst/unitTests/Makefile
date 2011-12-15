TOP=../..
PKG=${shell cd ${TOP};pwd}
SUITE=doRUnit.R
#R=${R_HOME}/bin/R
R=~/bin/R-devel

all: inst test

inst: # Install package
	cd ${TOP}/..;\
	##${R} CMD INSTALL ${PKG}
	/home/bst/student/rscharpf/bin/idevel ${PKG}

test: # Run unit tests
	export RCMDCHECK=FALSE;\
	export RUNITFILEPATTERN="$(file)";\
	cd ${TOP}/tests;\
	${R} --vanilla --slave < ${SUITE}