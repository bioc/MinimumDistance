make_test_Pedigree <- function(){
	new("Pedigree")
}

test_Pedigree_construction <- function(){
	##checkException(Pedigree(), silent=TRUE)
	checkTrue(validObject(Pedigree()))
	checkTrue(validObject(new("Pedigree")))
}
