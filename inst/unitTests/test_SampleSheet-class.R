test_SampleSheet_construction <- function(){
	##checkException(SampleSheet(), silent=TRUE)  ## use if expect an error
	checkTrue(validObject(SampleSheet()))
	checkTrue(validObject(new("SampleSheet")))
	checkTrue(validObject(SampleSheet(row.names=letters)))
	data(swiss)
	checkTrue(validObject(SampleSheet(swiss)))
}
