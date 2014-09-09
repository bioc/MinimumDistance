make_test_Pedigree <- function(){
  pedTest <- new("Pedigree", trios=data.frame(F=c("NA06993", "NA11881"),
                               M=c("NA06985", "NA11882"),
                               O=c("NA06991", "NA10859"),
                               stringsAsFactors=FALSE),
                 trioIndex=data.frame(individualId=c("NA06993", "NA11881",
                                        "NA06985", "NA11882", "NA06991", "NA10859"),
                   memberId=rep(c("F", "M", "O"), each=2),
                   index.in.pedigree=rep(1:2, 3),
                   stringsAsFactors=FALSE))
}

test_Pedigree_construction <- function(){
  library(oligoClasses)
  checkTrue(validObject(Pedigree()))
  path <- system.file("extdata", package="MinimumDistance")
  load(file.path(path, "pedigreeInfo.rda"))
  ped <- Pedigree(fatherIds=pedigreeInfo$F,
                  motherIds=pedigreeInfo$M,
                  offspringIds=pedigreeInfo$O)
  checkTrue(validObject(ped))

  ped2 <- ped[2, ]
  checkTrue(validObject(ped2))
  checkTrue(validObject(Pedigree(pedigreeInfo)))

  ##validObject(ped[1, ])
  ped2 <- make_test_Pedigree()
  checkIdentical(ped, ped2)

  ## can not have duplicate offspring identifiers
  checkException(Pedigree(data.frame(F=c("F0.txt", "F0.txt"),
                                     M=c("M0.txt", "M0.txt"),
                                     O=c("O0.txt", "O0.txt"))), silent=TRUE)

  checkTrue(validObject(Pedigree(data.frame(F=c("F0.txt", "F0.txt"),
                                            M=c("M0.txt", "M0.txt"),
                                            O=c("O0.txt", "O1.txt")))))

  trio <- trios(ped)
  trio$F[[1]] <- "badname"
  ped@trios <- trio
  checkException(validObject(ped), silent=TRUE)

  ##checkTrue(identical(sampleNames(ped), offspringNames(ped)))
}

test_subsetPedigree <- function(){
  object <- Pedigree(fatherIds=rep(letters[1:3], each=2),
                     motherIds=rep(letters[4:6], each=2),
                     offspringIds=letters[11:16])
  checkTrue(validObject(object[6, ]))
  ## offspring is both parent and an offspring in the same trio
  ## (should throw an error)
  checkException(Pedigree(fatherIds=rep(letters[1:2]),
                          motherIds=rep(letters[3:4]),
                          offspringIds=letters[c(1,5)]))
  ## offspring is a parent in one trio and an offspring in a
  ## different trio (should be valid)
  ped <- Pedigree(fatherIds=rep(letters[1:2]),
                  motherIds=rep(letters[3:4]),
                  offspringIds=letters[c(5,1)])
  checkTrue(validObject(ped))

  validObject(Pedigree(data.frame(F=letters[1],
                                  M=letters[2],
                                  O=letters[3])))
}
