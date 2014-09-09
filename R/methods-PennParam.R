#' Constructor for class \code{PennParam}
#'
#' Parameters for the PennCNV Hidden Markov model
#'
#' @param states character vector of possible trio states
#' @param referenceState the reference (normal) trio copy number state (typically '222')
#' @param prLessLikelyCN as defined in Wang et al. 2007, this is the probability of the less likely allele-specific copy numbers for the trio
#' @param prNonMendelian the prior probability of a non-Mendelian copy number alteration
#' @param prInitialStateNotDiploid initial probability for non-diploid states
#' @param prTransitionToNewState probability of transitioning to a new state
#' @param tauNM probability of transitioning from a Mendelian given previous event was non-Mendelian (and vice versa).
#' @references Wang et al., Genome Res. 2007 Nov;17(11):1665-74.  PMID: 17921354
#' @export
PennParam <- function(states, referenceState="222", prLessLikelyCN=0.0009,
                      prNonMendelian=1.5e-6,
                      ##minimum_distance_threshold=0.9,
                      prInitialStateNotDiploid=4/5, ## uniform
                      prTransitionToNewState=0.5,
                      tauNM=0.01){
                      ##minimum_MAD=0.1){
  if(missing(states)) states <- trioStates(0:4)
  state_names <- trioStateNames(states)
  rownames(states) <- state_names
  table1 <- readTable1(a=prLessLikelyCN)
  table1v <- vectorizeTable1(table1, states)
  path <- system.file("extdata", package="MinimumDistance", mustWork=TRUE)
  load(file.path(path, "pennCNV_MendelianProb.rda"))
  initial_probs <- setNames(.initialStateProbs(5L, normal.index=3, prInitialStateNotDiploid),
                            paste0("CN:", 0:4))
  transition_probs <- transitionProbability(5L, epsilon=prTransitionToNewState) ## 5 states

  trNM <- setNames(rep(NA, 4), c("NM=0,0", "NM=0,1", "NM=1,0", "NM=1,1"))
  trNM["NM=0,0"] <- (1-tauNM)*(1-prNonMendelian)
  trNM["NM=0,1"] <- (tauNM)*(1-prNonMendelian)
  trNM["NM=1,0"] <- (tauNM)*(prNonMendelian)
  trNM["NM=1,1"] <- (1-tauNM)*(prNonMendelian)

  new("PennParam",
      table1=table1v,
      table3=pennCNV_MendelianProb,
      transitionNM=trNM,
      states=states,
      names=state_names,
      referenceState=referenceState,
      ##minimum_distance_threshold=minimum_distance_threshold,
      prNonMendelian=prNonMendelian,
      initialStateProb=initial_probs,
      transitionProb=transition_probs)
      ##minimum_MAD=minimum_MAD)
}

setGeneric("transitionNM", function(object) standardGeneric("transitionNM"))
setMethod("transitionNM", "PennParam", function(object) object@transitionNM)
setMethod("transitionNM", "MinDistParam", function(object) penncnv(object)@transitionNM)

setValidity("PennParam", function(object){
  msg <- TRUE
  x <- state(object)
  if(nrow(x) != 121){
    return("There should be 121 trio states")
  }
  msg
})

.initialStateProbs <- function(nstates, normal.index=3, epsilon=0.01){
  initial.state.probs <- rep(epsilon/(nstates-1), nstates)
  initial.state.probs[normal.index] <- 1-epsilon
  initial.state.probs
}

setMethod("table1", "PennParam", function(object) object@table1)
setMethod("table3", "PennParam", function(object) object@table3)
setMethod("state", "PennParam", function(object) object@states)
setMethod("stateNames", "PennParam", function(object) object@names)
setMethod("referenceState", "PennParam", function(object) object@referenceState)
setMethod("prNonMendelian", "PennParam", function(object) object@prNonMendelian)
setMethod("minimum_distance_threshold", "PennParam", function(object) object@minimum_distance_threshold)
setMethod("initialStateProb", "PennParam", function(object) object@initialStateProb)

setReplaceMethod("initialStateProb", "PennParam", function(object, value) {
  object@initialStateProb <- value
  object
})

setMethod("transitionProb", "PennParam", function(object) object@transitionProb)
setMethod("minimum_MAD", "PennParam", function(object) object@minimum_MAD)

setMethod("show", "PennParam", function(object){
  cat("Object of class `PennParam'\n")
  cat("  ", nrow(state(object)), " trio states\n")
  names <- paste(head(stateNames(object)), collapse=",")
  cat("   state names:", names, "...\n")
  cat("   reference state:", referenceState(object), "\n")
  cat("   PennCNV Table 1: ", length(table1(object)), " vector\n")
  cat("   PennCNV Table 3: 5 x 5 x 5 x 5 x 5 x 5 array \n")
  cat("   probability non-Mendelian:", prNonMendelian(object), "\n")
  cat("   minimum distance threshold:", minimum_distance_threshold(object), "\n")
  pis <- paste(round(initialStateProb(object), 2), collapse=",")
  cat("   initial state probabilities: ", pis, "\n")
  cat("   transition prob: 5 x 5 matrix \n")
  cat("   See table1(), table3(), state(), stateNames(), referenceState()\n")
})
