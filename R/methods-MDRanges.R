.data_frame_posteriorSummaries <- function(g){
  L <- length(g)
  x <- matrix(NA, L, 4)
  colnames(x) <-  c("log_RR",
                    "log_odds",
                    "prob_MAP",
                    "prob_222")
  rownames(x) <- names(g)
  as.data.frame(x)
}

.mcols_posteriors <- function(g){
  L <- length(g)
  df <- DataFrame(log_RR=numeric(L),
                  log_odds=numeric(L),
                  prob_MAP=numeric(L),
                  prob_222=numeric(L),
                  prob_221=numeric(L))
}

.mdranges_mcols <- function(g){
  L <- length(g)
  df <- DataFrame(sample=character(),
                  numberProbes=integer(),
                  seg.mean=numeric(),
                  log_RR=numeric(L),
                  log_odds=numeric(L),
                  prob_MAP=numeric(L),
                  prob_222=numeric(L),
                  prob_221=numeric(L))
}

MDRanges <- function(..., posteriors){
  g <- GRanges(...)
  if(missing(posteriors)) posteriors <- .mdranges_mcols(g)
  g2 <- as(g, "MDRanges")
  mcols(g2) <- posteriors
  g2
}



setMethod("posteriorLogOdds", "MDRanges", function(object) object$log_odds)

posteriorOdds <- function(object) exp(posteriorLogOdds(object))
posteriorMAP <- function(object) object$prob_MAP
posteriorRR <- function(object) exp(posteriorLogRR(object))
posterior221 <- function(object) object$prob_221

setMethod("posteriorLogRR", "MDRanges", function(object) object$log_RR)

## #' @export
isDenovo <- function(states) (states %in% c(duplicationStates(), deletionStates())) & !is.na(states)

setMethod("is221", "MDRanges", function(object)  object$calls=="221" & !is.na(object$calls))

setMethod("is221", "GRangesList", function(object) lapply(object, is221))
