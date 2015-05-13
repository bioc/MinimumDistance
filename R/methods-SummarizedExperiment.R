## these methods are used internally...

#' @aliases father,RangedSummarizedExperiment-method
#' @rdname MinDistExperiment-class
setMethod("father", "RangedSummarizedExperiment", function(object) assays(object)[["father"]])

#' @aliases mother,RangedSummarizedExperiment-method
#' @rdname MinDistExperiment-class
setMethod("mother", "RangedSummarizedExperiment", function(object) assays(object)[["mother"]])

#' @aliases offspring,RangedSummarizedExperiment-method
#' @rdname MinDistExperiment-class
setMethod("offspring", "RangedSummarizedExperiment", function(object) assays(object)[["offspring"]])
