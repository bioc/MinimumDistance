## these methods are used internally...

#' @aliases father,SummarizedExperiment-method
#' @rdname MinDistExperiment-class
setMethod("father", "SummarizedExperiment", function(object) assays(object)[["father"]])

#' @aliases mother,SummarizedExperiment-method
#' @rdname MinDistExperiment-class
setMethod("mother", "SummarizedExperiment", function(object) assays(object)[["mother"]])

#' @aliases offspring,SummarizedExperiment-method
#' @rdname MinDistExperiment-class
setMethod("offspring", "SummarizedExperiment", function(object) assays(object)[["offspring"]])
