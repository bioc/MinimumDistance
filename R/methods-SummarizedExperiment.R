setMethod("father", "SummarizedExperiment", function(object) assays(object)[["father"]])
setMethod("mother", "SummarizedExperiment", function(object) assays(object)[["mother"]])
setMethod("offspring", "SummarizedExperiment", function(object) assays(object)[["offspring"]])
