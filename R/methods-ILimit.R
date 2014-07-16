setMethod("range", "ILimit", function(x, ...) c(start(x), end(x)))

setMethod("seq_along2", "ILimit", function(along.with){
  seq(start(along.with), end(along.with), 1)
})

ILimit <- function(...) as(IRanges(...), "ILimit")
