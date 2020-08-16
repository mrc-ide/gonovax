extract_value <- function(x, what, t) {
  sapply(x, function(x) sapply(x, function(x) {
    i <- which(x$t == t)
    sum(x[[what]][i, , , ])
  }))
}
