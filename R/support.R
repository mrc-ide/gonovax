extract_value <- function(x, what, t) {
  sapply(x, function(x) sapply(x, function(x) {
    i <- which(x$t == t)
    sum(x[[what]][i, , ])
  }))
}

sum_group_vax <- function(x, what) {
  sapply(x, function(x) apply(x[[what]], 1, sum))
}