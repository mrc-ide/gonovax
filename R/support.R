##' @name aggregate
##' @title aggregate model runs by vaccination group x vaccine strata
##' @param x a transformed model run output
##' @param what a character string of the cumulative trajectory to transform
##' @param as_incid logical specifying whether to convert to incidence
##' @param f function to apply across parameter sets (eg. mean / median)
##' @param ... named arguments to pass to f
##' @return a transformed time series / array
##' @export
aggregate <- function(x, what, as_incid = FALSE, f = identity, ...) {
  y <- sapply(x, function(x) apply(x[[what]], 1, sum))
  if (as_incid) {
    dt <- diff(x[[1]]$t)
    y <- apply(y, 2, diff) / dt
  }
  apply(y, 1, f, ...)
}

extract_flows <- function(y) {
  # extract cumulative flows
  flow_names <- grep(pattern = "^cum", names(y[[1]]), value = TRUE)
  cumulative_flows <- lapply(flow_names, function(x) t(aggregate(y, x))[-1, ])
  names(cumulative_flows) <- flow_names

  # extract flows
  flows <- lapply(flow_names, function(x) t(aggregate(y, x, TRUE)))
  names(flows) <- gsub("^cum_", "", flow_names)
  c(cumulative_flows, flows)
}
