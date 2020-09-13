##' @name aggregate
##' @title aggregate model runs by vaccination group x vaccine strata
##' @param x a transformed model run output
##' @param what a character string of the cumulative trajectory to transform
##' @param as_incid logical specifying whether to convert to incidence
##' @param stratum an integer denoting the stratum to aggregate over, if NULL
##' aggregation is over all strata
##' @param f function to apply across parameter sets (eg. mean / median)
##' @param ... named arguments to pass to f
##' @return a transformed time series / array
##' @export
aggregate <- function(x, what, as_incid = FALSE, stratum = NULL,
                      f = identity, ...) {
  strata <- stratum %||% seq_len(dim(x[[1]][[what]])[3])
  y <- sapply(x, function(x) apply(x[[what]][, , strata], 1, sum))
  if (as_incid) {
    dt <- diff(x[[1]]$t)
    y <- apply(y, 2, diff) / dt
  }
  apply(y, 1, f, ...)
}

extract_flows <- function(y) {

  # extract cumulative flows
  flow_names <- c("cum_incid", "cum_diag_a", "cum_diag_s",
                  "cum_treated", "cum_screened", "cum_vaccinated")
  cumulative_flows <- lapply(flow_names, function(x) t(aggregate(y, x)))
  names(cumulative_flows) <- flow_names

  # extract vaccinations and revaccinations separately
  cum_newly_vaccinated <- t(aggregate(y, "cum_vaccinated", stratum = 1))
  cumulative_flows$cum_revaccinated <-
    cumulative_flows$cum_vaccinated - cum_newly_vaccinated
  cumulative_flows$cum_vaccinated <- cum_newly_vaccinated

  # extract annual flows
  flows <- lapply(cumulative_flows, function(x) apply(x, 2, diff))
  names(flows) <- gsub("^cum_", "", names(cumulative_flows))

  # remove time 0 from cumulative flows
  cumulative_flows <- lapply(cumulative_flows, "[", -1, )

  c(cumulative_flows, flows)
}
