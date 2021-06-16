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

##' @name extract_flows
##' @title extract flows used for run_grid
##' @param y a transformed model run output
##' @return cumulative and incident flows
##' @export
extract_flows <- function(y) {

  # extract cumulative flows
  flow_names <- c("cum_diag_a", "cum_diag_s", "cum_treated", "cum_screened",
                  "cum_vaccinated")
  cumulative_flows <- lapply(flow_names, function(x) t(aggregate(y, x)))
  names(cumulative_flows) <- flow_names

  # extract vaccinations and revaccinations separately
  cum_newly_vaccinated <- t(aggregate(y, "cum_vaccinated", stratum = 1))
  cumulative_flows$cum_revaccinated <-
    cumulative_flows$cum_vaccinated - cum_newly_vaccinated

  # extract annual flows
  flows <- lapply(cumulative_flows, function(x) apply(x, 2, diff))
  names(flows) <- gsub("^cum_", "", names(cumulative_flows))

  # remove time 0 from cumulative flows and remove those not needed
  cumulative_flows <-
    lapply(cumulative_flows[c("cum_treated", "cum_vaccinated")], "[", -1, )

  c(cumulative_flows, flows)
}

##' Convert a year into the number of years after 2009
##'
##' @title Convert a year into the number of years after 2009
##'
##' @param year an integer year
##'
##' @return An integer, being the number of years after 2009
##' @export
##' @examples
##' gonovax_year(2019)
##' gonovax_year(c(2018, 2019))
gonovax_year <- function(year) {
  years_after_2009 <- as.numeric(year - 2009)
  if (any(years_after_2009 < 0)) {
    stop("Negative dates, gonovax_year likely applied twice")
  }
  years_after_2009
}

##' @title Convert a gonovax year into calendar years
##' @param gonovax_year an integer
##' @return An integer, being the calendar year
##' @export
##' @examples
##' gonovax_year_as_year(3)
gonovax_year_as_year <- function(gonovax_year) {
  assert_positive_integer(gonovax_year)
  2009 + gonovax_year
}


##' @title pdf of a betabinomial parametrised in terms of probability and
##' over-dispersion
##' @param x data
##' @param size integer of sample size
##' @param prob probability of observing a single success
##' @param rho overdispersion parameter, has support [0, 1]
##'  with low values being less overdispersion
##' @param log logical indicating whether to return log value
##'
##' @return probability of observing x
##' @export

dbetabinom <- function(x, size, prob, rho, log = FALSE) {

  ##comparison with rmutil::dbetabinom
  ## s = (1 / rho - 1), so that
  ## rho = 1 / (s + 1).

  a <- prob * (1 / rho - 1)
  b <- (1 - prob) * (1 / rho - 1)

  out <- lchoose(size, x) + lbeta(x + a, size - x + b) - lbeta(a, b)

  if (!log) {
    out <- exp(out)
  }

  out
}

name_outputs <- function(res, strata_names) {

  group_names <- c("L", "H")
  state_names <- c("U", "I", "A", "S", "T", "N",
                   "cum_incid", "cum_diag_a", "cum_diag_s",
                   "cum_treated", "cum_screened", "cum_vaccinated")

  for (nm in state_names) {
    dimnames(res[[nm]]) <- list(NULL, group_names, strata_names)
  }
  dimnames(res$lambda) <- list(NULL, group_names)
  dimnames(res$eta) <- list(NULL, group_names)

  res
}
