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

  # extract annual flows
  flows <- lapply(cumulative_flows, function(x) apply(x, 2, diff))
  names(flows) <- gsub("^cum_", "", names(cumulative_flows))

  # remove time 0 from cumulative flows
  cumulative_flows <- lapply(cumulative_flows, "[", -1, )

  c(cumulative_flows, flows)
}

##' Convert a year into the number of years after 2007
##'
##' @title Convert a year into the number of years after 2007
##'
##' @param year an integer year
##'
##' @return An integer, being the number of years after 2007
##' @export
##' @examples
##' gonovax_year(2019)
##' gonovax_year(c(2018, 2019))
gonovax_year <- function(year) {
  years_after_2007 <- as.numeric(year - 2007)
  if (any(years_after_2007 < 0)) {
    stop("Negative dates, gonovax_year likely applied twice")
  }
  years_after_2007
}

gonovax_year_as_year <- function(gonovax_year) {
  assert_positive_integer(gonovax_year)
  2007 + gonovax_year
}


##' pdf of a betabinomial parametrised in terms of mean and over-dispersion
##'
##' @title pdf of a betabinomial parametrised in terms of probability and
##' over-dospersion
##'
##' @param x data
##'
##' @param size integer of sample size
##'
##' @param prob probability of observing a single success
##'
##' @param rho overdispersion parameter
##'
##' @return probability of observing x
##' @export

dbetabinom <- function(x, size, prob, rho, log = FALSE) {

  ## here we have prob = a / (a + b)
  ## so and rho = (a + b + n) / (a + b + 1)
  a <- prob * (1 / rho - 1)
  b <- (1 - prob) * (1 / rho - 1)

  out <- lchoose(size, x) + lbeta(x + a, size - x + b) - lbeta(a, b)

  if (!log) {
    out <- exp(out)
  }

  out
}
