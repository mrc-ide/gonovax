##' @name aggregate
##' @title aggregate model runs by vaccination group x vaccine strata
##' @param x a transformed model run output
##' @param what a character string of the cumulative trajectory to transform
##' @param as_incid logical specifying whether to convert to incidence
##' @param stratum an integer denoting the stratum to aggregate over, if NULL
##' aggregation is over all strata
##' @param f function to apply across parameter sets (eg. mean / median)
##' @param stochastic indicates if function is being used on output of
##' stochastic or deterministic version of the trial model, only matters if
##' as_incid = TRUE
##' @param ... named arguments to pass to f
##' @return a transformed time series / array
##' @export
aggregate <- function(x, what, as_incid = FALSE, stratum = NULL,
                      f = identity, stochastic = FALSE, ...) {
  strata <- stratum %||% seq_len(dim(x[[1]][[what]])[3])
  y <- sapply(x, function(x) apply(x[[what]][, , strata], 1, sum))

  if (as_incid) {

    if (stochastic) {
      dt <- diff(x[[1]]$time)
    }else {
      dt <- diff(x[[1]]$t)
    }

    y <- apply(y, 2, diff) / dt
  }
  apply(y, 1, f, ...)
}

##' @name extract_flows_xpvwrh
##' @title extract flows used for run_grid when the branching  xpvwrh model
##' has been run
##' @param y a transformed model run output
##' @return cumulative and incident flows
##' @export
extract_flows_xpvwrh <- function(y) {
  # extract cumulative flows
  flow_names <- c("cum_diag_a", "cum_diag_s", "cum_treated", "cum_screened",
                  "cum_vaccinated", "cum_vbe")
  cumulative_flows <- lapply(flow_names, function(x) t(aggregate(y, x)))
  names(cumulative_flows) <- flow_names

  ## extract vaccinations and revaccinations separately
  # primary vaccinated = everyone vaccinated from X(1) regardless of # doses
  cumulative_flows$cum_primary_total <-
    t(aggregate(y, "cum_vaccinated", stratum = 1))

  # partial to full = everyone vaccinated from P(2)
  # in earlier models this will be V, we expect this to be 0
  cumulative_flows$cum_part_to_full <-
    t(aggregate(y, "cum_vaccinated", stratum = 2))

  # revaccinated = everyone vaccinated from W(4), note does not include
  # vaccination of individuals in X who have waned from P
  cumulative_flows$cum_revaccinated <-
    t(aggregate(y, "cum_vaccinated", stratum = 4))

  # extract annual flows
  flows <- lapply(cumulative_flows, function(x) apply(x, 2, diff))
  names(flows) <- gsub("^cum_", "", names(cumulative_flows))

  # remove time 0 from cumulative flows and remove those not needed
  cumulative_flows <-
    lapply(cumulative_flows[c("cum_treated", "cum_vaccinated")], "[", -1, )

  c(cumulative_flows, flows)
}


##' @name extract_flows
##' @title extract flows used for run_grid
##' @param y a transformed model run output
##' @return cumulative and incident flows
##' @export
extract_flows <- function(y) {

  # extract cumulative flows
  flow_names <- c("cum_diag_a", "cum_diag_s", "cum_treated", "cum_screened",
                  "cum_vaccinated", "cum_vbe")
  cumulative_flows <- lapply(flow_names, function(x) t(aggregate(y, x)))
  names(cumulative_flows) <- flow_names

  # extract vaccinations and revaccinations separately
  cum_newly_vaccinated <- t(aggregate(y, "cum_vaccinated", stratum = 1))
  cumulative_flows$cum_revaccinated <-
    cumulative_flows$cum_vaccinated - cum_newly_vaccinated

  # extract those offered vaccination for first time (not including VbE)
  cumulative_flows$cum_offered_primary <-
    t(aggregate(y, "cum_offered", stratum = 1)) -
    t(aggregate(y, "cum_offered_vbe", stratum = 1)) # i.e. non-hesitant

  # extract annual flows
  flows <- lapply(cumulative_flows, function(x) apply(x, 2, diff))
  names(flows) <- gsub("^cum_", "", names(cumulative_flows))

  # remove time 0 from cumulative flows and remove those not needed
  cumulative_flows <-
    lapply(cumulative_flows[c("cum_treated", "cum_vaccinated")], "[", -1, )

  c(cumulative_flows, flows)
}

##' @name extract_flows_trial
##' @title extract flows for the XVW trial model
##' @param y a transformed model run output
##' @return cumulative and incident flows
##' @export

extract_flows_trial <- function(y) {
  # extract cumulative flows (standard)
  flow_names <- c("cum_diag_a", "cum_diag_s", "cum_incid",
                  "cum_treated", "cum_screened",
                  "N")
  cumulative_flows <- lapply(flow_names, function(x) t(aggregate(y, x)))
  names(cumulative_flows) <- flow_names

  strata <- dimnames(y[[1]]$N)[[3]] # extract strata names
  n_diag_rec <- sum(grepl("W", strata)) # count diag hist categories
  n_erlang <- sum((grepl("V", strata))) / n_diag_rec # count erlang

  idx <- stratum_index_xvw_trial(n_erlang, n_diag_rec)
  idx$never_diag <- seq(idx$V[1], by = n_diag_rec, length.out = n_erlang + 1)

  ## extract strata separately
  #  asymptomatic diagnoses across all X
  cumulative_flows$cum_diag_a_X_all_diaghist <-
    t(aggregate(y, "cum_diag_a", stratum = idx$X))

  #  asymptomatic diagnoses across all V+W
  cumulative_flows$cum_diag_a_VW_all_diaghist <-
    t(aggregate(y, "cum_diag_a", stratum = c(idx$V, idx$W)))

  #  asymptomatic diagnoses in first diagnosis history strata of X only
  cumulative_flows$cum_diag_a_X_first_diag_hist <-
    t(aggregate(y, "cum_diag_a", stratum = idx$X[1]))

  #  asymptomatic diagnoses in first diagnosis history strata of V+W only
  cumulative_flows$cum_diag_a_VW_first_diag_hist <-
    t(aggregate(y, "cum_diag_a", stratum = idx$never_diag))

  #  symptomatic diagnoses across all X
  cumulative_flows$cum_diag_s_X_all_diaghist <-
    t(aggregate(y, "cum_diag_s", stratum = idx$X))

  #  symptomatic diagnoses across all V+W
  cumulative_flows$cum_diag_s_VW_all_diaghist <-
    t(aggregate(y, "cum_diag_s", stratum = c(idx$V, idx$W)))

  #  symptomatic diagnoses in first diagnosis history strata of X only
  cumulative_flows$cum_diag_s_X_first_diag_hist <-
    t(aggregate(y, "cum_diag_s", stratum = idx$X[1]))

  #  symptomatic diagnoses in first diagnosis history strata of V+W only
  cumulative_flows$cum_diag_s_VW_first_diag_hist <-
    t(aggregate(y, "cum_diag_s", stratum = idx$never_diag))

  # Now we want person years of exposure in each timestep
  # Person-years of exposure before first diagnosis
  ## firstly in X (note this isn't cumulative)

  # total number in X.I, and V.I&W.I at each time point
  N_diag_hist_X <- t(aggregate(y, "N", stratum = idx$X[1]))
  N_diag_hist_VW <- t(aggregate(y, "N", stratum = idx$never_diag))

  # number in X.I and V.I&W.I undergoing treatment i.e people 'exposed' after
  # diagnosis
  N_diag_hist_X_T <- t(aggregate(y, "T", stratum = idx$X[1]))
  N_diag_hist_VW_T <- t(aggregate(y, "T", stratum = idx$never_diag))

  # obtain person-years exposed before first diagnosis i.e person-years in
  # treatment don't count towards this
  N_person_years_exposed_X <- N_diag_hist_X - N_diag_hist_X_T
  N_person_years_exposed_VW <- N_diag_hist_VW - N_diag_hist_VW_T
  N_pr_yrs_exp_all <- list(N_person_yrs_exp_X.I = N_person_years_exposed_X,
                           N_person_yrs_exp_VW.I = N_person_years_exposed_VW)
  #remove time = 0
  N_pr_yrs_exp_all <- lapply(N_pr_yrs_exp_all, function(x) x[-1, ])

  #calculate cumulative person-years exposed over time
  cum_N_pyrs <- lapply(N_pr_yrs_exp_all, function(x) apply(x, 2, cumsum))

  #extract annual flows
  flows <- lapply(cumulative_flows, function(x) apply(x, 2, diff))
  names(flows) <- gsub("^cum_", "", names(cumulative_flows))
  flows$N <- NULL

  # remove time 0 from cumulative flows
  cumulative_flows <- lapply(cumulative_flows, function(x) x[-1, ])

  ret <- c(cumulative_flows, flows, cum_N_pyrs)

  ret
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
                   "cum_treated", "cum_screened", "cum_vaccinated",
                   "cum_offered", "cum_vbe")

  for (nm in state_names) {
    dimnames(res[[nm]]) <- list(NULL, group_names, strata_names)
  }
  dimnames(res$lambda) <- list(NULL, group_names)
  dimnames(res$eta) <- list(NULL, group_names)

  res
}
