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

  strata <- dimnames(y[[1]]$N)[[3]] # extract strata names
  n_erlang <- sum((grepl("V", strata))) # count erlang

  idx <- stratum_index_xpvwrh(n_erlang)

  ## extract vaccinations and revaccinations separately
  # primary vaccinated = everyone vaccinated from X(1) regardless of # doses
  cumulative_flows$cum_primary_total <-
    t(aggregate(y, "cum_vaccinated", stratum = idx$X))

  # partial to full = everyone vaccinated from P(2)
  # in earlier models this will be V, we expect this to be 0
  cumulative_flows$cum_part_to_full <-
    t(aggregate(y, "cum_vaccinated", stratum = idx$P))

  # revaccinated = everyone vaccinated from W(4), note does not include
  # vaccination of individuals in X who have waned from P
  cumulative_flows$cum_revaccinated <-
    t(aggregate(y, "cum_vaccinated", stratum = idx$W))

  ######################
  #X = 1, P = 2, V = 3, W = 4, R = 5, H = 6, (for n_erlang = 1)
  # extract asymptomatic and symptomatic diagnoses in vaccinated (P,V,R)
  # and unvaccinated strata separately (X,W)

  #asymptomatic, unvax
  cumulative_flows$cum_diag_a_xwh <-
    t(aggregate(y, "cum_diag_a", stratum = c(idx$X, idx$W, idx$H)))

  #symptomatic, unvax
  cumulative_flows$cum_diag_s_xwh <-
    t(aggregate(y, "cum_diag_s", stratum = c(idx$X, idx$W, idx$H)))

  #asymptomatic, all-vax
  cumulative_flows$cum_diag_a_pvr <-
    t(aggregate(y, "cum_diag_a", stratum = c(idx$P, idx$V, idx$R)))

  #symptomatic, all-vax
  cumulative_flows$cum_diag_s_pvr <-
    t(aggregate(y, "cum_diag_s", stratum = c(idx$P, idx$V, idx$R)))

  #total, unvax (cumulative treated will be a bit different to the sum of a & s)
  # unvax
  cumulative_flows$cum_diag_t_xwh <-
    cumulative_flows$cum_diag_a_xwh + cumulative_flows$cum_diag_s_xwh

  # all-vax
  cumulative_flows$cum_diag_t_pvr <-
    cumulative_flows$cum_diag_a_pvr + cumulative_flows$cum_diag_s_pvr

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


##' @name adjust_baseline
##' @title adjust model run in absence of vaccination so diagnoses are
##' spread across vaccine protected and non-vaccine protected strata as if
##' the rate of movement was the same as in an equivalent model run in the
##' presence of vaccination
##' @param baseline A model run in the absence of vaccine uptake
##' @param y A model run in the presence of vaccine uptake
##' @return An adjusted baseline in which 'a' and 's' diagnoses in the baseline
##' are divided across vaccine protected and non-vaccine protected strata
##' allowing comparison of model runs with of baselines by vaccine-status
##' @export
adjust_baseline <- function(baseline, y) {
  Map(adjust_baseline_one, baseline, y)
}

adjust_baseline_one <- function(baseline, y) {
  # create an array in the same dimensions as y$N, where each entry is the
  # modelled proportion in each stratum for each risk group at time t
  # i.e. summing across strata will give 1: apply(prop_N, c(1, 2), sum)
  prop_N <-  y$N / array(apply(y$N, c(1, 2), sum), dim(y$N))

  # identify states that are compartments [time, group, strata]
  idx_state <- which(lengths(lapply(baseline, dim)) == 3)
  # distribute individuals in the baseline across strata in the same proportion
  # as in the model run
  adjusted_states <- lapply(baseline[idx_state],
                            function(state) c(state[, , 1]) * prop_N)
  replace(baseline, idx_state, adjusted_states)

}

##' @name gen_labels
##' @title generates the appropriate strata labels for the number of strata
##' in the model, which depends on the value given to n_erlang and diagnosis
##' history levels desired (n_diag_rec)
##' @param n_erlang integer giving the number of transitions that need to be
##'  made through vaccine-protected strata until that protection has waned
##' @param n_diag_rec integer giving the number of levels of diagnosis history
##' for each X, V(*n_erlang), and W stratum
##' @return a character vector of length n_vax containing strata labels
##' @export
gen_labels <- function(n_erlang = 1, n_diag_rec = 1) {


  diag_hist <- paste0(".", as.roman(seq_len(n_diag_rec)))

  output <- c(paste0("X", diag_hist),
              paste0("V", rep(seq_len(n_erlang), each = n_diag_rec), diag_hist),
              paste0("W", diag_hist))

  output


}