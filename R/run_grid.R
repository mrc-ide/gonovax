##' @name run_grid
##' @title run model from equilibrium with single vaccine at the input
##' efficacy / duration grid locations for n parameter sets
##' @param n a value denoting the number of parameter sets to run
##' @param t an integer number of years at which impact is to be assessed
##' @param model a gonovax run function, by default `run_onevax`
##' @param eff numeric vector (between 0-1) of efficacy values of the vaccine
##' @param dur numeric vector of duration in years of the vaccine
##' @param ve single numeric indicating % of population vaccinated before entry
##'  (between 0-1)
##' @param strategy ve: vaccination before entry only (default),
##' vd: vaccination on diagnosis, va: vaccination on attendance,
##' vt: targeted vaccination (i.e. all diagnosed, group H on screening)
##' @param uptake numeric (0-1) of strategy uptake
##' @param baseline optional input of a baseline to compare to, must be a
##' gonovax_grid object if supplied
##' @param full_output logical indicating whether full results should be output
##' @return A `gonovax_grid` object
##' @import furrr
##' @export
run_grid  <- function(n, t, model = run_onevax,
                      eff, dur, ve = 0, strategy = "ve", uptake = 0,
                      baseline = NULL, full_output = FALSE) {

  nn <- seq_len(n)
  tt <- seq(0, t)
  l <- expand.grid(eff = eff, dur = dur)

  # verify baseline

  if (is.null(baseline)) {
    baseline <- novax_baseline(nn, tt[-1])
    baseline$vaccinated <- baseline$cum_vaccinated <- 0
    baseline <- rep(list(baseline), nrow(l))
  } else {
    baseline <- verify_baseline(baseline, l, nn, tt)
    baseline <- switch_levels(baseline[-1])
  }

  # set strategy
  prop_vax <- set_strategy(strategy, uptake)

  res <- furrr::future_pmap(.l = l,
                            .f = model,
                            n = nn,
                            tt = tt,
                            ve = ve,
                            vd = prop_vax$vd,
                            vs = prop_vax$vs,
                            equilib = TRUE)

  ret <- furrr::future_pmap(.l = list(y = res, baseline = baseline),
                            .f = compare_baseline)

  names(ret) <- sprintf("eff%.2f_dur%02d", l$eff, l$dur)
  ret <- switch_levels(ret)
  out <- c(list(inputs = list(n = nn, t = tt, ve = ve,
                            vd = prop_vax$vd, vs = prop_vax$vs, grid = l)),
            ret)
  if (full_output) out$full_results <- res
  class(out) <- "gonovax_grid"
  out
}

verify_baseline <- function(baseline, l, nn, tt) {
  if (!inherits(baseline, "gonovax_grid")) {
    stop("baseline must be a gonovax_grid object")
  }
  if (!identical(baseline$inputs$grid, l)) {
    stop("dur / eff parameters do not match baseline")
  }
  if (!identical(baseline$inputs$n, nn)) {
    stop("model parameters do not match baseline")
  }
  if (!identical(baseline$inputs$t, tt)) {
    stop("t does not match baseline")
  }
  baseline
}

compare_baseline <- function(y, baseline) {

  # extract cumulative flows of interest
  cum_incid      <- t(aggregate(y, "cum_incid"))[-1, ]
  cum_vaccinated <- t(aggregate(y, "cum_vaccinated"))[-1, ]
  cum_diag_a     <- t(aggregate(y, "cum_diag_a"))[-1, ]
  cum_diag_s     <- t(aggregate(y, "cum_diag_s"))[-1, ]
  cum_screened   <- t(aggregate(y, "cum_screened"))[-1, ]

  # extract flows of interest
  incid          <- t(aggregate(y, "cum_incid", TRUE))
  vaccinated     <- t(aggregate(y, "cum_vaccinated", TRUE))

  # compare to baseline
  inc_vaccinated <- vaccinated - baseline$vaccinated
  red_incid <- baseline$incid - incid


  # output results
  list(incid = incid,
       vaccinated = vaccinated,
       cum_incid = cum_incid,
       cum_vaccinated = cum_vaccinated,
       cum_diag_a = cum_diag_a,
       cum_diag_s = cum_diag_s,
       cum_screened = cum_screened,
       red_incid = red_incid,
       inc_vaccinated = inc_vaccinated,
       red_cum_incid = baseline$cum_incid - cum_incid,
       inc_cum_vaccinated = cum_vaccinated - baseline$cum_vaccinated,
       red_cum_diag_a = baseline$cum_diag_a - cum_diag_a,
       red_cum_diag_s = baseline$cum_diag_s - cum_diag_s,
       inc_cum_screened = cum_screened - baseline$cum_screened
  )
}


##' @name format_grid
##' @title format grid for heatmap plotting
##' @param grid a `gonovax_grid` object
##' @param disc_rate annual discount rate for cost-effectiveness calc,
##' e.g. 0.035 = 3.5pc, default is 0
##' @param f the function that should be used to summarise the runs
##' @return a dataframe with columns denoting:
##' eff: efficacy of vaccine (%)
##' dur: duration of vaccine (years)
##' red_incid: Reduction in incidence after t years,
##' inc_vaccinated: Additional courses of vaccine over t years,
##' tot_red_incid: Infections averted over t years,
##' cost_eff: Courses of vaccine per infection averted (B / C),
##' discounted at `disc_rate`
##' tot_red_diag_a: Reduction in asymptomatic diagnoses over t years
##' tot_red_diag_s: Reduction in symptomatic diagnoses over t years
##' tot_inc_screened: Increase in screening over t years
##' @export
format_grid <- function(grid, disc_rate = 0, f = mean) {

  stopifnot(inherits(grid, "gonovax_grid"))

  tt <- grid$inputs$t[-1]

  # discount flows by (1 + i) ^ -(t - dt/2) to find PV

  dt <- diff(tt)[1]
  pv <- (1 + disc_rate) ^ - (tt - dt / 2)

  # calculate cost-effectiveness
  pv_inc_cum_vaccinated <- calc_pv(grid$inc_vaccinated, pv)

  pv_red_cum_incid <- calc_pv(grid$red_incid, pv)
  cost_eff <- mapply(FUN = `/`, pv_inc_cum_vaccinated, pv_red_cum_incid,
                     SIMPLIFY = FALSE)

  list(t                  = tt,
       eff                = grid$inputs$grid$eff * 100,
       dur                = grid$inputs$grid$dur,
       red_incid          = summarise(grid$red_incid, f),
       tot_inc_vaccinated = summarise(grid$inc_cum_vaccinated, f),
       tot_red_incid      = summarise(grid$red_cum_incid, f),
       cost_eff           = summarise(cost_eff, f),
       tot_red_diag_a     = summarise(grid$red_cum_diag_a, f),
       tot_red_diag_s     = summarise(grid$red_cum_diag_s, f),
       tot_inc_screened   = summarise(grid$inc_cum_screened, f))
}

calc_pv <- function(x, pv) {
  lapply(x, function(x) apply(x * pv, 2, cumsum))
}

summarise <- function(x, f) {
  y <- lapply(x, function(x) apply(x, 1, f))
  as.data.frame(y)
}
