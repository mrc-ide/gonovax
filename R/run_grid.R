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

  # verify baseline
  nn <- seq_len(n)
  if (is.null(baseline)) {
    baseline <- novax_baseline(nn, seq(1, t))
    baseline$cum_vaccinated <- 0
  } else {
    baseline <- verify_baseline(baseline, l, nn, t)
  }
  
  l <- expand.grid(eff = eff, dur = dur)

  # set strategy
  prop_vax <- set_strategy(strategy, uptake)
  
  res <- furrr::future_pmap(.l = l,
                            .f = model,
                            n = nn,
                            tt = seq(0, t),
                            ve = ve,
                            vd = prop_vax$vd,
                            vs = prop_vax$vs,
                            equilib = TRUE)
  compare_baseline(res[[1]], baseline)
  ret <- furrr::future_pmap(.l = list(y = res),
                            .f = compare_baseline, 
                            baseline = baseline)

  out <- list(inputs = list(n = nn, t = t, ve = ve,
                            vd = prop_vax$vd, vs = prop_vax$vs, grid = l),
            ret)
  if (full_output) out$results <- res
  class(out) <- "gonovax_grid"
  out
}

verify_baseline <- function(baseline, l, nn, t) {
  if (!inherits(baseline, "gonovax_grid")) {
    stop("baseline must be a gonovax_grid object")
  }
  if (!identical(baseline$inputs$grid, l)) {
    stop("dur / eff parameters do not match baseline")
  }
  if (!identical(baseline$inputs$n, nn)) {
    stop("model parameters do not match baseline")
  }
  if (!identical(baseline$inputs$t, t)) {
    stop("t does not match baseline")
  }
  baseline
}

extract_value <- function(x, what, as_incid = FALSE) {
  i <- seq(2, length(x[[1]][[1]]$t)) - as_incid
  lapply(x, function(x) aggregate(x, what, as_incid)[, i])
}

##' @name format_grid
##' @title format grid for heatmap plotting
##' @param grid a `gonovax_grid` object
##' @return a dataframe with columns denoting:
##' eff: efficacy of vaccine (%)
##' dur: duration of vaccine (years)
##' a: Reduction in incidence after t years,
##' b: Courses of vaccine over t years,
##' c: Infections averted over t years,
##' d: Courses of vaccine per infection averted (B / C)
##' diag_a: Reduction in asymptomatic diagnoses over t years
##' diag_s: Reduction in symptomatic diagnoses over t years
##' @export
format_grid <- function(grid) {
  stopifnot(inherits(grid, "gonovax_grid"))
  data.frame(eff    = grid$inputs$grid$eff * 100,
             dur    = grid$inputs$grid$dur,
             a      = colMeans(grid$red_incid),
             b      = colMeans(grid$cum_vaccinated),
             c      = colMeans(grid$cum_red_incid),
             d      = colMeans(grid$cum_vaccinated / grid$cum_red_incid),
             diag_a = colMeans(grid$cum_red_diag_a),
             diag_s = colMeans(grid$cum_red_diag_s))
}

compare_baseline <- function(y, baseline) {
  cum_incid      <- t(aggregate(y, "cum_incid"))[-1, ]
  incid          <- aggregate(y, "cum_incid", TRUE)
  cum_vaccinated <- t(aggregate(y, "cum_vaccinated"))[-1, ]
  cum_diag_a     <- t(aggregate(y, "cum_diag_a"))[-1, ]
  cum_diag_s     <- t(aggregate(y, "cum_diag_s"))[-1, ]
  cum_screened   <- t(aggregate(y, "cum_screened"))[-1, ]
  
  list(incid = incid,
       cum_incid = cum_incid,
       cum_vaccinated = cum_vaccinated - baseline$cum_vaccinated,
       cum_diag_a = cum_diag_a,
       cum_diag_s = cum_diag_s,
       cum_screened = cum_screened,
       red_incid = baseline$incid - incid,
       red_cum_incid = baseline$cum_incid - cum_incid,
       red_cum_diag_a = baseline$cum_diag_a - cum_diag_a,
       red_cum_diag_s = baseline$cum_diag_s - cum_diag_s,
       red_cum_screened = baseline$cum_screened - cum_screened
       )
}
