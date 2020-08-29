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
  l <- expand.grid(eff = eff, dur = dur)
  nn <- seq_len(n)
  # set strategy
  prop_vax <- set_strategy(strategy, uptake)

  res <- furrr::future_pmap(.l = l,
                            .f = model,
                            n = nn,
                            tt = c(0, t - 1, t),
                            ve = ve,
                            vd = prop_vax$vd,
                            vs = prop_vax$vs,
                            equilib = TRUE)

  cum_incid <- extract_value(res, "cum_incid", t)
  incid <- cum_incid - extract_value(res, "cum_incid", t - 1)
  cum_vaccinated <- extract_value(res, "cum_vaccinated", t)
  cum_diag_a <- extract_value(res, "cum_diag_a", t)
  cum_diag_s <- extract_value(res, "cum_diag_s", t)

  # verify baseline

  if (is.null(baseline)) {
    baseline <- novax_baseline(nn, t)
    baseline$cum_vaccinated <- 0
  } else {
    baseline <- verify_baseline(baseline, l, nn, t)
  }

  out <- list(inputs = list(n = nn, t = t, ve = ve,
                            vd = prop_vax$vd, vs = prop_vax$vs, grid = l),
              incid = incid,
              cum_incid = cum_incid,
              cum_vaccinated = cum_vaccinated - baseline$cum_vaccinated,
              red_incid = baseline$incid - incid,
              cum_red_incid = baseline$cum_incid - cum_incid,
              cum_diag_a = cum_diag_a,
              cum_diag_s = cum_diag_s,
              cum_red_diag_a = baseline$cum_diag_a - cum_diag_a,
              cum_red_diag_s = baseline$cum_diag_s - cum_diag_s)
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

extract_value <- function(x, what, t) {
  sapply(x, function(x) aggregate(x, what)[, x[[1]]$t == t])
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
