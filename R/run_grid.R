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
##' @param t_stop time at which vaccination should stop (years)
##' @param baseline optional input of a baseline to compare to, must be a
##' gonovax_grid object if supplied
##' @param full_output logical indicating whether full results should be output
##' @return A `gonovax_grid` object
##' @import furrr
##' @export
run_grid  <- function(n, t, model = run_onevax,
                      eff, dur, ve = 0, strategy = "ve", uptake = 0,
                      t_stop = 99,
                      baseline = NULL, full_output = FALSE) {

  nn <- seq_len(n)
  tt <- seq(0, t)
  l <- expand.grid(eff = eff, dur = dur)

  # verify baseline

  if (is.null(baseline)) {
    baseline <- novax_baseline(nn, tt[-1])
    baseline <- rep(list(baseline), nrow(l))
  } else {
    baseline <- verify_baseline(baseline, l, nn, tt)
    baseline <- baseline[-1]
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
                            t_stop = t_stop,
                            equilib = TRUE)

  ret <- furrr::future_pmap(.l = list(y = res, baseline = baseline),
                            .f = compare_baseline)

  names(ret) <- sprintf("eff%.2f_dur%02d", l$eff, l$dur)

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

  # extract cumulative and annula flows
  ret <- extract_flows(y)

  # compare to baseline (most will be reductions)
  ret_vs_baseline <- Map(`-`, ret, baseline[names(ret)])
  names(ret_vs_baseline) <- paste0("inc_", names(ret))

  # output results
  c(ret, ret_vs_baseline)
}
