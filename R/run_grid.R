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
run_grid  <- function(t, gono_params, init_params, 
                      baseline, model,
                      eff, dur, ve = 0, strategy = "VbE", uptake_total = 0,
                      uptake_second_dose = uptake_total,
                      qaly_loss_per_diag_s,
                      unit_cost_manage_symptomatic,
                      unit_cost_manage_asymptomatic,
                      unit_cost_screen_uninfected,
                      t_stop = 99, full_output = FALSE, disc_rate = 0) {
  
  
  tt <- seq.int(init_params[[1]]$t, length.out = t + 1)
  l <- expand.grid(eff = eff, dur = dur)
  verify_baseline(baseline, l, gono_params, tt)
  
  
  # set strategy
  prop_vax <- lapply(uptake_total, set_strategy, strategy = strategy)
  
  vd <- lapply(prop_vax, "[[", "vd")
  vs <- lapply(prop_vax, "[[", "vs")
  
  # run model grid
  res <- furrr::future_pmap(.l = l, .f = model, tt = tt,
                            gono_params = gono_params,
                            init_params = init_params,
                            ve = ve, vd = vd, vs = vs,
                            t_stop = t_stop)
  
  # compare to baseline
  ret <-
    furrr::future_pmap(.l = list(y = res, baseline = baseline),
                       .f = compare_baseline,
                       uptake_second_dose = uptake_second_dose,
                       disc_rate = disc_rate,
                       qaly_loss_per_diag_s = qaly_loss_per_diag_s,
                       unit_cost_manage_symptomatic = unit_cost_manage_symptomatic,
                       unit_cost_manage_asymptomatic = unit_cost_manage_asymptomatic,
                       unit_cost_screen_uninfected = unit_cost_screen_uninfected)
  
  names(ret) <- sprintf("eff%.2f_dur%02d", l$eff, l$dur)
  
  # prepare results
  out <- list(inputs = list(t = tt, ve = ve, vd = vd, vs = vs, grid = l,
                            gono_params = gono_params,
                            init_params = init_params,
                            uptake = list(total = uptake_total,
                                          second_dose = uptake_second_dose)),
              results = ret)
  if (full_output) out$full_results <- res
  class(out) <- "gonovax_grid"
  
  out
}


verify_baseline <- function(baseline, l, gono_params, tt) {
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

compare_baseline <- function (y, baseline, disc_rate, uptake_second_dose,
                              qaly_loss_per_diag_s,
                              unit_cost_manage_symptomatic,
                              unit_cost_manage_asymptomatic,
                              unit_cost_screen_uninfected) {
  
  ## compare run to baseline
  ret <- gonovax:::extract_flows(y)
  ret_vs_baseline <- Map(`-`, ret, baseline[names(ret)])
  names(ret_vs_baseline) <- paste0("inc_", names(ret))
  
  ## calculate cases averted per dose, both with and without discounting
  ret$inc_doses <- calc_doses(ret_vs_baseline, uptake_second_dose, TRUE)
  ret$cases_averted_per_dose <-
    calc_cases_averted_per_dose(ret_vs_baseline, ret$inc_doses, 0)
  ret$cases_averted_per_dose_pv <-
    calc_cases_averted_per_dose(ret_vs_baseline, ret$inc_doses, disc_rate)
  
  ## calculate incremental net benefit per person with QALY cost - £20k / £30k
  ## both calcs allow for discounting (i.e. are present values as at 2022)
  
  ret$inc_net_benefit_pp_20k <-
    calc_incremental_net_benefit_per_person(cost_qaly = 2e5, ret_vs_baseline,
                                            ret$inc_doses, disc_rate,
                                            qaly_loss_per_diag_s,
                                            unit_cost_manage_symptomatic,
                                            unit_cost_manage_asymptomatic,
                                            unit_cost_screen_uninfected)
  ret$inc_net_benefit_pp_30k <-
    calc_incremental_net_benefit_per_person(cost_qaly = 3e5, ret_vs_baseline,
                                            ret$inc_doses, disc_rate,
                                            qaly_loss_per_diag_s,
                                            unit_cost_manage_symptomatic,
                                            unit_cost_manage_asymptomatic,
                                            unit_cost_screen_uninfected)
  c(ret, ret_vs_baseline)
}


