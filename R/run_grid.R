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
run_grid  <- function(t, gono_params, init_params, cost_params,
                      baseline, model,
                      eff, dur, ve = 0, strategy = "VbE", uptake_total = 0,
                      uptake_second_dose = uptake_total,
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
  ret <- furrr::future_pmap(.l = list(y = res, baseline = baseline),
                       .f = compare_baseline,
                       cost_params = cost_params,
                       uptake_second_dose = uptake_second_dose,
                       disc_rate = disc_rate)
  
  names(ret) <- sprintf("eff%.2f_dur%02d", l$eff, l$dur)
  
  # prepare results
  out <- list(inputs = list(t = tt, ve = ve, vd = vd, vs = vs, grid = l,
                            gono_params = gono_params,
                            init_params = init_params,
                            cost_params = cost_params,
                            disc_rate = disc_rate,
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

compare_baseline <- function(y, baseline, uptake_second_dose, cost_params,
                              disc_rate) {
  
  ## compare run to baseline
  ret <- extract_flows(y)
  ret_vs_baseline <- Map(`-`, ret, baseline[names(ret)])
  names(ret_vs_baseline) <- paste0("inc_", names(ret))
  ret <- c(ret, ret_vs_baseline)
  
  ## calculate cases averted per dose, both with and without discounting
  ret$inc_doses <- calc_doses(ret, uptake_second_dose, TRUE)
  ret$cases_averted_per_dose <- calc_cases_averted_per_dose(ret, 0)
  ret$cases_averted_per_dose_pv <- calc_cases_averted_per_dose(ret, disc_rate)
  
  ## calculate cost-effectiveness threshold price as incremental net benefit per
  ## person with QALY cost - £20k / £30k
  ## both calcs allow for discounting (i.e. are present values as at 2022)
  
  ret$cet_20k <- calc_cost_eff_threshold(2e5, ret, cost_params, disc_rate)
  ret$cet_30k <- calc_cost_eff_threshold(3e5, ret, cost_params, disc_rate)

  ret
}

calc_doses <- function(forecast, uptake_second_dose, revax_one_dose = TRUE) {
  
  n_primary_doses_pp <- 1 + 1 / uptake_second_dose
  
  if (revax_one_dose) {
    n_revax_doses_pp <- 1
  } else {
    n_revax_doses_pp <- n_primary_doses_pp
  }
  
  inc_primary_vaccination <- forecast$inc_vaccinated - forecast$inc_revaccinated
  inc_primary_doses <- t(inc_primary_vaccination) * n_primary_doses_pp
  inc_revax_doses <- t(forecast$inc_revaccinated) * n_revax_doses_pp
  
  inc_doses <- inc_primary_doses + inc_revax_doses
  
  t(inc_doses)
}

calc_pv <- function(x, disc_rate) {
  
  tt <- seq_len(nrow(x))
  pv <- (1 + disc_rate) ^ -(tt - 0.5)
  ret <- apply(x * pv, 2, cumsum)
  ret
}

calc_cases_averted_per_dose <- function(forecast, disc_rate = 0) {
  
  pv_inc_doses <- calc_pv(forecast$inc_doses, disc_rate)
  pv_cases_averted  <- calc_pv(-forecast$inc_treated, disc_rate)
  
  cases_averted_per_dose <- pv_cases_averted / pv_inc_doses
  cases_averted_per_dose[is.na(cases_averted_per_dose)] <- 0
  
  cases_averted_per_dose
}


calc_cost_eff_threshold <-
function(cost_qaly, forecast, cost_params, disc_rate) {
    
    qaly_cost_per_diag_s <- cost_params$qaly_loss_per_diag_s * cost_qaly
    
    red_diag_s <- -t(forecast$inc_diag_s)
    red_diag_a <- -t(forecast$inc_diag_a)
    inc_screen <- t(forecast$inc_screened)
    
    red_cost_diag_s <- red_diag_s * cost_params$unit_cost_manage_symptomatic
    red_cost_diag_a <- red_diag_a * cost_params$unit_cost_manage_asymptomatic
    inc_cost_screen <- inc_screen * cost_params$unit_cost_screen_uninfected
    red_qaly_cost   <- red_diag_s * qaly_cost_per_diag_s
    
    inc_net_benefit <-
      t(red_cost_diag_s + red_cost_diag_a - inc_cost_screen + red_qaly_cost)
    
    pv_inc_doses <- calc_pv(forecast$inc_doses, disc_rate)
    pv_inc_net_benefit <- calc_pv(inc_net_benefit, disc_rate)
    
    ## calc pv incremental net benefit per person protected
    pv_inc_net_benefit / pv_inc_doses
    
  }
