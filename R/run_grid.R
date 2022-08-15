##' @name run_grid
##' @title run model from equilibrium with single vaccine at the input
##' efficacy / duration grid locations for n parameter sets
##' @param gono_params gono params
##' @param init_params initial state parameters
##' @param cost_params cost effectiveness parameters
##' @param baseline optional input of a baseline to compare to, must be a
##' gonovax_grid object if supplied
##' @param model a gonovax run function, by default `run_onevax`
##' @param eff numeric vector (between 0-1) of vaccine efficacy against
##' acquisition
##' @param dur numeric vector of duration in years of the vaccine
##' @param vbe single numeric indicating % of population vaccinated before entry
##'  (between 0-1), default = 0
##' @param strategy, default is null, no vaccination
##' `VoD`: vaccination on diagnosis, `VoA`: vaccination on attendance,
##' `VoD(L)+VoA(H)`: targeted vaccination (i.e. all diagnosed plus group H on
##'  screening)
##' @param uptake_total numeric (0-1) of strategy uptake, default = 0
##' @param uptake_second_dose numeric (0-1) of strategy uptake,
##'  default = uptake_total
##' @param t_stop time at which vaccination should stop (years), default = 99
##' @param full_output logical indicating whether full results should be output
##' @param disc_rate discount rate for cost-effectiveness output
##'  per annum, default = 0
##' @return A `gonovax_grid` object
##' @import furrr
##' @export
run_grid  <- function(gono_params, init_params, cost_params,
                      baseline, model,
                      eff, dur, vbe = 0, strategy = NULL, uptake_total = 0,
                      uptake_second_dose = uptake_total,
                      t_stop = 99, full_output = FALSE, disc_rate = 0) {

  tt <- seq.int(init_params[[1]]$t,
                length.out = nrow((baseline[[1]][[1]])) + 1)
  l <- expand.grid(vea = eff, dur = dur)

  # run model grid
  res <- furrr::future_pmap(.l = l, .f = model, tt = tt,
                            gono_params = gono_params,
                            init_params = init_params,
                            vbe = vbe, uptake = uptake_total,
                            strategy = strategy,
                            t_stop = t_stop)

  # compare to baseline
  ret <- furrr::future_pmap(.l = list(y = res, baseline = baseline),
                            .f = compare_baseline,
                            cost_params = cost_params,
                          uptake_first_dose = uptake_total / uptake_second_dose,
                            uptake_second_dose = uptake_second_dose,
                            disc_rate = disc_rate)

  names(ret) <- sprintf("eff%.2f_dur%02d", l$vea, l$dur)

  # prepare results
  out <- list(inputs = list(t = tt, vbe = vbe, strategy = strategy, grid = l,
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

##' @name compare_baseline_xpvwrh
##' @title compare model runs with vaccination to a baseline run for the
##' branching XPVWRH model where both partial and full vaccination have a
##' given level of efficacy
##' @param y list of model runs, e.g. created by `run_onevax_xvwv`, each list
##' entry refers to a different parameter set
##' @param baseline list of baseline runs e.g. created by `run_onevax_xvwv`.
##' Should be same length as `y`
##' @param cost_params list of cost effectiveness parameters, containing entries
##' `qaly_loss_per_diag_s`, `unit_cost_manage_symptomatic`,
##' `unit_cost_manage_asymptomatic`, `unit_cost_screen_uninfected`. Each entry
##' length 1 or same as `y`
##' @param uptake_first_dose numeric (0-1) either length 1 or same as `y`
##' @param uptake_second_dose numeric (0-1) either length 1 or same as `y`
##' @param disc_rate discount rate for cost-effectiveness output per annum,
##'  default = 0
##' @return A list of matrices with rows = time point and columns
##' = parameter set
##' list entries are:
##' `cum_treated` = cumulative number treated in model run
##' `cum_vaccinated` = cumulative vaccinated treated in model run
##' `diag_a` = annual number of asymptomatic diagnoses in model run
##' `diag_s` = annual number of symptomatic diagnoses in model run
##' `treated` = annual number treated in model run
##' `screened` = annual number screened in model run
##' `vaccinated` = annual number vaccinated in model run
##' (includes those getting vbe and boosters)
##' `vbe` = annual number vaccinated before entry in model run
##' `revaccinated` = annual number receiving booster in model run
##' `inc_cum_treated` = cumulative number treated compared to baseline
##' `inc_cum_vaccinated` = cumulative number vaccinated compared to baseline
##' `inc_diag_a` = annual number of asymptomatic diagnoses compared to baseline
##' `inc_diag_s` = annual number of symptomatic diagnoses compared to baseline
##' `inc_treated` = annual number treated compared to baseline
##' `inc_screened` = annual number screened compared to baseline
##' `inc_vaccinated`  = annual number vaccinated compared to baseline
##' (includes those getting vbe and boosters)
##' `inc_vbe` = annual number vaccinated before entry compared to baseline
##' `inc_revaccinated` = annual number receiving booster compared to
##' baseline
##' `inc_cum_revaccinated` = cumulative number revaccinated by recieving booster
##' compared to baseline
##' `inc_primary` = annual number receiving primary vaccination (full or
##' partial) compared to baseline (does not include vbe)
##' `inc_cum_primary` = cumulative number of individuals receiving primary
##' vaccination (full or partial) compared to baseline (does not include vbe)
##' `inc_part_to_full`= annual number of individuals who have already
##' received their 1st dose, who receive their 2nd dose compared to baseline
##' `inc_cum_part_to_full` = cumulative number of individuals who have already
##' received their 1st dose, who receive their 2nd dose compared to baseline
##' `inc_doses` = number of doses compared to baseline (per year). Assumes
##' primary vaccination uses 2 doses, booster uses 1 dose.
##' `inc_cum_doses` = cumulative number of doses compared to baseline. Assumes
##' primary vaccination uses 2 doses, booster uses 1 dose.
##' `inc_primary_part_doses` = annual number of 1st doses administered as part
##' of primary vaccination which were not followed by a 2nd dose compared to
##' baseline.
##' `inc_cum_primary_part_doses` = cumulative number of 1st doses
##' administered as part of primary vaccination which were not followed by a 2nd
##' dose compared to baseline.
##' `inc_primary_full_doses` = annual number 1st and 2nd doses administered
##' together as part of primary vaccination compared to baseline.
##' `inc_cum_primary_full_doses` = cumulative number of 1st and 2nd doses
##' administered together as part of primary vaccination compared to baseline.
##' `inc_primary_total_doses` = annual number of primary doses administered
##' compared to baseline. Includes doses contributing to partial and full
##' vaccination
##' `inc_cum_primary_total_doses` = cumulative number of primary doses
##' administered compared to baseline. Includes doses contributing to partial
##' and full vaccination.
##' `inc_booster_doses` = annual number of booster doses
##' administered to waned individuals for revaccination compared to baseline,
##' where 1 dose is needed for full protection.
##' `inc_cum_booster_doses` = cumulative number of booster doses
##' administered to waned individuals for revaccination compared to baseline,
##'  where 1 dose is needed for full protection.
##'  `inc_vbe_doses` = annual number of vaccinations given before entry,
##'  where 2 doses give full protection.
##'  `inc_cum_vbe_doses` = cumulative number of vaccinations given before entry,
##'  where 2 doses give full protection.
##' `cases_averted_per_dose` = cumulative number of cases (i.e. diagnoses)
##' averted per dose of vaccine
##' `cases_averted_per_dose_pv` = present value of cases_averted_per_dose
##' (i.e. sum of annual numbers discounted at rate `disc_rate` to time 0)
##' `pv_inc_doses` = present value of number of doses compared to baseline
##' (i.e. sum of annual numbers discounted at rate `disc_rate` to time 0)
##' `pv_red_net_cost` = present value of the reduction in net costs compared to
##' baseline (i.e. sum of annual reduction in costs discounted at rate
##'  `disc_rate` to time 0)
##' `pv_qaly_gain` = present value of the QALYs gained due to vaccination
##' compared to baseline (i.e. QALYs gained each year discounted at rate
##'  `disc_rate` to time 0, and summed)
##' `cet_20k` = cost effectiveness threshold (i.e. £ value per dose at which
##' vaccination becomes cost-effective) calculated using £20,000 / QALY.
##' (pv_red_net_cost + pv_qaly_gain * £20,000) / pv_inc_doses
##' `cet_30k` = cost effectiveness threshold (i.e. £ value per dose at which
##' vaccination becomes cost-effective) calculated using £30,000 / QALY
##' (pv_red_net_cost + pv_qaly_gain * £30,000) / pv_inc_doses
##' `inc_costs_9` = present value of incremental costs assuming £9 / dose.
##' Incremental costs are calculated as: pv_inc_doses * £9 - pv_red_net_cost
##' `inc_costs_18` = present value of incremental costs assuming £18 / dose.
##' Incremental costs are calculated as: pv_inc_doses * £18 - pv_red_net_cost
##' `inc_costs_35` = present value of incremental costs assuming £35 / dose.
##' Incremental costs are calculated as: pv_inc_doses * £35 - pv_red_net_cost
##' `inc_costs_50` = present value of incremental costs assuming £50 / dose.
##' Incremental costs are calculated as: pv_inc_doses * £50 - pv_red_net_cost
##' `inc_costs_70` = present value of incremental costs assuming £70 / dose.
##' Incremental costs are calculated as: pv_inc_doses * £70 - pv_red_net_cost
##' `inc_costs_85` = present value of incremental costs assuming £85 / dose
##' Incremental costs are calculated as: pv_inc_doses * £85 - pv_red_net_cost
##' `vacprotec_full` = number of fully vaccine protected individuals in the
##'  population at a given timepoint
##' `vacprotec_part` = number of partially vaccine protected individuals in the
##'  population at a given timepoint
##' `vacprotec_total` = total number of vaccine protected individuals in the
##'  population (both partially and fully vaccinated) at a given timepoint
##' `vacprotec_full_prop` = proportion of the population experiencing full
##' vaccine protection at a given timepoint
##' `vacprotec_part_prop` = proportion of the population experiencing partial
##' vaccine protection at a given timepoint
##' `vacprotec_total_prop` = total proportion of the population experiencing
##' some form of vaccine protection at a given timepoint
##' @export
compare_baseline_xpvwrh <- function(y, baseline, uptake_first_dose,
                             uptake_second_dose, cost_params,
                             disc_rate) {
  ## compare run to baseline
  flows <- extract_flows_xpvwrh(y)
  ret <- Map(`-`, flows, baseline[names(flows)])
  names(ret) <- paste0("inc_", names(flows))
  ret <- c(flows, ret)

  ## extract number under vaccine protection
      vacsnap <- list()
     # fully vaccine protected, snapshot of N in: V(3) and R(5)
        vacsnap$vacprotec_full <-
        t(aggregate(y, "N", stratum = c(3, 5)))

    # partially vaccine protected, snapshot of N in: P(2)
        vacsnap$vacprotec_part <-
         t(aggregate(y, "N", stratum = 2))

    # vaccine protected total, snapshot of N in: P(2), V(3), W(4)
        vacsnap$vacprotec_total <-
        t(aggregate(y, "N", stratum = c(2, 3, 5)))

    # remove the 0
        vacsnap <- lapply(vacsnap, "[", -1, )

  ret <- c(vacsnap, ret)

  ## calculate proportion of the population under vaccine protection
    # get total pop size  (including H!) This should be the same across
    # all model runs for all timepoints
    N <- t(aggregate(y, "N"))[1, 1]

    ret$vacprotec_full_prop <- ret$vacprotec_full / N
    ret$vacprotec_part_prop <- ret$vacprotec_part / N
    ret$vacprotec_total_prop <- ret$vacprotec_total / N

  ## calculate number receiving primary vaccination
  ret$inc_primary <- ret$inc_primary_total - ret$inc_vbe
  ret$inc_cum_primary <- apply(ret$inc_primary, 2, cumsum)

  ## calculate number receiving partial to full vaccination
  ret$inc_cum_part_to_full <- apply(ret$inc_part_to_full,
                                    2, cumsum)

  ## cumulative individuals receiving booster vaccination (revaccination)
  ret$inc_cum_revaccinated <- apply(ret$inc_revaccinated, 2, cumsum)

  ## calculate cases averted per dose, both with and without discounting
  ## return incremental annual and cumulative doses

  # all vbe get two doses
  ret$inc_vbe_doses <- ret$inc_vbe * 2

  # calculate doses given per person offered primary vaccination
  ret$inc_primary_full_doses <- ret$inc_primary *  uptake_second_dose * 2
  ret$inc_primary_part_doses <- ret$inc_primary * (1 - uptake_second_dose)
  ret$inc_primary_total_doses <- ret$inc_primary_full_doses +
    ret$inc_primary_part_doses

  # partial-to-full and booster vaccination take a single dose so is simply the
  # number of people vaccinated from P and W respectively
  ret$inc_part_to_full_doses <- ret$inc_part_to_full
  ret$inc_booster_doses <- ret$inc_revaccinated

  # calculate cumulative doses
  ret$inc_cum_vbe_doses <- apply(ret$inc_vbe_doses, 2, cumsum)
  ret$inc_cum_primary_full_doses <- apply(ret$inc_primary_full_doses, 2, cumsum)
  ret$inc_cum_primary_part_doses <- apply(ret$inc_primary_part_doses, 2, cumsum)
  ret$inc_cum_primary_total_doses <- apply(ret$inc_primary_total_doses, 2,
                                           cumsum)
  ret$inc_cum_part_to_full_doses <- apply(ret$inc_part_to_full_doses, 2, cumsum)
  ret$inc_cum_booster_doses <- apply(ret$inc_booster_doses, 2, cumsum)

  ret$inc_doses <- ret$inc_primary_total_doses +
                    ret$inc_part_to_full_doses +
                    ret$inc_booster_doses +
                    ret$inc_vbe_doses
  ret$inc_cum_doses <- apply(ret$inc_doses, 2, cumsum)

  ret$cases_averted_per_dose <- calc_cases_averted_per_dose(ret, 0)
  ret$cases_averted_per_dose_pv <- calc_cases_averted_per_dose(ret, disc_rate)

  ## calculate costs
  costs <- calc_costs(ret, cost_params, disc_rate)
  ret <- c(ret, costs)

  ## calculate cost-effectiveness threshold price as incremental net benefit per
  ## person with QALY cost - £20k / £30k
  ## both calcs allow for discounting (i.e. are present values as at 2022)
  ret$cet_20k <- calc_cet(2e4, costs)
  ret$cet_30k <- calc_cet(3e4, costs)

  ## calculate incremental cost of vaccination assuming
  ## £9, £18, £50, £85 per dose
  ## both calcs allow for discounting (i.e. are present values as at 2022)
  ret$inc_costs_9 <- calc_inc_costs(9, costs)
  ret$inc_costs_18 <- calc_inc_costs(18, costs)
  ret$inc_costs_35 <- calc_inc_costs(35, costs)
  ret$inc_costs_50 <- calc_inc_costs(50, costs)
  ret$inc_costs_70 <- calc_inc_costs(70, costs)
  ret$inc_costs_85 <- calc_inc_costs(85, costs)

  ret
}

##' @name compare_baseline
##' @title compare model runs with vaccination to a baseline runs
##' @param y list of model runs, e.g. created by `run_onevax_xvwv`, each list
##' entry refers to a different parameter set
##' @param baseline list of baseline runs e.g. created by `run_onevax_xvwv`.
##' Should be same length as `y`
##' @param cost_params list of cost effectiveness parameters, containing entries
##' `qaly_loss_per_diag_s`, `unit_cost_manage_symptomatic`,
##' `unit_cost_manage_asymptomatic`, `unit_cost_screen_uninfected`. Each entry
##' length 1 or same as `y`
##' @param uptake_first_dose numeric (0-1) either length 1 or same as `y`
##' @param uptake_second_dose numeric (0-1) either length 1 or same as `y`
##' @param disc_rate discount rate for cost-effectiveness output per annum,
##'  default = 0
##' @return A list of matrices with rows = time point and columns
##' = parameter set
##' list entries are:
##' `cum_treated` = cumulative number treated in model run
##' `cum_vaccinated` = cumulative vaccinated treated in model run
##' `diag_a` = annual number of asymptomatic diagnoses in model run
##' `diag_s` = annual number of symptomatic diagnoses in model run
##' `treated` = annual number treated in model run
##' `screened` = annual number screened in model run
##' `vaccinated` = annual number vaccinated in model run
##' (includes those getting vbe and boosters)
##' `vbe` = annual number vaccinated before entry in model run
##' `revaccinated` = annual number receiving booster in model run
##' `offered_primary` = annual number offered primary vaccination
##' (does not include hesitant or vbe)
##' `inc_cum_treated` = cumulative number treated compared to baseline
##' `inc_cum_vaccinated` = cumulative number vaccinated compared to baseline
##' `inc_diag_a` = annual number of asymptomatic diagnoses compared to baseline
##' `inc_diag_s` = annual number of symptomatic diagnoses compared to baseline
##' `inc_treated` = annual number treated compared to baseline
##' `inc_screened` = annual number screened compared to baseline
##' `inc_vaccinated`  = annual number vaccinated compared to baseline
##' (includes those getting vbe and boosters)
##' `inc_vbe` = annual number vaccinated before entry compared to baseline
##' `inc_revaccinated` = annual number receiving booster compared to
##' baseline
##' `inc_cum_revaccinated` = cumulative number revaccinated by recieving booster
##' compared to baseline
##' `inc_offered_primary` = annual number offered primary vaccination compared
##' to baseline (does not include hesitant or vbe)
##' `inc_primary` = annual number receiving primary vaccination
##' compared to baseline (does not include vbe)
##' `inc_cum_primary` = cumulative number of individuals receiving primary
##' vaccination compared to baseline (does not include vbe)
##' `inc_doses` = number of doses compared to baseline (per year). Assumes
##' primary vaccination uses 2 doses, booster uses 1 dose.
##' `inc_cum_doses` = cumulative number of doses compared to baseline. Assumes
##' primary vaccination uses 2 doses, booster uses 1 dose.
##' `inc_primary_doses` = annual number of primary doses administered
##' compared to baseline, where 2 doses are needed for full protection.
##' `inc_cum_primary_doses` = cumulative number of primary doses administered
##' compared to baseline, where 2 doses are needed for full protection.
##' `inc_booster_doses` = annual number of booster doses
##' administered to waned individuals for revaccination compared to baseline,
##' where 1 dose is needed for full protection.
##' `inc_cum_booster_doses` = cumulative number of booster doses
##' administered to waned individuals for revaccination compared to baseline,
##'  where 1 dose is needed for full protection.
##'  `inc_vbe_doses` = annual number of vaccinations given before entry,
##'  where 2 doses give full protection.
##'  `inc_cum_vbe_doses` = cumulative number of vaccinations given before entry,
##'  where 2 doses give full protection.
##' `cases_averted_per_dose` = cumulative number of cases (i.e. diagnoses)
##' averted per dose of vaccine
##' `cases_averted_per_dose_pv` = present value of cases_averted_per_dose
##' (i.e. sum of annual numbers discounted at rate `disc_rate` to time 0)
##' `inc_doses_wasted` annual number of first doses given that are not followed
##' by a second dose, compared to baseline
##' `inc_cum_doses_wasted` cumulative number of first doses given that are not
##'  followed by a second dose, compared to baseline
##' `pv_inc_doses` = present value of number of doses compared to baseline
##' (i.e. sum of annual numbers discounted at rate `disc_rate` to time 0)
##' `pv_red_net_cost` = present value of the reduction in net costs compared to
##' baseline (i.e. sum of annual reduction in costs discounted at rate
##'  `disc_rate` to time 0)
##' `pv_qaly_gain` = present value of the QALYs gained due to vaccination
##' compared to baseline (i.e. QALYs gained each year discounted at rate
##'  `disc_rate` to time 0, and summed)
##' `cet_20k` = cost effectiveness threshold (i.e. £ value per dose at which
##' vaccination becomes cost-effective) calculated using £20,000 / QALY.
##' (pv_red_net_cost + pv_qaly_gain * £20,000) / pv_inc_doses
##' `cet_30k` = cost effectiveness threshold (i.e. £ value per dose at which
##' vaccination becomes cost-effective) calculated using £30,000 / QALY
##' (pv_red_net_cost + pv_qaly_gain * £30,000) / pv_inc_doses
##' `inc_costs_9` = present value of incremental costs assuming £9 / dose.
##' Incremental costs are calculated as: pv_inc_doses * £9 - pv_red_net_cost
##' `inc_costs_18` = present value of incremental costs assuming £18 / dose.
##' Incremental costs are calculated as: pv_inc_doses * £18 - pv_red_net_cost
##' `inc_costs_50` = present value of incremental costs assuming £50 / dose.
##' Incremental costs are calculated as: pv_inc_doses * £50 - pv_red_net_cost
##' `inc_costs_85` = present value of incremental costs assuming £85 / dose
##' Incremental costs are calculated as: pv_inc_doses * £85 - pv_red_net_cost
##' @export
compare_baseline <- function(y, baseline, uptake_first_dose,
                             uptake_second_dose, cost_params,
                             disc_rate) {

  ## compare run to baseline
  flows <- extract_flows(y)
  ret <- Map(`-`, flows, baseline[names(flows)])
  names(ret) <- paste0("inc_", names(flows))
  ret <- c(flows, ret)

  ## calculate number receiving primary vaccination
  ret$inc_primary <- ret$inc_vaccinated - ret$inc_revaccinated - ret$inc_vbe
  ret$inc_cum_primary <- apply(ret$inc_primary, 2, cumsum)

  ## cumulative individuals receiving booster vaccination (revaccination)
  ret$inc_cum_revaccinated <- apply(ret$inc_revaccinated, 2, cumsum)

  ## calculate cases averted per dose, both with and without discounting
  ## return incremental annual and cumulative doses
  # all vbe get two doses
  ret$inc_vbe_doses <- ret$inc_vbe * 2
  # calculate doses given per person offered primary vaccination
  ret$inc_primary_doses <- uptake_first_dose * (1 + uptake_second_dose) *
    ret$inc_offered_primary
  # booster vaccination takes a single dose so is simply the number of people
  # revaccinated
  ret$inc_booster_doses <- ret$inc_revaccinated

  # calculate cumulative doses
  ret$inc_cum_vbe_doses <- apply(ret$inc_vbe_doses, 2, cumsum)
  ret$inc_cum_primary_doses <- apply(ret$inc_primary_doses, 2, cumsum)
  ret$inc_cum_booster_doses <- apply(ret$inc_booster_doses, 2, cumsum)

  ret$inc_doses <- ret$inc_primary_doses + ret$inc_booster_doses +
    ret$inc_vbe_doses
  ret$inc_cum_doses <- apply(ret$inc_doses, 2, cumsum)

  ret$cases_averted_per_dose <- calc_cases_averted_per_dose(ret, 0)
  ret$cases_averted_per_dose_pv <- calc_cases_averted_per_dose(ret, disc_rate)

  ## calculate vaccine doses wasted
  ret$inc_doses_wasted <-
    ret$inc_offered_primary * uptake_first_dose * (1 - uptake_second_dose)
  ret$inc_cum_doses_wasted <- apply(ret$inc_doses_wasted, 2, cumsum)

  ## calculate costs
  costs <- calc_costs(ret, cost_params, disc_rate)
  ret <- c(ret, costs)

  ## calculate cost-effectiveness threshold price as incremental net benefit per
  ## person with QALY cost - £20k / £30k
  ## both calcs allow for discounting (i.e. are present values as at 2022)
  ret$cet_20k <- calc_cet(2e4, costs)
  ret$cet_30k <- calc_cet(3e4, costs)

  ## calculate incremental cost of vaccination assuming
  ## £9, £18, £50, £85 per dose
  ## both calcs allow for discounting (i.e. are present values as at 2022)
  ret$inc_costs_9 <- calc_inc_costs(9, costs)
  ret$inc_costs_18 <- calc_inc_costs(18, costs)
  ret$inc_costs_50 <- calc_inc_costs(50, costs)
  ret$inc_costs_85 <- calc_inc_costs(85, costs)

  ret
}


calc_pv <- function(x, disc_rate) {

  tt <- seq_len(nrow(x))
  pv <- (1 + disc_rate) ^ - (tt - 0.5)
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


calc_costs <- function(forecast, cost_params, disc_rate) {

  red_diag_s <- -t(forecast$inc_diag_s)
  red_diag_a <- -t(forecast$inc_diag_a)
  inc_screen <- t(forecast$inc_screened)

  red_cost_diag_s <- red_diag_s * cost_params$unit_cost_manage_symptomatic
  red_cost_diag_a <- red_diag_a * cost_params$unit_cost_manage_asymptomatic
  inc_cost_screen <- inc_screen * cost_params$unit_cost_screen_uninfected

  red_net_cost <- t(red_cost_diag_s + red_cost_diag_a - inc_cost_screen)

  qaly_gain <- t(red_diag_s * cost_params$qaly_loss_per_diag_s)

  ## calc pv incremental net benefit per person protected
  list(pv_inc_doses = calc_pv(forecast$inc_doses, disc_rate),
       pv_red_net_cost = calc_pv(red_net_cost, disc_rate),
       pv_qaly_gain = calc_pv(qaly_gain, disc_rate))
}

calc_cet <- function(cost_qaly, costs) {
  # PV inc net benefit = PV reduction in treatment costs + PV gain in QALYs
  pv_inc_net_benefit <- costs$pv_red_net_cost + costs$pv_qaly_gain * cost_qaly
  # value of vaccination = PV inc net benefit / PV inc doses
  pv_inc_net_benefit / costs$pv_inc_doses
}

calc_inc_costs <- function(price_per_dose, costs) {
  costs$pv_inc_doses * price_per_dose - costs$pv_red_net_cost
}
