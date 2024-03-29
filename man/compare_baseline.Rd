% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_grid.R
\name{compare_baseline}
\alias{compare_baseline}
\title{compare model runs with vaccination to a baseline runs}
\usage{
compare_baseline(
  y,
  baseline,
  uptake_first_dose,
  uptake_second_dose,
  cost_params,
  disc_rate
)
}
\arguments{
\item{y}{list of model runs, e.g. created by `run_onevax_xvwv`, each list
entry refers to a different parameter set}

\item{baseline}{list of baseline runs e.g. created by `run_onevax_xvwv`.
Should be same length as `y`}

\item{uptake_first_dose}{numeric (0-1) either length 1 or same as `y`}

\item{uptake_second_dose}{numeric (0-1) either length 1 or same as `y`}

\item{cost_params}{list of cost effectiveness parameters, containing entries
`qaly_loss_per_diag_s`, `unit_cost_manage_symptomatic`,
`unit_cost_manage_asymptomatic`, `unit_cost_screen_uninfected`. Each entry
length 1 or same as `y`}

\item{disc_rate}{discount rate for cost-effectiveness output per annum,
default = 0}
}
\value{
A list of matrices with rows = time point and columns
= parameter set
list entries are:
`cum_treated` = cumulative number treated in model run
`cum_vaccinated` = cumulative vaccinated treated in model run
`diag_a` = annual number of asymptomatic diagnoses in model run
`diag_s` = annual number of symptomatic diagnoses in model run
`treated` = annual number treated in model run
`screened` = annual number screened in model run
`vaccinated` = annual number vaccinated in model run
(includes those getting vbe and boosters)
`vbe` = annual number vaccinated before entry in model run
`revaccinated` = annual number receiving booster in model run
`offered_primary` = annual number offered primary vaccination
(does not include hesitant or vbe)
`inc_cum_treated` = cumulative number treated compared to baseline
`inc_cum_vaccinated` = cumulative number vaccinated compared to baseline
`inc_diag_a` = annual number of asymptomatic diagnoses compared to baseline
`inc_diag_s` = annual number of symptomatic diagnoses compared to baseline
`inc_treated` = annual number treated compared to baseline
`inc_screened` = annual number screened compared to baseline
`inc_vaccinated`  = annual number vaccinated compared to baseline
(includes those getting vbe and boosters)
`inc_vbe` = annual number vaccinated before entry compared to baseline
`inc_revaccinated` = annual number receiving booster compared to
baseline
`inc_cum_revaccinated` = cumulative number revaccinated by recieving booster
compared to baseline
`inc_offered_primary` = annual number offered primary vaccination compared
to baseline (does not include hesitant or vbe)
`inc_primary` = annual number receiving primary vaccination
compared to baseline (does not include vbe)
`inc_cum_primary` = cumulative number of individuals receiving primary
vaccination compared to baseline (does not include vbe)
`inc_doses` = number of doses compared to baseline (per year). Assumes
primary vaccination uses 2 doses, booster uses 1 dose.
`inc_cum_doses` = cumulative number of doses compared to baseline. Assumes
primary vaccination uses 2 doses, booster uses 1 dose.
`inc_primary_doses` = annual number of primary doses administered
compared to baseline, where 2 doses are needed for full protection.
`inc_cum_primary_doses` = cumulative number of primary doses administered
compared to baseline, where 2 doses are needed for full protection.
`inc_booster_doses` = annual number of booster doses
administered to waned individuals for revaccination compared to baseline,
where 1 dose is needed for full protection.
`inc_cum_booster_doses` = cumulative number of booster doses
administered to waned individuals for revaccination compared to baseline,
 where 1 dose is needed for full protection.
 `inc_vbe_doses` = annual number of vaccinations given before entry,
 where 2 doses give full protection.
 `inc_cum_vbe_doses` = cumulative number of vaccinations given before entry,
 where 2 doses give full protection.
`cases_averted_per_dose` = cumulative number of cases (i.e. diagnoses)
averted per dose of vaccine
`cases_averted_per_dose_pv` = present value of cases_averted_per_dose
(i.e. sum of annual numbers discounted at rate `disc_rate` to time 0)
`inc_doses_wasted` annual number of first doses given that are not followed
by a second dose, compared to baseline
`inc_cum_doses_wasted` cumulative number of first doses given that are not
 followed by a second dose, compared to baseline
`pv_inc_doses` = present value of number of doses compared to baseline
(i.e. sum of annual numbers discounted at rate `disc_rate` to time 0)
`pv_red_net_cost` = present value of the reduction in net costs compared to
baseline (i.e. sum of annual reduction in costs discounted at rate
 `disc_rate` to time 0)
`pv_qaly_gain` = present value of the QALYs gained due to vaccination
compared to baseline (i.e. QALYs gained each year discounted at rate
 `disc_rate` to time 0, and summed)
`cet_20k` = cost effectiveness threshold (i.e. £ value per dose at which
vaccination becomes cost-effective) calculated using £20,000 / QALY.
(pv_red_net_cost + pv_qaly_gain * £20,000) / pv_inc_doses
`cet_30k` = cost effectiveness threshold (i.e. £ value per dose at which
vaccination becomes cost-effective) calculated using £30,000 / QALY
(pv_red_net_cost + pv_qaly_gain * £30,000) / pv_inc_doses
`inc_costs_9` = present value of incremental costs assuming £9 / dose.
Incremental costs are calculated as: pv_inc_doses * £9 - pv_red_net_cost
`inc_costs_18` = present value of incremental costs assuming £18 / dose.
Incremental costs are calculated as: pv_inc_doses * £18 - pv_red_net_cost
`inc_costs_50` = present value of incremental costs assuming £50 / dose.
Incremental costs are calculated as: pv_inc_doses * £50 - pv_red_net_cost
`inc_costs_85` = present value of incremental costs assuming £85 / dose
Incremental costs are calculated as: pv_inc_doses * £85 - pv_red_net_cost
}
\description{
compare model runs with vaccination to a baseline runs
}
