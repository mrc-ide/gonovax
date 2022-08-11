test_that("calc_pv works as expected", {

  ## time is dimension 1
  x <- matrix(seq_len(6), 3, 2)
  y0 <- calc_pv(x, 0)

  expect_equal(y0, apply(x, 2, cumsum))

  d <- 0.1
  disc_fac <- (1 + d) ^ - (1:3 - 0.5)

  y1 <- calc_pv(x, d)
  expect_equal(y1[1, ], x[1, ] * disc_fac[1])
  expect_equal(y1[3, ], c(disc_fac %*% x))
})

test_that("compare baseline xpvwrh works as expected", {

  gp <- gono_params(1:2)
  ip <- lapply(run_onevax_xpvwrh(0:1, gp), restart_params)
  tt <- 1:4

  # baseline of 50% uptake vbe with a 50% eff vaccine lasting 1 year
  bl <- extract_flows_xpvwrh(run_onevax_xpvwrh(tt, gp, ip, vea = 0.5, dur_v = 1,
                                      vbe = 0.5))
  blv <- rep(list(bl), 4)
  cp <- list(qaly_loss_per_diag_s = c(0.002, 0.001),
             unit_cost_manage_symptomatic = c(98, 99),
             unit_cost_manage_asymptomatic = c(93, 94),
             unit_cost_screen_uninfected = c(70, 71))

  r2 <- 0.7
  r1 <- 0.9

  # compare to adding VoD with same vaccine, uptake = 63%
  y <- run_onevax_xpvwrh(tt, gp, ip, vea = 0.5, dur_v = 1, vbe = 0.5,
                       r1 = r1, r2 = r2, strategy = "VoD")
  yy <- extract_flows_xpvwrh(y)

  z <- compare_baseline_xpvwrh(y, bl, r1, r2, cp, 0)

  flownames <- names(yy)
  incnames <- paste0("inc_", flownames)
  compnames <- setdiff(names(z), c(flownames, incnames))

  # check values are carried over correctly
  expect_equal(z[flownames], yy)
  expect_equivalent(z[incnames],
                    Map(`-`, yy[flownames], bl[flownames]))

  f <- function(x) {
    pv <- x$vaccinated - x$revaccinated - x$vbe - x$part_to_full
    2 * (pv * r2) + (pv * (1 - r2)) + 1 * x$revaccinated +
                                            1 * x$part_to_full +
                                            2 * x$vbe
  }

  # check doses are calculated correctly
  # no need to worry about vbe here as it is the same in baseline and run
  n <- length(tt) - 1
  expect_equal(f(yy) - f(bl), z$inc_doses)
  expect_equal(z$inc_cum_doses[1, ], z$inc_doses[1, ])
  expect_equal(z$inc_cum_doses[n, ], colSums(z$inc_doses))
  expect_equal(sum(z$inc_vbe), 0)

  expect_equal(z$cases_averted_per_dose[1, ],
               -z$inc_treated[1, ] / z$inc_doses[1, ])
  expect_equal(z$cases_averted_per_dose[n, ],
               -colSums(z$inc_treated) / colSums(z$inc_doses))

  ## check against a baseline of no vaccination
  bl0 <- extract_flows_xpvwrh(run_onevax_xpvwrh(tt, gp, ip))
  blv0 <- rep(list(bl0), 4)

  z1 <- compare_baseline_xpvwrh(y, bl0, r1, r2, cp, 0)
  expect_equivalent(z1$inc_vbe, matrix(12000 * 0.5, 3, 2))
  expect_equal(f(yy) - f(bl0), z1$inc_doses)

  ## with zero disc rate
  expect_equal(z$cases_averted_per_dose, z$cases_averted_per_dose_pv)

  tmp <- calc_cases_averted_per_dose(z, 0.03)
  expect_equal(z$cases_averted_per_dose[1, ], tmp[1, ])
  expect_true(all(tmp[2, ] != z$cases_averted_per_dose[1, ]))

  ## when r2 = 1, r2_p = 0, booster_uptake = 0, total doses = 2 * # vaccinations
  y <- run_onevax_xpvwrh(tt, gp, ip,
                         r1 = 1, r2 = 1, booster_uptake = 0,
                         r2_p = 0, strategy = "VoD")
  z2 <- compare_baseline_xpvwrh(y, bl0, uptake_first_dose = 1,
                                 uptake_second_dose = 1, cp, 0)
  expect_equal(z2$inc_cum_doses, (2 * z2$inc_cum_vaccinated))

  ## when r2 = 0, vbe = 0, total doses = total vaccinations
  y <- run_onevax_xpvwrh(tt, gp, ip,
                         r1 = 1, r2 = 0, booster_uptake = 0,
                         r2_p = 0, vbe = 0, strategy = "VoD")
  z2 <- compare_baseline_xpvwrh(y, bl0, uptake_first_dose = 1,
                                uptake_second_dose = 0, cp, 0)
  expect_equal(z2$inc_cum_doses, z2$inc_cum_vaccinated)

  ## test cost eff threshold
  d <- 0
  cost <- calc_costs(z, cp, d)
  cet_20k <- calc_cet(2e4, cost)
  cet_30k <- calc_cet(3e4, cost)

  expect_equivalent(cet_20k,
                    matrix(c(0.762513348510373, 2.45450640004958,
                             4.99544894499636, 0.690487829848825,
                             2.29553458279502, 4.74762990475995),
                             nrow = 3))
  expect_equal(cet_30k,
               matrix(c(0.871045465836508, 2.80048568207085, 5.69605378633193,
                        0.748390761083548, 2.4868096465653, 5.14203978621594),
                        nrow = 3))

  expect_equal(z$cet_20k, cet_20k)
  expect_equal(z$cet_30k, cet_30k)

  calc_cost <- function(x, cp, qc) {
    Qu <- cp$unit_cost_screen_uninfected * t(x$screened)
    Qs <- (cp$unit_cost_manage_symptomatic + qc * cp$qaly_loss_per_diag_s) *
      t(x$diag_s)
    Qa <- cp$unit_cost_manage_asymptomatic * t(x$diag_a)

    t(Qu + Qs + Qa)
  }

  expect_equal(calc_pv(calc_cost(bl, cp, 2e4) - calc_cost(yy, cp, 2e4), 0) /
                 calc_pv(z$inc_doses, 0),
               z$cet_20k)

  expect_equal(calc_pv(calc_cost(bl, cp, 3e4) - calc_cost(yy, cp, 3e4), 0) /
                 calc_pv(z$inc_doses, 0),
               z$cet_30k)



  expect_equal(cet_20k * cost$pv_inc_doses - cost$pv_red_net_cost,
               cost$pv_qaly_gain * 20000)
  expect_equal(cet_30k * cost$pv_inc_doses - cost$pv_red_net_cost,
               cost$pv_qaly_gain * 30000)

  expect_equal(calc_inc_costs(85, cost) - calc_inc_costs(9, cost),
               cost$pv_inc_doses * (85 - 9))
  expect_equal(calc_inc_costs(85, cost) - calc_inc_costs(18, cost),
               cost$pv_inc_doses * (85 - 18))
  expect_equal(calc_inc_costs(85, cost) - calc_inc_costs(35, cost),
               cost$pv_inc_doses * (85 - 35))
  expect_equal(calc_inc_costs(85, cost) - calc_inc_costs(50, cost),
               cost$pv_inc_doses * (85 - 50))
  expect_equal(calc_inc_costs(85, cost) - calc_inc_costs(70, cost),
               cost$pv_inc_doses * (85 - 70))

  ## test that proportion of the population vaccine protected is calculated
  # correctly

   # when uptake = 0, vaccine protected raw and proportion = 0
   y <- run_onevax_xpvwrh(tt, gp,
                         r1 = 0,
                         r2 = 0, booster_uptake = 0,
                         strategy = "VoD")
   bl <- extract_flows_xpvwrh(run_onevax_xpvwrh(tt, gp, ip))

   z <- compare_baseline_xpvwrh(y, bl, uptake_first_dose = 0,
                               uptake_second_dose = 0, cp, 0)

    expect_true(all(z$inc_vacprotec_total == 0))
    expect_true(all(z$inc_vacprotec_total_prop == 0))

   # when r2 = 1, vac_protect_part = 0

   y <- run_onevax_xpvwrh(tt, gp,
                           r1 = 0.5,
                           r2 = 1, booster_uptake = 0,
                           strategy = "VoD")

   z <- compare_baseline_xpvwrh(y, bl, uptake_first_dose = 0.5,
                                 uptake_second_dose = 0.2, cp, 0)

  expect_true(all(z$inc_vacprotec_part == 0))
  expect_true(all(z$inc_vacprotec_part_prop == 0))
  expect_true(all(z$inc_vacprotec_total != 0))
  expect_true(all(z$inc_vacprotec_total_prop != 0))
  expect_equal(z$inc_vacprotec_total, z$inc_vacprotec_full)

  # proportion * 6e05 = raw number
  expect_equal(z$inc_vacprotec_total_prop * 6e05, z$inc_vacprotec_total)

  ## test that cumulative primary and booster vaccination calc correctly

  bl <- extract_flows_xpvwrh(run_onevax_xpvwrh(tt, gp))

  y <- run_onevax_xpvwrh(tt, gp, vea = 1, dur_v = 1,
                         r1 = 1,
                         r2 = 1, booster_uptake = 0,
                          strategy = "VoD")
  z <- compare_baseline_xpvwrh(y, bl, uptake_first_dose = 1,
                        uptake_second_dose = 1, cp, 0)

  # number undergoing primary vaccination is equal to total people
  # vaccinated when booster uptake = 0

  expect_equal(z$inc_primary[length(tt) - 1, 1],
               z$inc_vaccinated[length(tt) - 1, 1])

  # number undergoing primary vaccination over 3 years is equal to the
  # cumulative number for the 3 years

  expect_equal(sum(z$inc_primary),
               sum(z$inc_cum_primary[length(tt) - 1, ]))

   # number booster vaccinations is equal to
      # inc_vaccinated - inc_primary
      # sum re-vaccinated

  y <- run_onevax_xpvwrh(tt, gp, vea = 1, dur_v = 1,
                         r1 = 1,
                         r2 = 1, booster_uptake = 0.5,
                        strategy = "VoD")
  z <- compare_baseline_xpvwrh(y, bl, uptake_first_dose = 1,
                                      uptake_second_dose = 1,
                                      cp, 0)

  expect_equal(z$inc_vaccinated - z$inc_primary - z$inc_vbe, z$inc_revaccinated)

  # number undergoing booster vaccination over 3 years is equal to the
  # cumulative number for the 3 years

  expect_equal(sum(z$inc_revaccinated),
               sum(z$inc_cum_revaccinated[length(tt) - 1, ]))

  # cumulative primary vaccination + cumulative booster vaccination =
  # cumulative vaccinated overall when vbe = 0

expect_equal(z$inc_cum_primary + z$inc_cum_part_to_full +
               z$inc_cum_revaccinated, z$inc_cum_vaccinated)

  # cumulative doses of different types calculated correctly
  # sum of vbe, primary and revaccination doses = all doses

  y <- run_onevax_xvwrh(tt, gp, vea = 1, dur = 1,
                        primary_uptake = 1, booster_uptake = 0.5,
                        strategy = "VoD")
  z <- compare_baseline_xpvwrh(y, bl, uptake_first_dose = 1,
                                      uptake_second_dose = 1,
                                      cp, 0)

  s <- sum(z$inc_cum_primary_total_doses[length(tt) - 1, ] +
             z$inc_cum_part_to_full_doses[length(tt) - 1, ] +
             z$inc_cum_booster_doses[length(tt) - 1, ] +
             z$inc_cum_vbe_doses[length(tt) - 1, ])

  t <- sum(z$inc_cum_doses[length(tt) - 1, ])

  expect_equal(s, t)

})
###############################################################################
test_that("compare baseline works as expected", {

  gp <- gono_params(1:2)
  ip <- lapply(run_onevax_xvwv(0:1, gp), restart_params)
  tt <- 1:4

  # baseline of 50% uptake vbe with a 50% eff vaccine lasting 1 year
  bl <- extract_flows(run_onevax_xvwv(tt, gp, ip, vea = 0.5, dur = 1,
                                      vbe = 0.5))
  blv <- rep(list(bl), 4)
  cp <- list(qaly_loss_per_diag_s = c(0.002, 0.001),
             unit_cost_manage_symptomatic = c(98, 99),
             unit_cost_manage_asymptomatic = c(93, 94),
             unit_cost_screen_uninfected = c(70, 71))

  r2 <- 0.7
  r1 <- 0.9

  # compare to adding VoD with same vaccine, uptake = 63%
  y <- run_onevax_xvwv(tt, gp, ip, vea = 0.5, dur = 1, vbe = 0.5,
                       uptake = r1 * r2, strategy = "VoD")
  yy <- extract_flows(y)

  z <- compare_baseline(y, bl, r1, r2, cp, 0)

  flownames <- names(yy)
  incnames <- paste0("inc_", flownames)
  compnames <- setdiff(names(z), c(flownames, incnames))

  # check values are carried over correctly
  expect_equal(z[flownames], yy)
  expect_equivalent(z[incnames],
                    Map(`-`, yy[flownames], bl[flownames]))

  f <- function(x) {
    pv <- x$vaccinated - x$revaccinated - x$vbe
    pv * (1 + 1 / r2) + 1 * x$revaccinated + 2 * x$vbe
  }

  # check doses are calculated correctly
  # no need to worry about vbe here as it is the same in baseline and run
  n <- length(tt) - 1
  expect_equal(f(yy) - f(bl), z$inc_doses)
  expect_equal(z$inc_cum_doses[1, ], z$inc_doses[1, ])
  expect_equal(z$inc_cum_doses[n, ], colSums(z$inc_doses))
  expect_equal(sum(z$inc_vbe), 0)

  expect_equal(z$cases_averted_per_dose[1, ],
               -z$inc_treated[1, ] / z$inc_doses[1, ])
  expect_equal(z$cases_averted_per_dose[n, ],
               -colSums(z$inc_treated) / colSums(z$inc_doses))

  ## note that inc_offered_primary * r1 * (1- r2) = inc_doses_wasted
  ## and inc_primary = inc_offered_primary * r1 * r2
  expect_equal(z$inc_primary / r2, z$inc_doses_wasted / (1 - r2))

  ## check against a baseline of no vaccination
  bl0 <- extract_flows(run_onevax_xvwv(tt, gp, ip))
  blv0 <- rep(list(bl0), 4)

  z1 <- compare_baseline(y, bl0, r1, r2, cp, 0)
  expect_equivalent(z1$inc_vbe, matrix(12000 * 0.5, 3, 2))
  expect_equal(f(yy) - f(bl0), z1$inc_doses)

  ## with zero disc rate
  expect_equal(z$cases_averted_per_dose, z$cases_averted_per_dose_pv)

  tmp <- calc_cases_averted_per_dose(z, 0.03)
  expect_equal(z$cases_averted_per_dose[1, ], tmp[1, ])
  expect_true(all(tmp[2, ] != z$cases_averted_per_dose[1, ]))

  ## test cost eff threshold
  d <- 0
  cost <- calc_costs(z, cp, d)
  cet_20k <- calc_cet(2e4, cost)
  cet_30k <- calc_cet(3e4, cost)

  expect_equivalent(cet_20k,
                    matrix(c(0.534243384207646, 1.72471712423966,
                             3.52594421792319, 0.483751150425284,
                             1.61256157784351, 3.34884963770966), nrow = 3))
  expect_equal(cet_30k,
               matrix(c(0.610287112367831, 1.96785021889529,
                        4.02052903497914, 0.524318410369186,
                        1.74693578901651, 3.62708130277735), nrow = 3))

  expect_equal(z$cet_20k, cet_20k)
  expect_equal(z$cet_30k, cet_30k)

  calc_cost <- function(x, cp, qc) {
    Qu <- cp$unit_cost_screen_uninfected * t(x$screened)
    Qs <- (cp$unit_cost_manage_symptomatic + qc * cp$qaly_loss_per_diag_s) *
      t(x$diag_s)
    Qa <- cp$unit_cost_manage_asymptomatic * t(x$diag_a)

    t(Qu + Qs + Qa)
  }

  expect_equal(calc_pv(calc_cost(bl, cp, 2e4) - calc_cost(yy, cp, 2e4), 0) /
                 calc_pv(z$inc_doses, 0),
               z$cet_20k)

  expect_equal(calc_pv(calc_cost(bl, cp, 3e4) - calc_cost(yy, cp, 3e4), 0) /
                 calc_pv(z$inc_doses, 0),
               z$cet_30k)

  expect_equal(cet_20k * cost$pv_inc_doses - cost$pv_red_net_cost,
               cost$pv_qaly_gain * 20000)
  expect_equal(cet_30k * cost$pv_inc_doses - cost$pv_red_net_cost,
               cost$pv_qaly_gain * 30000)

  expect_equal(calc_inc_costs(85, cost) - calc_inc_costs(9, cost),
               cost$pv_inc_doses * (85 - 9))
  expect_equal(calc_inc_costs(85, cost) - calc_inc_costs(18, cost),
               cost$pv_inc_doses * (85 - 18))
  expect_equal(calc_inc_costs(85, cost) - calc_inc_costs(35, cost),
               cost$pv_inc_doses * (85 - 35))
  expect_equal(calc_inc_costs(85, cost) - calc_inc_costs(50, cost),
               cost$pv_inc_doses * (85 - 50))
  expect_equal(calc_inc_costs(85, cost) - calc_inc_costs(70, cost),
               cost$pv_inc_doses * (85 - 70))

  ## test doses wasted calculated correctly

  # no doses wasted when second dose uptake = 100%
  z <- compare_baseline(y, bl, r1 * r2, uptake_second_dose = 1, cp, 0)

  expect_true(all(z$inc_cum_doses_wasted == 0))
  expect_true(all(z$inc_doses_wasted == 0))

  # all doses wasted when second dose uptake = 0%
  y1 <- run_onevax_xvwv(tt, gp, ip, vea = 0.5, dur = 1, vbe = 0.5,
                        uptake = 0, strategy = "VoD")
  z <- compare_baseline(y1, bl, r1, uptake_second_dose = 0, cp, 0)
  expect_equal(z$inc_cum_doses_wasted, z$inc_cum_doses)

  # regardless of hesitancy (primary uptake = booster uptake = 1)

  bl_h <- extract_flows(run_onevax_xvwrh(tt, gp, hes = 0.3))
  y_h <- run_onevax_xvwrh(tt, gp, vea = 0.5, dur = 1, vbe = 0.5,
                          primary_uptake = r1, strategy = "VoD", hes = 0.3)

  z_h <- compare_baseline(y_h, bl_h,
                          uptake_first_dose = r1, uptake_second_dose = 1,
                          cp, 0)

  expect_true(all(z_h$inc_cum_doses_wasted == 0))
  expect_true(all(z_h$inc_doses_wasted == 0))

  # all primary doses wasted when second dose uptake = 0

  y_2 <- run_onevax_xvwrh(tt, gp, vea = 0.5, dur = 1,
                          primary_uptake = 0, strategy = "VoD", hes = 0.3)
  z_f2 <- compare_baseline(y_2, bl_h,
                           uptake_first_dose = r1, uptake_second_dose = 0,
                           cp, 0)
  expect_equal(z_f2$inc_cum_doses_wasted, z_f2$inc_cum_doses)

  # doses used + doses wasted = 2 * #primary doses + 1 * booster doses
  # + 2 * vbe doses

  y_3 <- run_onevax_xvwrh(tt, gp, vea = 0.5, dur = 1, vbe = 0.5,
                          primary_uptake = 0, strategy = "VoD", hes = 0.3)
  expect_equal(z_h$inc_doses + z_h$inc_doses_wasted,
               2 * z_h$inc_primary + 1 * z_h$inc_revaccinated +
                 2 * z_h$inc_vbe)

  ## test that cumulative primary and booster vaccination calc correctly

  bl <- extract_flows(run_onevax_xvwrh(tt, gp))

  y <- run_onevax_xvwrh(tt, gp, vea = 1, dur = 1,
                        primary_uptake = 1, booster_uptake = 0,
                        strategy = "VoD")
  z <- compare_baseline(y, bl, uptake_first_dose = 1,
                        uptake_second_dose = 1, cp, 0)

  # number undergoing primary vaccination is equal to total people
  # vaccinated when booster uptake = 0

  expect_equal(z$inc_primary[length(tt) - 1, 1],
               z$inc_vaccinated[length(tt) - 1, 1])

  # number undergoing primary vaccination over 3 years is equal to the
  # cumulative number for the 3 years
  expect_equal(sum(z$inc_primary),
               sum(z$inc_cum_primary[length(tt) - 1, ]))

  # number booster vaccinations is equal to
  # inc_vaccinated - inc_primary
  # sum re-vaccinated

  y <- run_onevax_xvwrh(tt, gp, vea = 1, dur = 1,
                        primary_uptake = 1, booster_uptake = 0.5,
                        strategy = "VoD")
  z <- compare_baseline(y, bl, uptake_first_dose = 1, uptake_second_dose = 1,
                        cp, 0)

  expect_equal(z$inc_vaccinated - z$inc_primary - z$inc_vbe, z$inc_revaccinated)

  # number undergoing booster vaccination over 3 years is equal to the
  # cumulative number for the 3 years

  expect_equal(sum(z$inc_revaccinated),
               sum(z$inc_cum_revaccinated[length(tt) - 1, ]))

  # cumulative primary vaccination + cumulative booster vaccination =
  # cumulative vaccinated overall when vbe = 0

  expect_equal(z$inc_cum_primary + z$inc_cum_revaccinated, z$inc_cum_vaccinated)

  # cumulative doses of different types calculated correctly
  # sum of vbe, primary and revaccination doses = all doses

  y <- run_onevax_xvwrh(tt, gp, vea = 1, dur = 1,
                        primary_uptake = 1, booster_uptake = 0.5,
                        strategy = "VoD")
  z <- compare_baseline(y, bl, uptake_first_dose = 1, uptake_second_dose = 1,
                        cp, 0)

  s <- sum(z$inc_cum_primary_doses[length(tt) - 1, ] +
             z$inc_cum_booster_doses[length(tt) - 1, ] +
             z$inc_cum_vbe_doses[length(tt) - 1, ])

  t <- sum(z$inc_cum_doses[length(tt) - 1, ])

  expect_equal(s, t)

})


test_that("run_grid works as expected", {
  ## test with vbe only

  gp <- gono_params(1:2)
  ip <- lapply(run_onevax_xvwv(0:1, gp, vea = 0, dur = 1), restart_params)
  tt <- 1:3

  bl <- extract_flows(run_onevax_xvwv(tt, gp, ip, vea = 0, dur = 1))
  blv <- rep(list(bl), 4)
  cp <- list(qaly_loss_per_diag_s = c(0.002, 0.001),
             unit_cost_manage_symptomatic = c(98, 99),
             unit_cost_manage_asymptomatic = c(93, 94),
             unit_cost_screen_uninfected = c(70, 71))

  y <- run_onevax_xvwv(tt, gp, ip, vea = 0, dur = 1, vbe = 0.5)
  z <- compare_baseline(y, bl, 0.8, 0.7, cp, 0)

  zz <- run_grid(gp, ip, cp, blv,
                model = run_onevax_xvwv,
                eff = c(0, 1), dur = c(1, 2), vbe = 0.5,
                uptake_total = 1,
                full_output = TRUE)
  expect_equal(zz$results$eff0.00_dur01$inc_treated, matrix(0, 2, 2),
               tolerance = 0.1)
  expect_equal(zz$results$eff0.00_dur02$inc_treated,
               zz$results$eff0.00_dur01$inc_treated,
               tolerance = 0.1)
  expect_true(all(unlist(zz$results$inc_treated$eff0.10_dur01) > -0.1))
  expect_equal(zz$results$eff0.00_dur01$inc_cum_treated,  matrix(0, 2, 2),
               tolerance = 0.1)
  expect_equal(zz$results$eff0.00_dur02$inc_cum_treated,  matrix(0, 2, 2),
               tolerance = 0.1)
  expect_true(all(abs(unlist(zz$results$vaccinated) - 1200 * 10 * 0.5) < 1e-5))
  expect_equal(length(zz$full_results), nrow(zz$inputs$grid))

  ## test with vod only
  zz2 <- run_grid(gp, ip, cp, blv,
                 model = run_onevax_xvwv,
                 strategy = "VoD",
                 eff = c(0, 1), dur = c(1, 2), vbe = 0.5,
                 uptake_total = 1)
  expect_equal(zz2$results$eff0.00_dur01$inc_treated, matrix(0, 2, 2),
               tolerance = 0.1)
  expect_equal(zz2$results$eff0.00_dur02$inc_treated,
               zz2$results$eff0.00_dur01$inc_treated,
               tolerance = 0.1)
  expect_true(all(unlist(zz2$results$inc_treated) > -0.1))
  expect_equal(zz2$results$eff0.00_dur01$inc_cum_treated, matrix(0, 2, 2),
               tolerance = 0.1)
  expect_equal(zz2$results$eff0.00_dur02$inc_cum_treated,
               zz2$results$eff0.00_dur01$inc_cum_treated,
               tolerance = 0.1)

  ## test with va only
  zz3 <- run_grid(gp, ip, cp, blv,
                  model = run_onevax_xvwv,
                  strategy = "VoA",
                  eff = c(0, 1), dur = c(1, 2), vbe = 0.5,
                  uptake_total = 1, full_output = TRUE)
  expect_equal(zz3$results$eff0.00_dur01$inc_treated, matrix(0, 2, 2),
               tolerance = 0.1)
  expect_equal(zz3$results$eff0.00_dur02$inc_treated,
               zz3$results$eff0.00_dur01$inc_treated,
               tolerance = 0.1)
  expect_true(all(unlist(zz3$results$inc_treated) > -0.1))
  expect_equal(zz3$results$eff0.00_dur01$inc_cum_treated, matrix(0, 2, 2),
               tolerance = 0.1)
  expect_equal(zz3$results$eff0.00_dur02$inc_cum_treated,
               zz3$results$eff0.00_dur01$inc_cum_treated,
               tolerance = 0.1)

  ## test with vt only
  zz5 <- run_grid(gp, ip, cp, blv,
                 model = run_onevax_xvwv,
                 strategy = "VoD(L)+VoA(H)",
                 eff = c(0, 1), dur = c(1, 99), vbe = 0.5,
                 uptake_total = 1, full_output = TRUE)
  expect_equal(zz5$results$eff0.00_dur01$inc_treated, matrix(0, 2, 2),
               tolerance = 0.1)
  expect_equal(zz5$results$eff0.00_dur99$inc_treated,
               zz5$results$eff0.00_dur01$inc_treated,
               tolerance = 0.1)
  expect_true(all(unlist(zz5$results$inc_treated) > -0.1))
  expect_equal(zz5$results$eff0.00_dur01$inc_cum_treated, matrix(0, 2, 2),
               tolerance = 0.1)
  expect_equal(zz5$results$eff0.00_dur99$inc_cum_treated,
               zz5$results$eff0.00_dur01$inc_cum_treated,
               tolerance = 0.1)
  expect_true(all(zz5$full_results[[1]][[1]]$cum_vaccinated[, 1, ] <=
           zz3$full_results[[1]][[1]]$cum_vaccinated[, 1, ]))

  tmp <- run_onevax_xvwv(tt, gp[1], init_params = ip[1], vea = 1, dur = 99,
                         vbe = 0.5, strategy = "VoD(L)+VoA(H)",
                         uptake = 1, t_stop = 99)
  tmp2 <- run_onevax_xvwv(tt, gp[1], init_params = ip[1], vea = 0, dur = 1,
                         vbe = 0.5, strategy = "VoD(L)+VoA(H)",
                         uptake = 1, t_stop = 99)
  expect_equal(tmp[[1]], zz5$full_results[[4]][[1]])
  expect_equal(tmp2[[1]], zz5$full_results[[1]][[1]])

})
