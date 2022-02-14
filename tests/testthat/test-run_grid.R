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

test_that("compare baseline works as expected", {

  gp <- gono_params(1:2)
  ip <- lapply(run_onevax_xvwv(0:1, gp, vea = 0, dur = 1), restart_params)
  tt <- 1:4

  bl <- extract_flows(run_onevax_xvwv(tt, gp, ip, vea = 0, dur = 1, vbe = 0.5))
  blv <- rep(list(bl), 4)
  cp <- list(qaly_loss_per_diag_s = c(0.002, 0.001),
             unit_cost_manage_symptomatic = c(98, 99),
             unit_cost_manage_asymptomatic = c(93, 94),
             unit_cost_screen_uninfected = c(70, 71))
  p <- 0.7

  y <- run_onevax_xvwv(tt, gp, ip, vea = 0.5, dur = 1, vbe = 0.5,
                       uptake = 0.75, strategy = "VoD")
  yy <- extract_flows(y)

  z <- compare_baseline(y, bl, uptake_second_dose = p, cp, 0)

  flownames <- names(yy)
  incnames <- paste0("inc_", flownames)
  compnames <- setdiff(names(z), c(flownames, incnames))

  # check values are carried over correctly
  expect_equal(z[flownames], yy)
  expect_equivalent(z[incnames],
                    Map(`-`, yy[flownames], bl[flownames]))

  f <- function(x) {
    pv <- x$vaccinated - x$revaccinated - x$vbe
    pv * (1 + 1 / p) + x$revaccinated + 2 * x$vbe
  }

  # check doses are calculated correctly
  n <- length(tt) - 1
  expect_equal(f(yy) - f(bl), z$inc_doses)
  expect_equal(z$inc_cum_doses[1, ], z$inc_doses[1, ])
  expect_equal(z$inc_cum_doses[n, ], colSums(z$inc_doses))

  expect_equal(z$cases_averted_per_dose[1, ],
               -z$inc_treated[1, ] / z$inc_doses[1, ])
  expect_equal(z$cases_averted_per_dose[n, ],
               -colSums(z$inc_treated) / colSums(z$inc_doses))

  ## check correct with double dose at revaccination
  expect_equal(calc_doses(z, uptake_second_dose = p, revax_one_dose = FALSE),
               z$inc_vaccinated * (1 + 1 / p))

  ## check against a baseline of no vaccination
  bl0 <- extract_flows(run_onevax_xvwv(tt, gp, ip, vea = 0, dur = 1, vbe = 0))
  blv0 <- rep(list(bl0), 4)

  z1 <- compare_baseline(y, bl0, uptake_second_dose = p, cp, 0)
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
                    matrix(c(0.793537914637929, 2.44111479409769,
                             4.7985238475747, 0.762611923450001,
                             2.42625463089821, 4.85085764438457), nrow = 3))
  expect_equal(cet_30k,
               matrix(c(0.906458889656841, 2.78510557322961,
                        5.47131813910508, 0.826552593480719,
                        2.62838651575972, 5.25377713095235), nrow = 3))

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

  expect_equal(calc_inc_costs(85, cost) - calc_inc_costs(18, cost),
               cost$pv_inc_doses * (85 - 18))

  ## test doses wasted calculated correctly

      # no doses wasted when second dose uptake dose 1
  z <- compare_baseline(y, bl, uptake_second_dose = 1, cp, 0)

  expect_equal(z$inc_cum_doses_wasted, matrix(rep(0, 6), 3, 2))
  expect_equal(z$inc_doses_wasted, matrix(rep(0, 6), 3, 2))

        # regardless of first dose uptake

        # first dose = 0.1, second dose = 1
  gp <- gonovax::gono_params(1:2)
  ip <- lapply(run_onevax_xvwrh(0:1, gp, vea = 0, dur = 1), restart_hes)
  tt <- 1:4

  bl <- extract_flows(run_onevax_xvwrh(tt, gp, ip, vea = 0, dur = 1))

  y_f <- run_onevax_xvwrh(tt, gp, ip, vea = 0.5, dur = 1, vbe = 0,
                        primary_uptake = 0.5,
                        strategy = "VoD")

  z_f <- compare_baseline(y_f, bl, uptake_second_dose = 1, cp, 0)

  expect_equal(z_f$inc_cum_doses_wasted, matrix(rep(0, 6), 3, 2))
  expect_equal(z_f$inc_doses_wasted, matrix(rep(0, 6), 3, 2))

        # regardless of hesitancy (primary uptake = booster uptake = 1)
  ip_hes <- lapply(run_onevax_xvwrh(0:1, gp, vea = 0, dur = 1),
               restart_hes, hes = 0.3)

  bl_hes <- extract_flows(run_onevax_xvwrh(tt, gp, ip_hes, vea = 0, dur = 1))

  y_h <- run_onevax_xvwrh(tt, gp, ip_hes, vea = 0.5, dur = 1, vbe = 0.5,
                        primary_uptake = 1, strategy = "VoD", hes = 0.3)

  z_h <- compare_baseline(y, bl, uptake_second_dose = 1, cp, 0)

  expect_equal(z_h$inc_cum_doses_wasted, matrix(rep(0, 6), 3, 2))
  expect_equal(z_h$inc_doses_wasted, matrix(rep(0, 6), 3, 2))

      # all primary doses wasted when second dose uptake = 0

  z_f2 <- compare_baseline(y_f, bl, uptake_second_dose = 0, cp, 0)

      # (when second_uptake = 0, primary doses pp is infinitely large)
      # i.e any and all possible numbers of doses are wasted

  expect_equal(z_f2$inc_cum_doses_wasted, matrix(rep(Inf, 6), 3, 2))

    # doses used + doses wasted = 2 * #primary doses + 1 * booster doses
                                                    # + 2 * vbe doses

  a <- z_h$inc_cum_doses + z_h$inc_cum_doses_wasted
  b <- 2 * z_h$inc_cum_primary + 1 * z_h$inc_cum_revaccinated +
       2 * z_h$inc_cum_vbe
  expect_equal(a, b)

  ## test that cumulative primary and booster vaccination calc correctly
  gp <- gono_params(1:2)

  cp <- list(qaly_loss_per_diag_s = c(0.002, 0.001),
             unit_cost_manage_symptomatic = c(98, 99),
             unit_cost_manage_asymptomatic = c(93, 94),
             unit_cost_screen_uninfected = c(70, 71))

  ip <- lapply(run_onevax_xvwrh(0:1, gp, vea = 0, dur = 1), restart_hes)
  tt <- 1:4

  bl <- extract_flows(run_onevax_xvwrh(tt, gp, ip, vea = 0, dur = 1e50))

  y <- run_onevax_xvwrh(tt, gp, ip, vea = 1, dur = 1, vbe = 0,
                          primary_uptake = 1, booster_uptake = 0,
                          strategy = "VoD")
  z <- compare_baseline(y, bl, uptake_second_dose = 1, cp, 0)

  # number undergoing primary vaccination is equal to total people
  # vaccinated when booster uptake = 0

    expect_equal(z$inc_cum_primary[length(tt) - 1, 1],
               z$inc_cum_vaccinated[length(tt) - 1, 1])

  # number booster vaccinations is equal to
      # inc_cum_vaccinated - inc_cum_primary
      # sum re-vaccinated

  y <- run_onevax_xvwrh(tt, gp, ip, vea = 1, dur = 1, vbe = 0,
                        primary_uptake = 1, booster_uptake = 0.5,
                        strategy = "VoD")
  z <- compare_baseline(y, bl, uptake_second_dose = 1, cp, 0)
  
  expect_equal(z$inc_cum_vaccinated[length(tt) -1, 1] -
                 z$inc_cum_primary[length(tt) - 1, 1],
               z$inc_cum_revaccinated[length(tt) - 1, 1])
  
  expect_equal(sum(z$inc_revaccinated[, 1]),
               z$inc_cum_revaccinated[length(tt) - 1, 1])

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
  z <- compare_baseline(y, bl, 0.7, cp, 0)

  zz <- run_grid(gp, ip, cp, blv,
                model = run_onevax_xvwv,
                strategy = "VbE",
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
