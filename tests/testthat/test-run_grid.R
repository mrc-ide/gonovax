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

test_that("set_strategy works as expected", {

  i <- 0.234
  expect_equal(set_strategy("VbE", i), list(vd = 0, vs = 0))
  expect_equal(set_strategy("VoD", i), list(vd = i, vs = 0))
  expect_equal(set_strategy("VoD(H)", i), list(vd = c(0, i), vs = 0))
  expect_equal(set_strategy("VoA", i), list(vd = i, vs = i))
  expect_equal(set_strategy("VoA(H)", i), list(vd = c(0, i), vs = c(0, i)))
  expect_equal(set_strategy("VoD(L)+VoA(H)", i), list(vd = i, vs = c(0, i)))
  expect_equal(set_strategy("VoS", i), list(vd = 0, vs = i))

  expect_error(set_strategy("hello", i), "strategy not recognised")
  expect_error(set_strategy("VbE", c(i, i)), "uptake must be length 1")

})


test_that("compare baseline works as expected", {

  gp <- gono_params(1:2)
  ip <- lapply(run_onevax_xvwv(0:1, gp, eff = 0, dur = 1), restart_params)
  tt <- 1:4

  bl <- extract_flows(run_onevax_xvwv(tt, gp, ip, eff = 0, dur = 1, ve = 0.5))
  blv <- rep(list(bl), 4)
  cp <- list(qaly_loss_per_diag_s = c(0.002, 0.001),
             unit_cost_manage_symptomatic = c(98, 99),
             unit_cost_manage_asymptomatic = c(93, 94),
             unit_cost_screen_uninfected = c(70, 71))
  p <- 0.7

  y <- run_onevax_xvwv(tt, gp, ip, eff = 0.5, dur = 1, ve = 0.5,
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
    (x$vaccinated - x$revaccinated) * (1 + 1 / p) + x$revaccinated
  }

  # check doses are calculated correctly
  expect_equal(f(yy) - f(bl), z$inc_doses)
  expect_equal(z$cases_averted_per_dose[1, ],
               -z$inc_treated[1, ] / z$inc_doses[1, ])
  expect_equal(z$cases_averted_per_dose[length(tt) - 1, ],
               -colSums(z$inc_treated) / colSums(z$inc_doses))
  ## with zero disc rate
  expect_equal(z$cases_averted_per_dose, z$cases_averted_per_dose_pv)

  tmp <- calc_cases_averted_per_dose(z, 0.03)
  expect_equal(z$cases_averted_per_dose[1, ], tmp[1, ])
  expect_true(all(tmp[2, ] != z$cases_averted_per_dose[1, ]))

  ## test cost eff threshold
  expect_equal(z$cet_20k, calc_cost_eff_threshold(20000, z, cp, 0))
  expect_equal(z$cet_30k, calc_cost_eff_threshold(30000, z, cp, 0))

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

})


test_that("run_grid works as expected", {
  ## test with ve only

  gp <- gono_params(1:2)
  ip <- lapply(run_onevax_xvwv(0:1, gp, eff = 0, dur = 1), restart_params)
  tt <- 1:3

  bl <- extract_flows(run_onevax_xvwv(tt, gp, ip, eff = 0, dur = 1))
  blv <- rep(list(bl), 4)
  cp <- list(qaly_loss_per_diag_s = c(0.002, 0.001),
             unit_cost_manage_symptomatic = c(98, 99),
             unit_cost_manage_asymptomatic = c(93, 94),
             unit_cost_screen_uninfected = c(70, 71))

  y <- run_onevax_xvwv(tt, gp, ip, eff = 0, dur = 1, ve = 0.5)
  z <- compare_baseline(y, bl, 0.7, cp, 0)

  zz <- run_grid(gp, ip, cp, blv,
                model = run_onevax_xvwv,
                strategy = "VbE",
                eff = c(0, 1), dur = c(1, 2), ve = 0.5,
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

  ## test with vd only
  zz2 <- run_grid(gp, ip, cp, blv,
                 model = run_onevax_xvwv,
                 strategy = "VoD",
                 eff = c(0, 1), dur = c(1, 2), ve = 0.5,
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
                  eff = c(0, 1), dur = c(1, 2), ve = 0.5,
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
                 eff = c(0, 1), dur = c(1, 99), ve = 0.5,
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

  tmp <- run_onevax_xvwv(tt, gp[1], init_params = ip[1], eff = 1, dur = 99,
                         ve = 0.5, strategy = "VoD(L)+VoA(H)",
                         uptake = 1, t_stop = 99)
  tmp2 <- run_onevax_xvwv(tt, gp[1], init_params = ip[1], eff = 0, dur = 1,
                         ve = 0.5, strategy = "VoD(L)+VoA(H)",
                         uptake = 1, t_stop = 99)
  expect_equal(tmp[[1]], zz5$full_results[[4]][[1]])
  expect_equal(tmp2[[1]], zz5$full_results[[1]][[1]])

})
