test_that("calc_pv works as expected", {

  ## time is dimension 1
  x <- matrix(seq_len(12), 3, 2)
  y0 <- calc_pv(x, 0)

  expect_equal(y0, apply(x, 2, cumsum))

  d <- 0.1
  disc_fac <- (1 + d) ^ -(1:3 - 0.5)

  y1 <- calc_pv(x, d)
  expect_equal(y1[1, ], x[1, ] * disc_fac[1])
  expect_equal(y1[3, ], c(disc_fac %*% x))
})

test_that("set_strategy works as expected", {

  i <- 0.234
  expect_equal(set_strategy("VbE", i), list(vd = 0, vs = 0))
  expect_equal(set_strategy("VoD(all)", i), list(vd = i, vs = 0))
  expect_equal(set_strategy("VoD(H)", i), list(vd = c(0, i), vs = 0))
  expect_equal(set_strategy("VoA(all)", i), list(vd = i, vs = i))
  expect_equal(set_strategy("VoA(H)", i), list(vd = c(0, i), vs = c(0, i)))
  expect_equal(set_strategy("VoD(L)+VoA(H)", i), list(vd = i, vs = c(0, i)))
  
  expect_error(set_strategy("hello", i), "strategy not recognised")
  expect_error(set_strategy("VbE", c(i, i)), "uptake must be length 1")

  
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
  z
  
  zz <- run_grid(t = 2, gp, ip, cp, blv,
                model = run_onevax_xvwv,
                strategy = "VbE", 
                eff = c(0, 1), dur = c(1, 2), ve = 0.5,
                uptake_total = 1,
                full_output = TRUE)
  expect_equal(zz$results$eff0.00_dur01$inc_incid, matrix(0, 2, 2), tolerance = 0.1)
  expect_equal(zz$results$eff0.00_dur02$inc_incid, zz$results$eff0.00_dur01$inc_incid,
               tolerance = 0.1)
  expect_true(all(unlist(zz$results$inc_incid$eff0.10_dur01) > -0.1))
  expect_equal(zz$results$eff0.00_dur01$inc_cum_incid,  matrix(0, 2, 2), tolerance = 0.1)
  expect_equal(zz$results$eff0.00_dur02$inc_cum_incid,  matrix(0, 2, 2), tolerance = 0.1)
  expect_true(all(abs(unlist(zz$results$vaccinated) - 1200 * 10 * 0.5) < 1e-5))
  expect_equal(length(zz$full_results), nrow(zz$inputs$grid))

  ## test with vd only
  zz2 <- run_grid(t = 2, gp, ip, cp, blv,
                 model = run_onevax_xvwv,
                 strategy = "VoD(all)", 
                 eff = c(0, 1), dur = c(1, 2), ve = 0.5,
                 uptake_total = 1)
  expect_equal(zz2$results$eff0.00_dur01$inc_incid, matrix(0, 2, 2),
               tolerance = 0.1)
  expect_equal(zz2$results$eff0.00_dur02$inc_incid,
               zz2$results$eff0.00_dur01$inc_incid,
               tolerance = 0.1)
  expect_true(all(unlist(zz2$results$inc_incid) > -0.1))
  expect_equal(zz2$results$eff0.00_dur01$inc_cum_incid, matrix(0, 2, 2),
               tolerance = 0.1)
  expect_equal(zz2$results$eff0.00_dur02$inc_cum_incid,
               zz2$results$eff0.00_dur01$inc_cum_incid,
               tolerance = 0.1)

  ## test with va only
  zz3 <- run_grid(t = 2, gp, ip, cp, blv,
                  model = run_onevax_xvwv,
                  strategy = "VoA(all)", 
                  eff = c(0, 1), dur = c(1, 2), ve = 0.5,
                  uptake_total = 1, full_output = TRUE)
  expect_equal(zz3$results$eff0.00_dur01$inc_incid, matrix(0, 2, 2),
               tolerance = 0.1)
  expect_equal(zz3$results$eff0.00_dur02$inc_incid,
               zz3$results$eff0.00_dur01$inc_incid,
               tolerance = 0.1)
  expect_true(all(unlist(zz3$results$inc_incid) > -0.1))
  expect_equal(zz3$results$eff0.00_dur01$inc_cum_incid, matrix(0, 2, 2),
               tolerance = 0.1)
  expect_equal(zz3$results$eff0.00_dur02$inc_cum_incid,
               zz3$results$eff0.00_dur01$inc_cum_incid,
               tolerance = 0.1)

  ## test with vt only
  zz5 <- run_grid(t = 2, gp, ip, cp, blv,
                 model = run_onevax_xvwv,
                 strategy = "VoD(all)_VoA(H)", 
                 eff = c(0, 1), dur = c(1, 2), ve = 0.5,
                 uptake_total = 1, full_output = TRUE)
  expect_equal(zz5$results$eff0.00_dur01$inc_incid, matrix(0, 2, 2),
               tolerance = 0.1)
  expect_equal(zz5$results$eff0.00_dur02$inc_incid,
               zz5$results$eff0.00_dur01$inc_incid,
               tolerance = 0.1)
  expect_true(all(unlist(zz5$results$inc_incid) > -0.1))
  expect_equal(zz5$results$eff0.00_dur01$inc_cum_incid, matrix(0, 2, 2),
               tolerance = 0.1)
  expect_equal(zz5$results$eff0.00_dur02$inc_cum_incid,
               zz5$results$eff0.00_dur01$inc_cum_incid,
               tolerance = 0.1)
  expect_true(all(zz5$results$full_results[[1]][[1]]$cum_vaccinated[, 1, ] <=
           zz3$full_results[[1]][[1]]$cum_vaccinated[, 1, ]))

## TODO
  # check error cases
  # expect_error(run_grid(n = 3, t = 2, eff = c(0, 1), dur = c(1, 2),
  #                              ve = 0.5, baseline = y0),
  #              "model parameters do not match baseline")
  # expect_error(run_grid(n = 2, t = 2, eff = c(0, 1), dur = c(1, 3),
  #                              ve = 0.5, baseline = y0),
  #              "dur / eff parameters do not match baseline")
  # expect_error(run_grid(n = 2, t = 3, eff = c(0, 1), dur = c(1, 2),
  #                              ve = 0.5, baseline = y0),
  #              "t does not match baseline")
  # class(y0) <- NULL
  # expect_error(run_grid(n = 2, t = 2, eff = c(0, 1), dur = c(1, 2),
  #                              ve = 0.5, baseline = y0),
  #              "baseline must be a gonovax_grid object")

})

