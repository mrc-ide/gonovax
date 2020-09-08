
test_that("run_grid works as expected", {
  ## test with ve only
  y <- run_grid(n = 2, t = 2, eff = c(0, 1), dur = c(1, 2), ve = 0.5,
                       full_output = TRUE)
  expect_equal(y$red_incid$eff0.00_dur01, matrix(0, 2, 2), tol = 0.1)
  expect_equal(y$red_incid$eff0.00_dur02, y$red_incid$eff0.00_dur01, tol = 0.1)
  expect_true(all(unlist(y$red_incid) > -0.1))
  expect_equal(y$red_cum_incid$eff0.00_dur01,  matrix(0, 2, 2), tol = 0.1)
  expect_equal(y$red_cum_incid$eff0.00_dur02,  matrix(0, 2, 2), tol = 0.1)
  expect_true(all(abs(unlist(y$vaccinated) - 1200 * 10 * 0.5) < 1e-5))
  expect_equal(length(y$full_results), nrow(y$inputs$grid))

  ## test with vd only
  y2 <- run_grid(n = 2, t = 2, eff = c(0, 1), dur = c(1, 2),
                 strategy = "vd", uptake = 1)
  expect_equal(y2$red_incid$eff0.00_dur01, matrix(0, 2, 2), tol = 0.1)
  expect_equal(y2$red_incid$eff0.00_dur02, y2$red_incid$eff0.00_dur01,
               tol = 0.1)
  expect_true(all(unlist(y2$red_incid) > -0.1))
  expect_equal(y2$red_cum_incid$eff0.00_dur01, matrix(0, 2, 2), tol = 0.1)
  expect_equal(y2$red_cum_incid$eff0.00_dur02, y2$red_cum_incid$eff0.00_dur01,
               tol = 0.1)

  ## test with va only
  y3 <- run_grid(n = 2, t = 2, eff = c(0, 1), dur = c(1, 2),
                 strategy = "va", uptake = 1, full_output = TRUE)
  expect_equal(y3$red_incid$eff0.00_dur01, matrix(0, 2, 2), tol = 0.1)
  expect_equal(y3$red_incid$eff0.00_dur02, y3$red_incid$eff0.00_dur01,
               tol = 0.1)
  expect_true(all(unlist(y3$red_incid) > -0.1))
  expect_equal(y3$red_cum_incid$eff0.00_dur01, matrix(0, 2, 2), tol = 0.1)
  expect_equal(y3$red_cum_incid$eff0.00_dur02, y3$red_cum_incid$eff0.00_dur01,
               tol = 0.1)

  ## test with vt only
  y5 <- run_grid(n = 2, t = 2, eff = c(0, 1), dur = c(1, 2),
                 strategy = "vt", uptake = 1, full_output = TRUE)
  expect_equal(y5$red_incid$eff0.00_dur01, matrix(0, 2, 2), tol = 0.1)
  expect_equal(y5$red_incid$eff0.00_dur02, y5$red_incid$eff0.00_dur01,
               tol = 0.1)
  expect_true(all(unlist(y5$red_incid) > -0.1))
  expect_equal(y5$red_cum_incid$eff0.00_dur01, matrix(0, 2, 2), tol = 0.1)
  expect_equal(y5$red_cum_incid$eff0.00_dur02, y5$red_cum_incid$eff0.00_dur01,
               tol = 0.1)
  expect_true(all(y5$full_results[[1]][[1]]$cum_vaccinated[, 1, ] <=
           y3$full_results[[1]][[1]]$cum_vaccinated[, 1, ]))

  # test with with user baseline
  y0 <- run_grid(n = 2, t = 2, eff = c(0, 1), dur = c(1, 2))
  y4 <- run_grid(n = 2, t = 2, eff = c(0, 1), dur = c(1, 2), ve = 0.5,
                        baseline = y0)
  expect_equal(y4$red_incid, y$red_incid, tol = 0.1)
  expect_equal(y4$cum_red_incid, y$cum_red_incid, tol = 0.1)

  # check error cases
  expect_error(run_grid(n = 3, t = 2, eff = c(0, 1), dur = c(1, 2),
                               ve = 0.5, baseline = y0),
               "model parameters do not match baseline")
  expect_error(run_grid(n = 2, t = 2, eff = c(0, 1), dur = c(1, 3),
                               ve = 0.5, baseline = y0),
               "dur / eff parameters do not match baseline")
  expect_error(run_grid(n = 2, t = 3, eff = c(0, 1), dur = c(1, 2),
                               ve = 0.5, baseline = y0),
               "t does not match baseline")
  class(y0) <- NULL
  expect_error(run_grid(n = 2, t = 2, eff = c(0, 1), dur = c(1, 2),
                               ve = 0.5, baseline = y0),
               "baseline must be a gonovax_grid object")

})

test_that("run_grid works with onevax_waning model", {
  y <- run_grid(n = 2, t = 2, model = run_onevax_waning,
                eff = c(0, 1), dur = c(1, 2), ve = 0.5,
                full_output = TRUE)
  expect_equal(y$red_incid$eff0.00_dur01, matrix(0, 2, 2), tol = 0.1)
  expect_equal(y$red_incid$eff0.00_dur02, y$red_incid$eff0.00_dur01, tol = 0.1)
  expect_true(all(unlist(y$red_incid) > -0.1))
  expect_equal(y$red_cum_incid$eff0.00_dur01,  matrix(0, 2, 2), tol = 0.1)
  expect_equal(y$red_cum_incid$eff0.00_dur02,  matrix(0, 2, 2), tol = 0.1)
  expect_true(all(abs(unlist(y$vaccinated) - 1200 * 10 * 0.5) < 1e-5))
  expect_equal(length(y$full_results), nrow(y$inputs$grid))

})

test_that("gonovax_grid format method works as expected", {
  t <- 2
  y <- run_grid(n = 2, t = t, model = run_onevax,
                eff = c(0.1, 1), dur = c(1, 2), ve = 0.5)
  z <- format_grid(y, f = mean)
  y0 <- run_novax(n = 1:2, tt = 0:2, equilib = TRUE)

  cum_incid0 <- sapply(y0, function(x) rowSums(x$cum_incid))
  cum_diag_a0 <- sapply(y0, function(x) rowSums(x$cum_diag_a))
  cum_diag_s0 <- sapply(y0, function(x) rowSums(x$cum_diag_s))
  cum_screened0 <- sapply(y0, function(x) rowSums(x$cum_screened))
  incid0 <- apply(cum_incid0, 2, diff)

  # check heatmaps are compiled correctly
  expect_equivalent(z$ts$eff0.10_dur01$red_incid,
                    rowMeans(incid0 - y$incid[[1]]), tol = 0.1)
  expect_equivalent(z$ts$eff0.10_dur01$inc_cum_vaccinated,
                    rowMeans(y$cum_vaccinated$eff0.10_dur01))
  tot_red_incid <- cum_incid0[-1, ] - y$cum_incid[[1]]
  expect_equivalent(z$ts$eff0.10_dur01$red_cum_incid,
                    rowMeans(tot_red_incid), tol = 0.1)
  expect_equivalent(z$ts$eff0.10_dur01$cost_eff,
               rowMeans(y$cum_vaccinated[[1]] / tot_red_incid),
               tol = 0.1)
  expect_equivalent(z$ts$eff0.10_dur01$red_cum_diag_a,
                    rowMeans(cum_diag_a0[-1, ] - y$cum_diag_a[[1]]))
  expect_equivalent(z$ts$eff0.10_dur01$red_cum_diag_s,
                    rowMeans(cum_diag_s0[-1, ] - y$cum_diag_s[[1]]))

  # check discount rate
  z1 <- format_grid(y, 0.05, mean)
  w <- which(names(z) == "cost_eff")
  expect_equal(z$ts$eff0.10_dur01[-w], z1$ts$eff0.10_dur01[-w])

  pv <- 1.05 ^ -c(0.5, 1.5)
  pv_inc_vaccinated <- sapply(y$inc_vaccinated, function(x) colSums(x * pv))
  pv_red_incid <- sapply(y$red_incid, function(x) colSums(x * pv))
  expect_equivalent(z1$ts$eff0.10_dur01$cost_eff[1],
                    z$ts$eff0.10_dur01$cost_eff[1])
  expect_equivalent(z1$ts$eff0.10_dur01$cost_eff[2],
                    colMeans(pv_inc_vaccinated / pv_red_incid)[1])

  # check error case
  class(y) <- NULL
  expect_error(format_grid(y, 0, mean))
})
