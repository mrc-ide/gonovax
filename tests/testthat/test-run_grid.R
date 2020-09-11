
test_that("run_grid works as expected", {
  ## test with ve only
  y <- run_grid(n = 2, t = 2, eff = c(0, 1), dur = c(1, 2), ve = 0.5,
                       full_output = TRUE)
  expect_equal(y$eff0.00_dur01$inc_incid, matrix(0, 2, 2), tol = 0.1)
  expect_equal(y$eff0.00_dur02$inc_incid, y$eff0.00_dur01$inc_incid, tol = 0.1)
  expect_true(all(unlist(y$inc_incid$eff0.10_dur01) > -0.1))
  expect_equal(y$eff0.00_dur01$inc_cum_incid,  matrix(0, 2, 2), tol = 0.1)
  expect_equal(y$eff0.00_dur02$inc_cum_incid,  matrix(0, 2, 2), tol = 0.1)
  expect_true(all(abs(unlist(y$vaccinated) - 1200 * 10 * 0.5) < 1e-5))
  expect_equal(length(y$full_results), nrow(y$inputs$grid))

  ## test with vd only
  y2 <- run_grid(n = 2, t = 2, eff = c(0, 1), dur = c(1, 2),
                 strategy = "vd", uptake = 1)
  expect_equal(y2$eff0.00_dur01$inc_incid, matrix(0, 2, 2), tol = 0.1)
  expect_equal(y2$eff0.00_dur02$inc_incid, y2$eff0.00_dur01$inc_incid,
               tol = 0.1)
  expect_true(all(unlist(y2$inc_incid) > -0.1))
  expect_equal(y2$eff0.00_dur01$inc_cum_incid, matrix(0, 2, 2), tol = 0.1)
  expect_equal(y2$eff0.00_dur02$inc_cum_incid, y2$eff0.00_dur01$inc_cum_incid,
               tol = 0.1)

  ## test with va only
  y3 <- run_grid(n = 2, t = 2, eff = c(0, 1), dur = c(1, 2),
                 strategy = "va", uptake = 1, full_output = TRUE)
  expect_equal(y3$eff0.00_dur01$inc_incid, matrix(0, 2, 2), tol = 0.1)
  expect_equal(y3$eff0.00_dur02$inc_incid, y3$eff0.00_dur01$inc_incid,
               tol = 0.1)
  expect_true(all(unlist(y3$inc_incid) > -0.1))
  expect_equal(y3$eff0.00_dur01$inc_cum_incid, matrix(0, 2, 2), tol = 0.1)
  expect_equal(y3$eff0.00_dur02$inc_cum_incid, y3$eff0.00_dur01$inc_cum_incid,
               tol = 0.1)

  ## test with vt only
  y5 <- run_grid(n = 2, t = 2, eff = c(0, 1), dur = c(1, 2),
                 strategy = "vt", uptake = 1, full_output = TRUE)
  expect_equal(y5$eff0.00_dur01$inc_incid, matrix(0, 2, 2), tol = 0.1)
  expect_equal(y5$eff0.00_dur02$inc_incid, y5$eff0.00_dur01$inc_incid,
               tol = 0.1)
  expect_true(all(unlist(y5$inc_incid) > -0.1))
  expect_equal(y5$eff0.00_dur01$inc_cum_incid, matrix(0, 2, 2), tol = 0.1)
  expect_equal(y5$eff0.00_dur02$inc_cum_incid, y5$eff0.00_dur01$inc_cum_incid,
               tol = 0.1)
  expect_true(all(y5$full_results[[1]][[1]]$cum_vaccinated[, 1, ] <=
           y3$full_results[[1]][[1]]$cum_vaccinated[, 1, ]))

  # test with with user baseline
  y0 <- run_grid(n = 2, t = 2, eff = c(0, 1), dur = c(1, 2))
  y4 <- run_grid(n = 2, t = 2, eff = c(0, 1), dur = c(1, 2), ve = 0.5,
                        baseline = y0)
  expect_equivalent(y4, y[-6], tol = 0.1)



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
  expect_equal(y$eff0.00_dur01$inc_incid, matrix(0, 2, 2), tol = 0.1)
  expect_equal(y$eff0.00_dur02$inc_incid, y$eff0.00_dur01$inc_incid, tol = 0.1)
  expect_true(all(unlist(y$inc_incid) > -0.1))
  expect_equal(y$eff0.00_dur01$inc_cum_incid,  matrix(0, 2, 2), tol = 0.1)
  expect_equal(y$eff0.00_dur02$inc_cum_incid,  matrix(0, 2, 2), tol = 0.1)
  expect_true(all(abs(unlist(y$vaccinated) - 1200 * 10 * 0.5) < 1e-5))
  expect_equal(length(y$full_results), nrow(y$inputs$grid))

})
