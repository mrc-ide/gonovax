
test_that("run_grid works as expected", {
  y <- run_grid(n = 2, t = 2, eff = c(0, 1), dur = c(1, 2), ve = 0.5,
                       full_output = TRUE)
  expect_equal(y$red_incid[, 1], c(0, 0), tol = 0.1)
  expect_equal(y$red_incid[, 1], y$red_incid[, 3], tol = 0.1)
  expect_true(all(y$red_incid > -0.1))
  expect_equal(y$cum_red_incid[, 1], c(0, 0), tol = 0.1)
  expect_equal(y$red_incid[, 1], y$red_incid[, 3], tol = 0.1)
  expect_true(all(abs(y$cum_vaccinated - 1200 * 10 * 0.5 * 2) < 1e-5))
  expect_equal(length(y$results), nrow(y$inputs$grid))
  plot(y)

  y2 <- run_grid(n = 2, t = 2, eff = c(0, 1), dur = c(1, 2),
                 strategy = "vd", uptake = 1)
  expect_equal(y2$red_incid[, 1], c(0, 0), tol = 0.1)
  expect_equal(y2$red_incid[, 1], y2$red_incid[, 3], tol = 0.1)
  expect_true(all(y2$red_incid > -0.1))
  expect_equal(y2$cum_red_incid[, 1], c(0, 0), tol = 0.1)
  expect_equal(y2$red_incid[, 1], y2$red_incid[, 3], tol = 0.1)
  plot(y2)
  ### this looks fine
  y3 <- run_grid(n = 2, t = 2, eff = c(0, 1), dur = c(1, 2),
                 strategy = "va", uptake = 1, full_output = TRUE)
  expect_equal(y3$red_incid[, 1], c(0, 0), tol = 0.1)
  expect_equal(y3$red_incid[, 1], y3$red_incid[, 3], tol = 0.1)
  expect_equal(y3$cum_red_incid[, 1], c(0, 0), tol = 0.1)
  expect_equal(y3$cum_red_incid[, 1], y3$cum_red_incid[, 3], tol = 0.1)
  expect_true(all(y3$red_incid > -0.1))
  plot(y3)

  y5 <- run_grid(n = 2, t = 2, eff = c(0, 1), dur = c(1, 2),
                 strategy = "vt", uptake = 1, full_output = TRUE)
  expect_equal(y5$red_incid[, 1], c(0, 0), tol = 0.1)
  expect_equal(y5$red_incid[, 1], y5$red_incid[, 3], tol = 0.1)
  expect_equal(y5$cum_red_incid[, 1], c(0, 0), tol = 0.1)
  expect_equal(y5$cum_red_incid[, 1], y5$cum_red_incid[, 3], tol = 0.1)
  expect_true(all(y5$red_incid > -0.1))
  expect_true(all(y5$results[[1]][[1]]$cum_vaccinated[, 1, ] <=
           y3$results[[1]][[1]]$cum_vaccinated[, 1, ]))

  y0 <- run_grid(n = 2, t = 2, eff = c(0, 1), dur = c(1, 2))
  # runs ok with user baseline
  y4 <- run_grid(n = 2, t = 2, eff = c(0, 1), dur = c(1, 2), ve = 0.5,
                        baseline = y0)
  expect_equal(y4$red_incid, y$red_incid)
  expect_equal(y4$cum_red_incid, y$cum_red_incid)

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
  expect_equal(y$red_incid[, 1], c(0, 0), tol = 0.1)
  expect_equal(y$red_incid[, 1], y$red_incid[, 3], tol = 0.1)
  expect_true(all(y$red_incid > -0.1))
  expect_equal(y$cum_red_incid[, 1], c(0, 0), tol = 0.1)
  expect_equal(y$red_incid[, 1], y$red_incid[, 3], tol = 0.1)
  expect_true(all(abs(y$cum_vaccinated - 1200 * 10 * 0.5 * 2) < 1e-5))
  expect_equal(length(y$results), nrow(y$inputs$grid))
  plot(y)
})
