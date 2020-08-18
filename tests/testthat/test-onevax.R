test_that("run_onevax works correctly", {
  tt <- seq(0, 5)
  y1 <- run_onevax(1, tt, eff = 0, dur = 1e3)[[1]]
  # check no-one is vaccinated with v switched off
  expect_error(run_onevax_int(1:2, tt, equilib = TRUE),
               label = "if length(n) > 1, equilib must be FALSE")
  expect_true(all(y1$cum_vaccinated == 0))
  y2 <- run_onevax(1, tt, eff = 0, dur = 1e3, ve = 1)[[1]]
  # check 100% vbe vaccinates all new entrants
  expect_equal(diff(rowSums(y2$cum_vaccinated[, , , 2])), rep(12e3, max(tt)))
  expect_equal(sum(y2$cum_vaccinated[, , , 1]), 0)

  # check can run from equilib
  y2e <- run_onevax(1, tt, eff = 0, dur = 1e3, ve = 1, equilib = TRUE)[[1]]

  expect_equal(y2$cum_vaccinated, y2e$cum_vaccinated)
  expect_equal(apply(y2e$N, 1, sum), rep(6e5, 6), tol = 1e-5)

  # check vaccination on treatment is working correctly
  y3e <- run_onevax(1, tt, eff = 1, dur = 1e3, vd = 1, equilib = TRUE)[[1]]
  expect_equal(y3e$cum_vaccinated[, , , 1], y3e$cum_treated[, , , 1])
  expect_equal(sum(y3e$cum_treated[, , , 2]), 0)
  expect_equal(apply(y3e$N, 1, sum), rep(6e5, 6), tol = 1e-5)

  # check vaccination on screening is working correctly
  y4e <- run_onevax(1, tt, eff = 1, dur = 1e3, vs = 1, equilib = TRUE)[[1]]
  expect_equal(y4e$cum_vaccinated[, , , 1], y4e$cum_screened[, , , 1])
  expect_equal(apply(y4e$N, 1, sum), rep(6e5, 6), tol = 1e-5)
})

test_that("run_grid_onevax works as expected", {
  y <- run_grid_onevax(n = 2, t = 2, eff = c(0, 1), dur = c(1, 2), ve = 0.5)
  expect_equal(y$red_incid[, 1], c(0, 0), tol = 0.1)
  expect_equal(y$red_incid[, 1], y$red_incid[, 3], tol = 0.1)
  expect_true(all(y$red_incid > -0.1))
  expect_equal(y$cum_red_incid[, 1], c(0, 0), tol = 0.1)
  expect_equal(y$red_incid[, 1], y$red_incid[, 3], tol = 0.1)
  expect_true(all(abs(y$cum_vaccinated - 1200 * 10 * 0.5 * 2) < 1e-5))
  plot(y)
  y2 <- run_grid_onevax(n = 2, t = 2, eff = c(0, 1), dur = c(1, 2), vd = 1)
  expect_equal(y2$red_incid[, 1], c(0, 0), tol = 0.1)
  expect_equal(y2$red_incid[, 1], y2$red_incid[, 3], tol = 0.1)
  expect_true(all(y2$red_incid > -0.1))
  expect_equal(y2$cum_red_incid[, 1], c(0, 0), tol = 0.1)
  expect_equal(y2$red_incid[, 1], y2$red_incid[, 3], tol = 0.1)
  plot(y2)
   ### this looks fine
  y3 <- run_grid_onevax(n = 2, t = 2, eff = c(0, 1), dur = c(1, 2), vs = 1)
  expect_equal(y3$red_incid[, 1], c(0, 0), tol = 0.1)
  expect_equal(y3$red_incid[, 1], y3$red_incid[, 3], tol = 0.1)
  expect_equal(y3$cum_red_incid[, 1], c(0, 0), tol = 0.1)
  expect_equal(y3$cum_red_incid[, 1], y3$cum_red_incid[, 3], tol = 0.1)
  expect_true(all(y3$red_incid > -0.1))
  plot(y3)
})
