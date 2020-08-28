test_that("run_onevax_waning works correctly", {
  tt <- seq(0, 5)
  y1 <- run_onevax_waning(1, tt, eff = 0, dur = 1)[[1]]
  expect_error(run_onevaxw_int(1:2, tt),
               label = "length(n) must equal 1")
  # check no-one is vaccinated with v switched off
  expect_true(all(y1$cum_vaccinated == 0))
  y2 <- run_onevax_waning(1, tt, eff = 0, dur = 1, ve = 1)[[1]]
  # check 100% vbe vaccinates all new entrants
  expect_equal(diff(rowSums(y2$cum_vaccinated[, , 2])), rep(12e3, max(tt)))
  expect_equal(sum(y2$cum_vaccinated[, , 1]), 0)

  # check can run from equilib
  y2e <- run_onevax_waning(1, tt, eff = 0, dur = 1, ve = 1,
                           equilib = TRUE)[[1]]

  expect_equal(y2$cum_vaccinated, y2e$cum_vaccinated)
  expect_equal(apply(y2e$N, 1, sum), rep(6e5, 6), tol = 1e-5)

  # check vaccination on treatment is working correctly
  y3e <- run_onevax_waning(1, tt, eff = 1, dur = 1, vd = 1,
                           equilib = TRUE)[[1]]
  expect_equal(y3e$cum_vaccinated[, , 1], y3e$cum_treated[, , 1])
  expect_equal(sum(y3e$cum_treated[, , 2]), 0)
  expect_equal(apply(y3e$N, 1, sum), rep(6e5, 6), tol = 1e-5)

  # check vaccination on screening is working correctly
  y4e <- run_onevax_waning(1, tt, eff = 1, dur = 1, vs = 1,
                           equilib = TRUE)[[1]]
  expect_equal(y4e$cum_vaccinated[, , 1], y4e$cum_screened[, , 1])
  expect_equal(apply(y4e$N, 1, sum), rep(6e5, 6), tol = 1e-5)
  # check that number of vaccinations is lower with no revaccination
  # this is fine for vs but, more complex for vd as incidence affects vax cov
  z <- run_onevax(1, tt, eff = 1, dur = 1, vs = 1,
                  equilib = TRUE)[[1]]
  expect_true(all(z$cum_vaccinated[, , 1] >= y4e$cum_vaccinated[, , 1]))

  # check waned vaccines move into correct compartment
  y5e <- run_onevax_waning(1, tt, eff = 1, dur = 1, ve = 1,
                           equilib = TRUE)[[1]]
  expect_equal(apply(y5e$N, 1, sum), rep(6e5, 6), tol = 1e-5)
  expect_true(all(apply(y5e$N[-1, , 5], 1, sum) > 0))
  # check that equivalent to waning into U
  z <- run_onevax(1, tt, eff = 1, dur = 1, ve = 1,
                         equilib = TRUE)[[1]]
  expect_equal(z$N[, , 1], y5e$N[, , 1] + y5e$N[, , 5], tol = 0.1)
  expect_equal(z$cum_incid[, , 1], y5e$cum_incid[, , 1] + y5e$cum_incid[, , 5],
               tol = 0.1)

  # check vaccination targeting
  y6e <- run_onevax_waning(1, tt, eff = 1, dur = 1,
                    vs = c(0, 1), vd = c(0, 1), equilib = TRUE)[[1]]
  expect_equal(y6e$N[, 1, 2], rep(0, 6))
  expect_true(all(y6e$N[-1, 2, 2] > 0))
  expect_true(all(y6e$N[-1, 2, 5] > 0))
})
