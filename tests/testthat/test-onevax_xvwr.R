context("onevax_xvwr")

test_that("run_onevax_xvwr works correctly", {
  tt <- seq(0, 5)
  gp <- gono_params(1:2)
  y1 <- run_onevax_xvwr(tt, gp, vea = 0, dur = 1e3)[[1]]

  # check no-one is vaccinated with v switched off
  expect_true(all(y1$cum_vaccinated == 0))
  y2 <- run_onevax_xvwr(tt, gp, vea = 0, dur = 1e3, vbe = 1)
  # check 100% vbe vaccinates all new entrants
  expect_equal(diff(rowSums(y2[[1]]$cum_vaccinated[, , 1])), rep(12e3, max(tt)))
  # and no-one else
  expect_equal(sum(y2[[1]]$cum_vaccinated[, , 2:3]), 0)

  # check can restart
  init_params <- lapply(y2, restart_params)
  y3 <- run_onevax_xvwr(seq(max(tt), length.out = 2, by = 1),
                        gp, init_params,
                        vea = 0, dur = 1e3, vbe = 1)
  for (i in seq_along(y2)) {
    expect_equal(y2[[i]]$U[length(tt), , ], y3[[i]]$U[1, , ])
    expect_equal(y2[[i]]$I[length(tt), , ], y3[[i]]$I[1, , ])
    expect_equal(y2[[i]]$A[length(tt), , ], y3[[i]]$A[1, , ])
    expect_equal(y2[[i]]$S[length(tt), , ], y3[[i]]$S[1, , ])
    expect_equal(y2[[i]]$T[length(tt), , ], y3[[i]]$T[1, , ])
  }

  primary_uptake <- booster_uptake <- 0.5
  # check VoD is working correctly
  y3e <- run_onevax_xvwr(tt, gp, vea = 1, dur = 1e3, strategy = "VoD",
                         primary_uptake = primary_uptake,
                         booster_uptake = booster_uptake)

  for (i in seq_along(y3e)) {
    # no-one in stratum V is vaccinated again
    expect_equal(sum(y3e[[i]]$cum_vaccinated[, , 2]), 0)
    # uptake % of treated are vaccinated
    expect_equal(y3e[[i]]$cum_vaccinated[, , -2] / primary_uptake,
                 y3e[[i]]$cum_treated[, , -2])
    # efficacy is perfect
    expect_equal(sum(y3e[[i]]$cum_treated[, , 2]), 0)
    # no-one is lost
    expect_equal(apply(y3e[[i]]$N, 1, sum), rep(6e5, 6), tolerance = 1e-5)
  }


  # check VoA is working correctly
  y4e <- run_onevax_xvwr(tt, gp, vea = 1, dur = 1e3, strategy = "VoA",
                         primary_uptake = primary_uptake,
                         booster_uptake = booster_uptake)

  for (i in seq_along(y4e)) {
    # no-one in stratum V is vaccinated again
  expect_equal(sum(y4e[[i]]$cum_vaccinated[, , 2]), 0)
  expect_equal(y4e[[i]]$cum_vaccinated[, , c(1, 3)] / primary_uptake,
               y4e[[i]]$cum_screened[, , c(1, 3)] +
                 y4e[[i]]$cum_treated[, , c(1, 3)])
  expect_equal(apply(y4e[[i]]$N, 1, sum), rep(6e5, 6), tolerance = 1e-5)
  }
  # check vaccination targeting
  y5e <- run_onevax_xvwr(tt, gp, vea = 1, dur = 1e3, strategy = "VoD(L)+VoA(H)",
                         primary_uptake = primary_uptake,
                         booster_uptake = booster_uptake)

  for (i in seq_along(y5e)) {
    # no-one in stratum V is vaccinated again
    expect_equal(sum(y5e[[i]]$cum_vaccinated[, , 2]), 0)
    # only treated L are vaccinated
    expect_equal(y5e[[i]]$cum_vaccinated[, 1, c(1, 3)] / primary_uptake,
                 y5e[[i]]$cum_treated[, 1, c(1, 3)])
    # all attending H are vaccinated
    expect_equal(y5e[[i]]$cum_vaccinated[, 2, c(1, 3)] / primary_uptake,
                 y5e[[i]]$cum_treated[, 2, c(1, 3)] +
                   y5e[[i]]$cum_screened[, 2, c(1, 3)])
  }

  # check length of uptake vector must be 1 or length(gp)
  expect_error(run_onevax_xvwr(tt, gp, vea = 1, dur = 1e3, strategy = "VbE",
                  primary_uptake = c(0, 0.5, 1)))
  # check length of vea must be 1
  expect_error(run_onevax_xvwr(tt, gp, vea = c(0, 1, 1), dur = 1e3,
                               strategy = "VbE", primary_uptake = 1))
  expect_error(run_onevax_xvwr(tt, gp, vea = 1, dur = c(0, 1e2, 1e3),
                               strategy = "VbE", primary_uptake = 1))

  ## test revax is working

  y6e <- run_onevax_xvwr(tt, gp, vea = 1, dur = 1e-10, dur_revax = 1e10,
                         strategy = "VoD(L)+VoA(H)",
                         primary_uptake = primary_uptake,
                         booster_uptake = booster_uptake)

  for (i in seq_along(y6e)) {
    # no-one in stratum V is vaccinated again
    expect_equal(sum(y6e[[i]]$cum_vaccinated[, , 2]), 0)
    # only treated L are vaccinated
    expect_equal(y6e[[i]]$cum_vaccinated[, 1, c(1, 3)] / primary_uptake,
                 y6e[[i]]$cum_treated[, 1, c(1, 3)])
    # all attending H are vaccinated
    expect_equal(y6e[[i]]$cum_vaccinated[, 2, c(1, 3)] / primary_uptake,
                 y6e[[i]]$cum_treated[, 2, c(1, 3)] +
                   y6e[[i]]$cum_screened[, 2, c(1, 3)])
    # stratum V empties immediately
    expect_equal(sum(y6e[[i]]$N[, , 2]), 0, tolerance = 1e-5)
    expect_true(all(y6e[[i]]$cum_vaccinated[, , 3] <=
                      y6e[[i]]$cum_vaccinated[, , 1]))
    # efficacy is perfect
    expect_equal(sum(y6e[[i]]$cum_incid[, , c(2, 4)]), 0)
    expect_true(all(y6e[[i]]$cum_incid[-1, , c(1, 3)] > 0))
  }

  y7e <- run_onevax_xvwr(tt, gp, vea = 0, vea_revax = 1, dur = 1,
                         strategy = "VoD(L)+VoA(H)",
                         primary_uptake = primary_uptake,
                         booster_uptake = booster_uptake)

  for (i in seq_along(y7e)) {
    # no-one in stratum V is vaccinated again
    expect_equal(sum(y7e[[i]]$cum_vaccinated[, , 2]), 0)
    # only treated L are vaccinated
    expect_equal(y7e[[i]]$cum_vaccinated[, 1, c(1, 3)] / primary_uptake,
                 y7e[[i]]$cum_treated[, 1, c(1, 3)])
    # all attending H are vaccinated
    expect_equal(y7e[[i]]$cum_vaccinated[, 2, c(1, 3)] / primary_uptake,
                 y7e[[i]]$cum_treated[, 2, c(1, 3)] +
                   y7e[[i]]$cum_screened[, 2, c(1, 3)])
    # efficacy is perfect in R
    expect_equal(sum(y7e[[i]]$cum_incid[, , 4]), 0)
    # but not in XWV
    expect_true(all(y7e[[i]]$cum_incid[-1, , -4] > 0))
  }
})
