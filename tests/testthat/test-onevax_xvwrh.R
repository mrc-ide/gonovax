#testing xvwrh model

test_that("run_onevax_xvwrh works correctly", {
    tt <- seq(0, 5)
    gp <- gonovax::gono_params(1:2)
    y1 <- run_onevax_xvwrh(tt, gp, vea = 0, dur = 1e3)[[1]]

    # check no-one is vaccinated with v switched off
    expect_true(all(y1$cum_vaccinated == 0))

    # check 100% vbe vaccinates all new entrants
    y2 <- run_onevax_xvwrh(tt, gp, vea = 0, dur = 1e3, vbe = 1)

    # and no-one else
    expect_equal(sum(y2[[1]]$cum_vaccinated[, , 2:5]), 0)

    expect_equal(diff(rowSums(y2[[1]]$cum_vaccinated[, , 1])), rep(12e3,
                                                                   max(tt)))
    # check this is still the case when hesitancy > 0
    y2.1 <- run_onevax_xvwrh(tt, gp, vea = 0, dur = 1e3, vbe = 1, hes = 0.5)

    expect_equal(diff(rowSums(y2.1[[1]]$cum_vaccinated[, , 1])), rep((12e3) / 2,
                                                                   max(tt)))
    # check no vaccination in hesitant entrants for vbe = 100%

    expect_equal((rowSums(y2.1[[1]]$cum_vaccinated[, , 5])), rep(0,
                                                                 max(tt) + 1))

    # check can restart
    init_params <- lapply(y2, restart_params)
    y3 <- run_onevax_xvwrh(seq(max(tt), length.out = 2, by = 1),
                          gp, init_params,
                          vea = 0, dur = 1e3, vbe = 1)
    for (i in seq_along(y2)) {
      expect_equal(y2[[i]]$U[length(tt), , ], y3[[i]]$U[1, , ])
      expect_equal(y2[[i]]$I[length(tt), , ], y3[[i]]$I[1, , ])
      expect_equal(y2[[i]]$A[length(tt), , ], y3[[i]]$A[1, , ])
      expect_equal(y2[[i]]$S[length(tt), , ], y3[[i]]$S[1, , ])
      expect_equal(y2[[i]]$T[length(tt), , ], y3[[i]]$T[1, , ])
    }

    y_h <- run_onevax_xvwrh(tt, gp, vea = 0, dur = 1e3, vbe = 1, hes = 0.3)

    # initial population split between non-vaccinated and hesitant only
    # other stratum empty

    for (i in seq_along(y_h)) {
    expect_true(all(y_h[[i]]$N[1, , 5] > 0))
    expect_true(all(y_h[[i]]$N[1, , 1] > 0))
    expect_true(all(y_h[[i]]$N[1, , 2] == 0))
    expect_true(all(y_h[[i]]$N[1, , 3] == 0))
    expect_true(all(y_h[[i]]$N[1, , 4] == 0))
    }

    # incidence in hesitant population increasing

    for (i in seq_along(y_h)) {
      expect_true(all(y_h[[i]]$cum_incid[-1, , 5] > 0))
    }

    # yearly population entrants enter X and H strata
    # in accordance with assigned proportion of hesitancy 'hes'

      # H and X stratum sizes remain constant in time

    y_h2 <- run_onevax_xvwrh(tt, gp, vea = 0, dur = 1e3, hes = 0.3)

    expect_equal(y_h2[[1]]$N[1, 1, 1], y_h2[[1]]$N[5, 1, 1])
    expect_equal(y_h2[[1]]$N[1, 2, 1], y_h2[[1]]$N[5, 2, 1])
    expect_equal(y_h2[[1]]$N[1, 1, 5], y_h2[[1]]$N[5, 1, 5])
    expect_equal(y_h2[[1]]$N[1, 2, 5], y_h2[[1]]$N[5, 2, 5])

      # H and X stratum equal when no vaccination and hes = 0.5

    y_h3 <- run_onevax_xvwrh(tt, gp, vea = 0, dur = 1e3, hes = 0.5)

    expect_equal(y_h3[[1]]$N[, , 1], y_h3[[1]]$N[, , 5])

      # Number of infections in X and H equal for no vaccination and hes = 0.5

    expect_equal(y_h3[[1]]$cum_incid[, , 1], y_h3[[1]]$cum_incid[, , 5])

    # if proportion hesitant is 0%, = outputs same as xvwr model

    y_h4 <- run_onevax_xvwrh(tt, gp, vea = 0, dur = 1e3, vbe = 1, hes = 0)
    y_xvwr <- run_onevax_xvwr(tt, gp, vea = 0, dur = 1e3, vbe = 1)

    expect_equal(rowSums(y_h4[[1]]$N[, , 1]), rowSums(y_xvwr[[1]]$N[, , 1]))
    expect_equal(rowSums(y_h4[[1]]$N[, , 2]), rowSums(y_xvwr[[1]]$N[, , 2]))
    expect_equal(rowSums(y_h4[[1]]$N[, , 3]), rowSums(y_xvwr[[1]]$N[, , 3]))

    expect_equal(rowSums(y_h4[[1]]$U[, , 1]), rowSums(y_xvwr[[1]]$U[, , 1]))
    expect_equal(rowSums(y_h4[[1]]$U[, , 2]), rowSums(y_xvwr[[1]]$U[, , 2]))
    expect_equal(rowSums(y_h4[[1]]$U[, , 3]), rowSums(y_xvwr[[1]]$U[, , 3]))

    uptake <- c(0.5, 1)

    # check VoD is working correctly
    y3e <- run_onevax_xvwrh(tt, gp, vea = 1, dur = 1e3, strategy = "VoD",
                           uptake = uptake)

    for (i in seq_along(y3e)) {
      # no-one in stratum V is vaccinated again
      expect_equal(sum(y3e[[i]]$cum_vaccinated[, , 2]), 0)
      # uptake % of treated are vaccinated
      expect_equal(y3e[[i]]$cum_vaccinated[, , -2] / uptake[i],
                   y3e[[i]]$cum_treated[, , -2])
      # efficacy is perfect
      expect_equal(sum(y3e[[i]]$cum_treated[, , 2]), 0)
      # no-one is lost
      expect_equal(apply(y3e[[i]]$N, 1, sum), rep(6e5, 6), tolerance = 1e-5)
    }


    # check VoA is working correctly
    y4e <- run_onevax_xvwrh(tt, gp, vea = 1, dur = 1e3, strategy = "VoA",
                           uptake = uptake)
    for (i in seq_along(y4e)) {
      # no-one in stratum V is vaccinated again
      expect_equal(sum(y4e[[i]]$cum_vaccinated[, , 2]), 0)
      expect_equal(y4e[[i]]$cum_vaccinated[, , c(1, 3)] / uptake[i],
                   y4e[[i]]$cum_screened[, , c(1, 3)] +
                     y4e[[i]]$cum_treated[, , c(1, 3)])
      expect_equal(apply(y4e[[i]]$N, 1, sum), rep(6e5, 6), tolerance = 1e-5)
    }
    # check vaccination targeting
    y5e <- run_onevax_xvwrh(tt, gp, vea = 1, dur = 1e3,
                            strategy = "VoD(L)+VoA(H)",
                           uptake = uptake)
    for (i in seq_along(y5e)) {
      # no-one in stratum V is vaccinated again
      expect_equal(sum(y5e[[i]]$cum_vaccinated[, , 2]), 0)
      # only treated L are vaccinated
      expect_equal(y5e[[i]]$cum_vaccinated[, 1, c(1, 3)] / uptake[i],
                   y5e[[i]]$cum_treated[, 1, c(1, 3)])
      # all attending H are vaccinated
      expect_equal(y5e[[i]]$cum_vaccinated[, 2, c(1, 3)] / uptake[i],
                   y5e[[i]]$cum_treated[, 2, c(1, 3)] +
                     y5e[[i]]$cum_screened[, 2, c(1, 3)])
    }

    # check length of uptake vector must be 1 or length(gp)
    expect_error(run_onevax_xvwrh(tt, gp, vea = 1, dur = 1e3, strategy = "VbE",
                                 uptake = c(0, 0.5, 1)))
    # check length of vea must be 1
    expect_error(run_onevax_xvwrh(tt, gp, vea = c(0, 1, 1), dur = 1e3,
                                 strategy = "VbE", uptake = 1))
    expect_error(run_onevax_xvwrh(tt, gp, vea = 1, dur = c(0, 1e2, 1e3),
                                 strategy = "VbE", uptake = 1))

    ## test revax is working

    y6e <- run_onevax_xvwrh(tt, gp, vea = 1, dur = 1e-10, dur_revax = 1e10,
                           strategy = "VoD(L)+VoA(H)",
                           uptake = uptake)
    for (i in seq_along(y6e)) {
      # no-one in stratum V is vaccinated again
      expect_equal(sum(y6e[[i]]$cum_vaccinated[, , 2]), 0)
      # only treated L are vaccinated
      expect_equal(y6e[[i]]$cum_vaccinated[, 1, c(1, 3)] / uptake[i],
                   y6e[[i]]$cum_treated[, 1, c(1, 3)])
      # all attending H are vaccinated
      expect_equal(y6e[[i]]$cum_vaccinated[, 2, c(1, 3)] / uptake[i],
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

    y7e <- run_onevax_xvwrh(tt, gp, vea = 0, vea_revax = 1, dur = 1,
                           strategy = "VoD(L)+VoA(H)",
                           uptake = uptake, hes = 0.3)
    for (i in seq_along(y7e)) {
      # no-one in stratum V is vaccinated again
      expect_equal(sum(y7e[[i]]$cum_vaccinated[, , 2]), 0)
      # only treated L are vaccinated
      expect_equal(y7e[[i]]$cum_vaccinated[, 1, c(1, 3)] / uptake[i],
                   y7e[[i]]$cum_treated[, 1, c(1, 3)])
      # all attending H are vaccinated
      expect_equal(y7e[[i]]$cum_vaccinated[, 2, c(1, 3)] / uptake[i],
                   y7e[[i]]$cum_treated[, 2, c(1, 3)] +
                     y7e[[i]]$cum_screened[, 2, c(1, 3)])
      # efficacy is perfect in R
      expect_equal(sum(y7e[[i]]$cum_incid[, , 4]), 0)
      # but not in XWV
      expect_true(all(y7e[[i]]$cum_incid[-1, , -c(4)] > 0))
    }

  ## test restart with hesitancy is working

  y8 <- run_onevax_xvwrh(tt, gp, vea = 0, dur = 1e3)

  i_p <- lapply(y8, restart_hes, n_vax = 5, hes = 0.5)
  y_hesres <- run_onevax_xvwrh(tt, gp, init_params = i_p, vea = 0, dur = 1e3,
                               hes = 0.5)

  # final timepoint y8 run = 2 * first timepoint of y_hesres run (as hes = 0.5)
  # (for XVWR)
  for (i in seq_along(y8)) {
    expect_equal(y8[[i]]$U[length(tt), , 1:4], y_hesres[[i]]$U[1, , 1:4] * 2)
    expect_equal(y8[[i]]$I[length(tt), , 1:4], y_hesres[[i]]$I[1, , 1:4] * 2)
    expect_equal(y8[[i]]$A[length(tt), , 1:4], y_hesres[[i]]$A[1, , 1:4] * 2)
    expect_equal(y8[[i]]$S[length(tt), , 1:4], y_hesres[[i]]$S[1, , 1:4] * 2)
    expect_equal(y8[[i]]$T[length(tt), , 1:4], y_hesres[[i]]$T[1, , 1:4] * 2)
  }

  # restart_hes moves X -> H only, and correctly

  for (i in seq_along(y_hesres)) {
    expect_equal(y_hesres[[i]]$U[, , 1], y_hesres[[i]]$U[, , 5])
    expect_equal(y_hesres[[i]]$I[, , 1], y_hesres[[i]]$I[, , 5])
    expect_equal(y_hesres[[i]]$A[, , 1], y_hesres[[i]]$A[, , 5])
    expect_equal(y_hesres[[i]]$S[, , 1], y_hesres[[i]]$S[, , 5])
    expect_equal(y_hesres[[i]]$T[, , 1], y_hesres[[i]]$T[, , 5])

    expect_equal(rowSums(y8[[i]]$U[, , 2:4]), rep(0, length(tt)))
    expect_equal(rowSums(y8[[i]]$I[, , 2:4]), rep(0, length(tt)))
    expect_equal(rowSums(y8[[i]]$A[, , 2:4]), rep(0, length(tt)))
    expect_equal(rowSums(y8[[i]]$S[, , 2:4]), rep(0, length(tt)))
    expect_equal(rowSums(y8[[i]]$T[, , 2:4]), rep(0, length(tt)))

  }


})
