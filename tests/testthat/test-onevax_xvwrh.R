#testing xvwrh model

test_that("run_onevax_xvwrh works correctly", {
  tt <- seq(0, 5)
  gp <- gono_params(1:2)
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

  # expect error if inputs are not of length 1 or equal to length of params

  uptake <- c(0, 2.5, 0.5, 0.75, 1)
  expect_error(y <- run_onevax_xvwrh(tt, gp, vea = 0, dur = 1e3, vbe = 1,
                                     primary_uptake = uptake))

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
  init_pop <- 600000 * c(0.85, 0.15)
  for (i in seq_along(y_h)) {
    expect_equivalent(y_h[[i]]$N[1, , 1], init_pop * 0.7)
    expect_equivalent(y_h[[i]]$N[1, , 5], init_pop * 0.3)
    expect_true(all(y_h[[i]]$N[1, , 2:4] == 0))
  }

  # incidence in hesitant population increasing
  for (i in seq_along(y_h)) {
    expect_true(all(y_h[[i]]$cum_incid[-1, , 5] > 0))
  }

  y_h2 <- run_onevax_xvwrh(tt, gp, vea = 0, dur = 1e3, vbe = 0, hes = 0.3)
  # yearly population entrants enter X and H strata
  # in accordance with assigned proportion of hesitancy 'hes'
  # H and X stratum sizes remain constant in time
  expect_equal(y_h2[[1]]$N[1, , ], y_h2[[1]]$N[5, , ])

  # H and X stratum equal when no vaccination and hes = 0.5
  y_h3 <- run_onevax_xvwrh(tt, gp, vea = 0, dur = 1e3, hes = 0.5)
  expect_equal(y_h3[[1]]$N[, , 1], y_h3[[1]]$N[, , 5])

  # Number of infections in X and H equal for no vaccination and hes = 0.5
  expect_equal(y_h3[[1]]$cum_incid[, , 1], y_h3[[1]]$cum_incid[, , 5])

  # if proportion hesitant is 0%, = outputs same as xvwr model
  # choose a difficult case where there are very few zero outputs.
  y_h4 <- run_onevax_xvwrh(tt, gp, vea = 0.5, dur = 1, vbe = 0.8, hes = 0,
                           primary_uptake = 0.5, booster_uptake = 0.3,
                           strategy = "VoD(L)+VoA(H)")
  y_xvwr <- run_onevax_xvwr(tt, gp, vea = 0.5, dur = 1, vbe = 0.8,
                            primary_uptake = 0.5, booster_uptake = 0.3,
                            strategy = "VoD(L)+VoA(H)")

  expect_equal(y_h4[[1]]$N[, , 1:4], y_xvwr[[1]]$N)
  expect_equal(y_h4[[1]]$U[, , 1:4], y_xvwr[[1]]$U)
  expect_equal(y_h4[[1]]$cum_incid[, , 1:4], y_xvwr[[1]]$cum_incid)

  primary_uptake <- c(0.25, 0.5)
  booster_uptake <- c(0.3, 0.6)

  # check VoD is working correctly
  y3e <- run_onevax_xvwrh(tt, gp, vea = 0.5, dur = 1, strategy = "VoD",
                          primary_uptake = primary_uptake,
                          booster_uptake = booster_uptake)

  for (i in seq_along(y3e)) {

    ## all treated in X or W are offered vaccination
    expect_equal(y3e[[i]]$cum_offered[, , c(1, 3)],
                 y3e[[i]]$cum_treated[, , c(1, 3)])
    # and no-one else
    expect_equal(sum(y3e[[i]]$cum_offered[, , -c(1, 3)]), 0)

    # uptake % of offered are vaccinated
    expect_equal(y3e[[i]]$cum_offered[, , 1] * primary_uptake[i],
                 y3e[[i]]$cum_vaccinated[, , 1])
    expect_equal(y3e[[i]]$cum_offered[, , 3] * booster_uptake[i],
                 y3e[[i]]$cum_vaccinated[, , 3])
    # and no-one else
    expect_equal(sum(y3e[[i]]$cum_vaccinated[, , -c(1, 3)]), 0)

    # no-one is lost
    expect_equal(apply(y3e[[i]]$N, 1, sum), rep(6e5, 6), tolerance = 1e-5)
  }


  # check VoA is working correctly
  y4e <- run_onevax_xvwrh(tt, gp, vea = 0.5, dur = 1, strategy = "VoA",
                          primary_uptake = primary_uptake,
                          booster_uptake = booster_uptake)

  for (i in seq_along(y4e)) {

    ## all treated or screened in X or W are offered vaccination
    expect_equal(y4e[[i]]$cum_offered[, , c(1, 3)],
                 y4e[[i]]$cum_treated[, , c(1, 3)] +
                   y4e[[i]]$cum_screened[, , c(1, 3)])
    # and no-one else
    expect_equal(sum(y4e[[i]]$cum_offered[, , -c(1, 3)]), 0)

    # uptake % of offered are vaccinated
    expect_equal(y4e[[i]]$cum_offered[, , 1] * primary_uptake[i],
                 y4e[[i]]$cum_vaccinated[, , 1])
    expect_equal(y4e[[i]]$cum_offered[, , 3] * booster_uptake[i],
                 y4e[[i]]$cum_vaccinated[, , 3])
    # and no-one else
    expect_equal(sum(y4e[[i]]$cum_vaccinated[, , -c(1, 3)]), 0)

    # no-one is lost
    expect_equal(apply(y4e[[i]]$N, 1, sum), rep(6e5, 6), tolerance = 1e-5)
  }
  # check vaccination targeting
  y5e <- run_onevax_xvwrh(tt, gp, vea = 0.5, dur = 1,
                          strategy = "VoD(L)+VoA(H)",
                          primary_uptake = primary_uptake,
                          booster_uptake = booster_uptake)

  for (i in seq_along(y5e)) {
    ## L who are treated in X or W are offered vaccination
    expect_equal(y5e[[i]]$cum_offered[, 1, c(1, 3)],
                 y5e[[i]]$cum_treated[, 1, c(1, 3)])
    ## H who are treated or screened in X or W are offered vaccination
    expect_equal(y5e[[i]]$cum_offered[, 2, c(1, 3)],
                 y5e[[i]]$cum_treated[, 2,  c(1, 3)] +
                   y5e[[i]]$cum_screened[, 2,  c(1, 3)])
    # and no-one else
    expect_equal(sum(y5e[[i]]$cum_offered[, , -c(1, 3)]), 0)

    # uptake % of offered are vaccinated
    expect_equal(y5e[[i]]$cum_offered[, , 1] * primary_uptake[i],
                 y5e[[i]]$cum_vaccinated[, , 1])
    expect_equal(y5e[[i]]$cum_offered[, , 3] * booster_uptake[i],
                 y5e[[i]]$cum_vaccinated[, , 3])
    # and no-one else
    expect_equal(sum(y5e[[i]]$cum_vaccinated[, , -c(1, 3)]), 0)

    # no-one is lost
    expect_equal(apply(y5e[[i]]$N, 1, sum), rep(6e5, 6), tolerance = 1e-5)
  }

  # check length of uptake vector must be 1 or length(gp)
  expect_error(run_onevax_xvwrh(tt, gp, vea = 1, dur = 1e3, strategy = "VbE",
                                primary_uptake = c(0, 0.5, 1)))
  # check length of vea must be 1
  expect_error(run_onevax_xvwrh(tt, gp, vea = c(0, 1, 1), dur = 1e3,
                                strategy = "VbE", primary_uptake = 1))
  # check length of dur must be 1
  expect_error(run_onevax_xvwrh(tt, gp, vea = 1, dur = c(0, 1e2, 1e3),
                                strategy = "VbE", primary_uptake = 1))

  ## test revax is working
  y6e <- run_onevax_xvwrh(tt, gp, vea = 1, dur = 1e-10, dur_revax = 1e10,
                          strategy = "VoD(L)+VoA(H)",
                          primary_uptake = primary_uptake,
                          booster_uptake = booster_uptake)
  for (i in seq_along(y6e)) {
    ## L who are treated in X or W are offered vaccination
    expect_equal(y6e[[i]]$cum_offered[, 1, c(1, 3)],
                 y6e[[i]]$cum_treated[, 1, c(1, 3)])
    ## H who are treated or screened in X or W are offered vaccination
    expect_equal(y6e[[i]]$cum_offered[, 2, c(1, 3)],
                 y6e[[i]]$cum_treated[, 2,  c(1, 3)] +
                   y6e[[i]]$cum_screened[, 2,  c(1, 3)])
    # and no-one else
    expect_equal(sum(y6e[[i]]$cum_offered[, , -c(1, 3)]), 0)

    # uptake % of offered are vaccinated
    expect_equal(y6e[[i]]$cum_offered[, , 1] * primary_uptake[i],
                 y6e[[i]]$cum_vaccinated[, , 1])
    expect_equal(y6e[[i]]$cum_offered[, , 3] * booster_uptake[i],
                 y6e[[i]]$cum_vaccinated[, , 3])

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
                          primary_uptake = primary_uptake,
                          booster_uptake = booster_uptake,
                          hes = 0.3)

  for (i in seq_along(y7e)) {
    ## L who are treated in X or W are offered vaccination
    expect_equal(y7e[[i]]$cum_offered[, 1, c(1, 3)],
                 y7e[[i]]$cum_treated[, 1, c(1, 3)])
    ## H who are treated or screened in X or W are offered vaccination
    expect_equal(y7e[[i]]$cum_offered[, 2, c(1, 3)],
                 y7e[[i]]$cum_treated[, 2,  c(1, 3)] +
                   y7e[[i]]$cum_screened[, 2,  c(1, 3)])
    # and no-one else
    expect_equal(sum(y7e[[i]]$cum_offered[, , -c(1, 3)]), 0)

    # uptake % of offered are vaccinated
    expect_equal(y7e[[i]]$cum_offered[, , 1] * primary_uptake[i],
                 y7e[[i]]$cum_vaccinated[, , 1])
    expect_equal(y7e[[i]]$cum_offered[, , 3] * booster_uptake[i],
                 y7e[[i]]$cum_vaccinated[, , 3])
    # efficacy is perfect in R
    expect_equal(sum(y7e[[i]]$cum_incid[, , 4]), 0)
    # but not in XWVH
    expect_true(all(y7e[[i]]$cum_incid[-1, , -c(4)] > 0))
  }

  ## test restart with hesitancy is working

  y8 <- run_onevax_xvwrh(tt, gp, vea = 0, dur = 1e3)

  i_p <- lapply(y8, restart_hes, hes = 0.5)
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

  # restart_hes gives error if baseline model run provided already has hesitancy

  y9 <- run_onevax_xvwrh(tt, gp, vea = 0, dur = 1e3, hes = 0.5)

  expect_error(lapply(y9, restart_hes, n_vax = 5, hes = 0.5))

  # restart_hes gives error if baseline model run provided contains vaccinated

  y10 <- run_onevax_xvwrh(tt, gp, vea = 0, dur = 1e3, vbe = 1)

  expect_error(lapply(y10, restart_hes, n_vax = 5, hes = 0.5))

  # check booster_uptake defaults to primary_uptake
  y_prim_only <- run_onevax_xvwrh(tt, gp, vea = 0, dur = 1,
                                  primary_uptake = 0.75,
                                  strategy = "VoD(L)+VoA(H)")
  y_prim_boost <- run_onevax_xvwrh(tt, gp, vea = 0, dur = 1,
                                   primary_uptake = 0.75, booster_uptake = 0.75,
                                   strategy = "VoD(L)+VoA(H)")

  expect_equal(y_prim_only, y_prim_boost)
})
