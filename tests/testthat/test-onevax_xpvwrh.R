#testing xpvwrh model

test_that("run_onevax_xpvwrh works correctly", {
  tt <- seq(0, 5)
  gp <- gono_params(1:2)
  y1 <- run_onevax_xpvwrh(tt, gp, vea = 0, dur_v = 1e3)[[1]]

  # check no-one is vaccinated with v switched off
  expect_true(all(y1$cum_vaccinated == 0))

  # check 100% vbe vaccinates all new entrants
  y2 <- run_onevax_xpvwrh(tt, gp, vea = 0, dur_v = 1e10, vbe = 1)

  # cum_vaccinated = 12000 each year = number of entrants
  expect_equal(diff(rowSums(y2[[1]]$cum_vaccinated[, , 1])), rep(12e3,
                                                                 max(tt)))

  # this also means there should be people moving into 'V'
  expect_true(sum(y2[[1]]$N[, , 3]) > 0)

  # but not P
  expect_equal(sum(y2[[1]]$N[, , 2]), 0)

  # other compartments empty
  expect_equal(sum(y2[[1]]$N[, , -c(1, 3, 4)]), 0)

  # and no-one else
  expect_equal(sum(y2[[1]]$cum_vaccinated[, , 2:5]), 0)
  expect_equal(sum(y2[[1]]$N[, , -c(1, 3, 4)]), 0)

  # check population is constant in size over time

  for (i in seq_along(y2))  {
  expect_equal(rowSums(y2[[i]]$N[, , ]), rep(6e+05, max(tt) + 1))
  }

  # check this is still the case when hesitancy > 0
  y2.1 <- run_onevax_xpvwrh(tt, gp, vea = 0, dur_v = 1e3, vbe = 1, hes = 0.5)

  for (i in seq_along(y2.1)) {
  expect_equal(diff(rowSums(y2.1[[i]]$cum_vaccinated[, , 1])), rep((12e3) / 2,
                                                                   max(tt)))
  expect_true(sum(y2[[i]]$N[, , 3]) > 0)
  # check no vaccination in hesitant entrants for vbe = 100%

  expect_equal((rowSums(y2.1[[i]]$cum_vaccinated[, , 5])), rep(0,
                                                               max(tt) + 1))
}
  # expect error if inputs are not of length 1 or equal to length of params

  uptake <- c(0, 2.5, 0.5, 0.75, 1)
  expect_error(y <- run_onevax_xpvwrh(tt, gp, vea = 0, dur_v = 1e3, vbe = 1,
                                     primary_uptake = uptake))

  # check can restart
  init_params <- lapply(y2, restart_params)
  y3 <- run_onevax_xpvwrh(seq(max(tt), length.out = 2, by = 1),
                         gp, init_params,
                         vea = 0, dur_v = 1e3, vbe = 1)
  for (i in seq_along(y2)) {
    expect_equal(y2[[i]]$U[length(tt), , ], y3[[i]]$U[1, , ])
    expect_equal(y2[[i]]$I[length(tt), , ], y3[[i]]$I[1, , ])
    expect_equal(y2[[i]]$A[length(tt), , ], y3[[i]]$A[1, , ])
    expect_equal(y2[[i]]$S[length(tt), , ], y3[[i]]$S[1, , ])
    expect_equal(y2[[i]]$T[length(tt), , ], y3[[i]]$T[1, , ])
  }

  y_h <- run_onevax_xpvwrh(tt, gp, vea = 0, dur_v = 1e3, vbe = 1, hes = 0.3)

  # initial population split between non-vaccinated and hesitant only
  # other stratum empty
  init_pop <- 600000 * c(0.85, 0.15)
  for (i in seq_along(y_h)) {
    expect_equivalent(y_h[[i]]$N[1, , 1], init_pop * 0.7)
    expect_equivalent(y_h[[i]]$N[1, , 6], init_pop * 0.3)
    expect_true(all(y_h[[i]]$N[1, , 2:5] == 0))
  }

  # incidence in hesitant population increasing
  for (i in seq_along(y_h)) {
    expect_true(all(y_h[[i]]$cum_incid[-1, , 6] > 0))
  }

  y_h2 <- run_onevax_xpvwrh(tt, gp, vea = 0, dur_v = 1e3, vbe = 0, hes = 0.3)

  for (i in seq_along(y_h2)) {
  # yearly population entrants enter X and H strata
  # in accordance with assigned proportion of hesitancy 'hes'
  # H and X stratum sizes remain constant in time
  expect_equal(y_h2[[i]]$N[1, , ], y_h2[[i]]$N[6, , ])
}
  # H and X stratum equal when no vaccination and hes = 0.5
  y_h3 <- run_onevax_xpvwrh(tt, gp, vea = 0, dur_v = 1e3, hes = 0.5)

  for (i in seq_along(y_h3)) {
  expect_equal(y_h3[[i]]$N[, , 1], y_h3[[i]]$N[, , 6])

  # Number of infections in X and H equal for no vaccination and hes = 0.5
  expect_equal(y_h3[[i]]$cum_incid[, , 1], y_h3[[i]]$cum_incid[, , 6])
}
  # if proportion hesitant is 0%, = outputs same as xvwr model
  # choose a difficult case where there are very few zero outputs.
  y_h4 <- run_onevax_xpvwrh(tt, gp, vea = 0.5, dur_v = 1, vbe = 0.8, hes = 0,
                           r2 = 1, r1 = 0.5, booster_uptake = 0.3,
                           strategy = "VoD(L)+VoA(H)")
  y_xvwr <- run_onevax_xvwr(tt, gp, vea = 0.5, dur = 1, vbe = 0.8,
                            primary_uptake = 0.5, booster_uptake = 0.3,
                            strategy = "VoD(L)+VoA(H)")

  for (i in seq_along(y_h4)) {
  expect_equal(y_h4[[i]]$N[, , c(1, 3:5)], y_xvwr[[i]]$N)
  expect_equal(y_h4[[i]]$U[, , c(1, 3:5)], y_xvwr[[i]]$U)
  expect_equal(y_h4[[i]]$cum_incid[, , c(1, 3:5)], y_xvwr[[i]]$cum_incid)

}

  # Test the vaccination maps are being generated as expected

  pars <- lapply(gp[1], model_params)
  vbe <- 1
  p <- set_strategy(strategy = "VoD(L)+VoA(H)", vbe > 0)
  i_eligible <- c(1, 1, 2, 4)
  i_p <- c(2, 3, 3, 5)

  vod_map <- create_vax_map_branching(n_vax = 6, p$vod, i_eligible, i_p)
  vos_map <- create_vax_map_branching(n_vax = 6, p$vos, i_eligible, i_p)
  vbe_map <- create_vax_map_branching(n_vax = 6, p$vbe, i_eligible, i_p,
                                      set_vbe = TRUE)

  # for vod, expect:
  expect_true(unique(vod_map[, 1, 1] == c(1, 1)))
  expect_true(unique(vod_map[, 2, 1] == c(-1, -1)))
  expect_true(unique(vod_map[, 3, 1] == c(-1, -1)))
  expect_true(unique(vod_map[, 2, 2] == c(1, 1)))
  expect_true(unique(vod_map[, 3, 2] == c(-1, -1)))
  expect_true(unique(vod_map[, 4, 4] == c(1, 1)))
  expect_true(unique(vod_map[, 5, 4] == c(-1, -1)))

  expect_equal(sum(vod_map[, -c(1, 2, 3), 1]), 0)
  expect_equal(sum(vod_map[, -c(2, 3), 2]), 0)
  expect_equal(sum(vod_map[, -c(4, 5), 4]), 0)

  expect_equal(sum(vod_map[, , c(3, 6, 5)]), 0)

  # for vos, expect:
  expect_true(unique(vos_map[, 1, 1] == c(0, 1)))
  expect_true(unique(vos_map[, 2, 1] == c(0, -1)))
  expect_true(unique(vos_map[, 3, 1] == c(0, -1)))
  expect_true(unique(vos_map[, 2, 2] == c(0, 1)))
  expect_true(unique(vos_map[, 3, 2] == c(0, -1)))
  expect_true(unique(vos_map[, 4, 4] == c(0, 1)))
  expect_true(unique(vos_map[, 5, 4] == c(0, -1)))

  expect_equal(sum(vos_map[, -c(1, 2, 3), 1]), 0)
  expect_equal(sum(vos_map[, -c(2, 3), 2]), 0)
  expect_equal(sum(vos_map[, -c(4, 5), 4]), 0)

  expect_equal(sum(vos_map[, , c(3, 6, 5)]), 0)

  # for vbe, expect:
  expect_true(unique(vbe_map[, 1, 1] == c(1, 1)))
  expect_true(unique(vbe_map[, 3, 1] == c(-1, -1)))
  expect_equal(sum(vbe_map[, -c(1, 3), 1]), 0)
  expect_equal(sum(vbe_map[, , 2:6]), 0)

  # test uptake maps are generated as expected
  r1 <- c(0.25, 0.5)
  r2 <- c(0.5, 0.75)
  r2_p <- c(0.4, 0.8)
  booster_uptake <- c(0.3, 0.75)

  for (i in seq_along(r1)) {
  u <- create_uptake_map(vod_map, r1[i], r2[i], r2_p[i], booster_uptake[i])

  acc_vax <- u * vod_map

  expect_true(unique(acc_vax[, 1, 1] == c(r1[i], r1[i])))
  expect_true(unique(acc_vax[, 2, 1] == c(- (r1[i] * (1 - r2[i])),
                                          - (r1[i] * (1 - r2[i])))))
  expect_true(unique(acc_vax[, 3, 1] == c(-r1[i] * r2[i], -r1[i] * r2[i])))

  expect_true(unique(acc_vax[, 2, 2] == c(r2_p[i], r2_p[i])))
  expect_true(unique(acc_vax[, 3, 2] == -c(r2_p[i], r2_p[i])))

  expect_true(unique(acc_vax[, 4, 4] == c(booster_uptake[i],
                                          booster_uptake[i])))
  expect_true(unique(acc_vax[, 5, 4] == -c(booster_uptake[i],
                                          booster_uptake[i])))
  }


  # check VoD is working correctly --> producing the correct outputs
  # + the vaxmaps*uptakemaps generated are doing what they think they are
  # supposed to be doing

  y3e <- run_onevax_xpvwrh(tt, gp, vea = 0.5, dur_v = 1, strategy = "VoD",
                          r1 = r1, r2 = r2, r2_p = r2_p,
                          booster_uptake = booster_uptake)

  for (i in seq_along(y3e)) {

    ## all treated in X, p or W are offered vaccination
    expect_equal(y3e[[i]]$cum_offered[, , c(1, 2, 4)],
                 y3e[[i]]$cum_treated[, , c(1, 2, 4)])

    # and no-one else
    expect_equal(sum(y3e[[i]]$cum_offered[, , -c(1, 2, 4)]), 0)

    # uptake % of offered are vaccinated
    expect_equal(rowSums(y3e[[i]]$cum_offered[, , 1] * r1[i]),
                 rowSums(y3e[[i]]$cum_vaccinated[, , 1]))

    # uptake % of 1st dose offered 2nd dose are vaccinated
    expect_equal(rowSums(y3e[[i]]$cum_offered[, , 2] * r2_p[i]),
                 rowSums(y3e[[i]]$cum_vaccinated[, , 2]))

    # uptake % of offered booster are vaccinated
    expect_equal(rowSums(y3e[[i]]$cum_offered[, , 4] * booster_uptake[i]),
                 rowSums(y3e[[i]]$cum_vaccinated[, , 4]))

    # and no-one else
    expect_equal(sum(y3e[[i]]$cum_vaccinated[, , -c(1, 2, 4)]), 0)

    # no-one is lost
    expect_equal(apply(y3e[[i]]$N, 1, sum), rep(6e5, 6), tolerance = 1e-5)
  }

  # check VoA is working correctly
  y4e <- run_onevax_xpvwrh(tt, gp, vea = 0.5, dur_v = 1, strategy = "VoA",
                          r1 = r1, r2 = r2, r2_p = r2_p,
                          booster_uptake = booster_uptake)

  for (i in seq_along(y4e)) {

    ## all treated or screened in X, P, or W are offered vaccination
    expect_equal(y4e[[i]]$cum_offered[, , c(1, 2, 4)],
                 y4e[[i]]$cum_treated[, , c(1, 2, 4)] +
                   y4e[[i]]$cum_screened[, , c(1, 2, 4)])
    # and no-one else
    expect_equal(sum(y4e[[i]]$cum_offered[, , -c(1, 2, 4)]), 0)

    # uptake % of offered are vaccinated
    expect_equal(rowSums(y4e[[i]]$cum_offered[, , 1] * r1[i]),
                 rowSums(y4e[[i]]$cum_vaccinated[, , 1]))

    # uptake % of 1st dose offered 2nd dose are vaccinated
    expect_equal(rowSums(y4e[[i]]$cum_offered[, , 2] * r2_p[i]),
                 rowSums(y4e[[i]]$cum_vaccinated[, , 2]))

    # uptake % of offered booster are vaccinated
    expect_equal(rowSums(y4e[[i]]$cum_offered[, , 4] * booster_uptake[i]),
                 rowSums(y4e[[i]]$cum_vaccinated[, , 4]))
    # and no-one else
    expect_equal(sum(y4e[[i]]$cum_vaccinated[, , -c(1, 2, 4)]), 0)

    # no-one is lost
    expect_equal(apply(y4e[[i]]$N, 1, sum), rep(6e5, 6), tolerance = 1e-5)
  }

  # check vaccination targeting
  y5e <- run_onevax_xpvwrh(tt, gp, vea = 0.5, dur_v = 1,
                          strategy = "VoD(L)+VoA(H)",
                          r1 = r1, r2 = r2, r2_p = r2_p,
                          booster_uptake = booster_uptake)

  for (i in seq_along(y5e)) {
    ## L who are treated in X, P, or W are offered vaccination
    expect_equal(y5e[[i]]$cum_offered[, 1, c(1, 2, 4)],
                 y5e[[i]]$cum_treated[, 1, c(1, 2, 4)])
    ## H who are treated or screened in X, P, or W are offered vaccination
    expect_equal(y5e[[i]]$cum_offered[, 2, c(1, 2, 4)],
                 y5e[[i]]$cum_treated[, 2,  c(1, 2, 4)] +
                   y5e[[i]]$cum_screened[, 2,  c(1, 2, 4)])
    # and no-one else
    expect_equal(sum(y5e[[i]]$cum_offered[, , -c(1, 2, 4)]), 0)

    # uptake % of offered are vaccinated
    expect_equal(y5e[[i]]$cum_offered[, , 1] * r1[i],
                 y5e[[i]]$cum_vaccinated[, , 1])

    # uptake % of 1st dose offered 2nd dose are vaccinated
    expect_equal(rowSums(y5e[[i]]$cum_offered[, , 2] * r2_p[i]),
                 rowSums(y5e[[i]]$cum_vaccinated[, , 2]))

    # uptake % of offered booster are vaccinated
    expect_equal(y5e[[i]]$cum_offered[, , 4] * booster_uptake[i],
                 y5e[[i]]$cum_vaccinated[, , 4])

    # and no-one else
    expect_equal(sum(y5e[[i]]$cum_vaccinated[, , -c(1, 2, 4)]), 0)

    # no-one is lost
    expect_equal(apply(y5e[[i]]$N, 1, sum), rep(6e5, 6), tolerance = 1e-5)
  }

  # check length of uptake vector must be 1 or length(gp)
  expect_error(run_onevax_xpvwrh(tt, gp, vea = 1, dur_v = 1e3, strategy = "VbE",
                                r1 = c(0, 0.5, 1)))
  expect_error(run_onevax_xpvwrh(tt, gp, vea = 1, dur_v = 1e3, strategy = "VbE",
                                r1 = 1, r2 = c(0, 0.5, 1)))
  expect_error(run_onevax_xpvwrh(tt, gp, vea = 1, dur_v = 1e3, strategy = "VbE",
                                 r1 = 1, r2 = 1, r2_p = c(0, 0.5, 1)))
  expect_error(run_onevax_xpvwrh(tt, gp, vea = 1, dur_v = 1e3, strategy = "VbE",
                                 r1 = 1, r2 = 1, booster_uptake = c(0, 0.5, 1)))

  # check length of vea must be 1
  expect_error(run_onevax_xpvwrh(tt, gp, vea = c(0, 1, 1), dur_v = 1e3,
                                strategy = "VbE", r1 = 1))
  # check length of dur_v must be 1
  expect_error(run_onevax_xpvwrh(tt, gp, vea = 1, dur_v = c(0, 1e2, 1e3),
                                strategy = "VbE", r1 = 1))
  # check length of dur_p must be 1
  expect_error(run_onevax_xpvwrh(tt, gp, vea = 1, dur_v = 1e3,
                                 dur_p = c(0, 1e2, 1e3),
                                 strategy = "VbE", r1 = 1))
  ## test revax is working
  y6e <- run_onevax_xpvwrh(tt, gp, vea = 1, dur_v = 1e-10, dur_revax = 1e10,
                          strategy = "VoD(L)+VoA(H)",
                          r1 = r1, r2 = r2, r2_p = r2_p,
                          booster_uptake = booster_uptake)
  for (i in seq_along(y6e)) {
    ## L who are treated in X, P, or W are offered vaccination
    expect_equal(y6e[[i]]$cum_offered[, 1, c(1, 2, 4)],
                 y6e[[i]]$cum_treated[, 1, c(1, 2, 4)])
    ## H who are treated or screened in X, P, or W are offered vaccination
    expect_equal(y6e[[i]]$cum_offered[, 2, c(1, 2, 4)],
                 y6e[[i]]$cum_treated[, 2,  c(1, 2, 4)] +
                   y6e[[i]]$cum_screened[, 2,  c(1, 2, 4)])
    # and no-one else
    expect_equal(sum(y6e[[i]]$cum_offered[, , -c(1, 2, 4)]), 0)

    # uptake % of offered are vaccinated
    expect_equal(y6e[[i]]$cum_offered[, , 1] * r1[i],
                 y6e[[i]]$cum_vaccinated[, , 1])
    expect_equal(y6e[[i]]$cum_offered[, , 2] * r2_p[i],
                 y6e[[i]]$cum_vaccinated[, , 2])
    expect_equal(y6e[[i]]$cum_offered[, , 4] * booster_uptake[i],
                 y6e[[i]]$cum_vaccinated[, , 4])

    # stratum V empties immediately
    expect_equal(sum(y6e[[i]]$N[, , 3]), 0, tolerance = 1e-5)

    # efficacy is perfect
    expect_equal(sum(y6e[[i]]$cum_incid[, , c(2, 3, 5)]), 0)
    expect_true(all(y6e[[i]]$cum_incid[-1, , c(1, 4)] > 0))
  }

  y7e <- run_onevax_xpvwrh(tt, gp, vea = 0, vea_revax = 1, dur_v = 1,
                          strategy = "VoD(L)+VoA(H)",
                          r1 = r1, r2 = r2, r2_p = r2_p,
                          booster_uptake = booster_uptake,
                          hes = 0.3)

  for (i in seq_along(y7e)) {
    ## L who are treated in X, P, or W are offered vaccination
    expect_equal(y7e[[i]]$cum_offered[, 1, c(1, 2, 4)],
                 y7e[[i]]$cum_treated[, 1, c(1, 2, 4)])
    ## H who are treated or screened in X, P or W are offered vaccination
    expect_equal(y7e[[i]]$cum_offered[, 2, c(1, 2, 4)],
                 y7e[[i]]$cum_treated[, 2,  c(1, 2, 4)] +
                   y7e[[i]]$cum_screened[, 2,  c(1, 2, 4)])
    # and no-one else
    expect_equal(sum(y7e[[i]]$cum_offered[, , -c(1, 2, 4)]), 0)

    # uptake % of offered are vaccinated
    expect_equal(y7e[[i]]$cum_offered[, , 1] * r1[i],
                 y7e[[i]]$cum_vaccinated[, , 1])
    expect_equal(y7e[[i]]$cum_offered[, , 2] * r2_p[i],
                 y7e[[i]]$cum_vaccinated[, , 2])
    expect_equal(y7e[[i]]$cum_offered[, , 4] * booster_uptake[i],
                 y7e[[i]]$cum_vaccinated[, , 4])

    # efficacy is perfect in R
    expect_equal(sum(y7e[[i]]$cum_incid[, , 5]), 0)
    # but not in XPWVH
    expect_true(all(y7e[[i]]$cum_incid[-1, , -c(2, 5)] > 0))
  }

  ## test restart with hesitancy is working

  y8 <- run_onevax_xpvwrh(tt, gp, vea = 0, dur_v = 1e3)

  i_p <- lapply(y8, restart_hes, n_vax = 6, hes = 0.5, branching = TRUE)
  y_hesres <- run_onevax_xpvwrh(tt, gp, init_params = i_p, vea = 0, dur_v = 1e3,
                               hes = 0.5)

  # final timepoint y8 run = 2 * first timepoint of y_hesres run (as hes = 0.5)
  # (for XVWR)
  for (i in seq_along(y8)) {
    expect_equal(y8[[i]]$U[length(tt), , 1:5], y_hesres[[i]]$U[1, , 1:5] * 2)
    expect_equal(y8[[i]]$I[length(tt), , 1:5], y_hesres[[i]]$I[1, , 1:5] * 2)
    expect_equal(y8[[i]]$A[length(tt), , 1:5], y_hesres[[i]]$A[1, , 1:5] * 2)
    expect_equal(y8[[i]]$S[length(tt), , 1:5], y_hesres[[i]]$S[1, , 1:5] * 2)
    expect_equal(y8[[i]]$T[length(tt), , 1:5], y_hesres[[i]]$T[1, , 1:5] * 2)
  }

  # restart_hes moves X -> H only, and correctly
  for (i in seq_along(y_hesres)) {
    expect_equal(y_hesres[[i]]$U[, , 1], y_hesres[[i]]$U[, , 6])
    expect_equal(y_hesres[[i]]$I[, , 1], y_hesres[[i]]$I[, , 6])
    expect_equal(y_hesres[[i]]$A[, , 1], y_hesres[[i]]$A[, , 6])
    expect_equal(y_hesres[[i]]$S[, , 1], y_hesres[[i]]$S[, , 6])
    expect_equal(y_hesres[[i]]$T[, , 1], y_hesres[[i]]$T[, , 6])

    expect_equal(rowSums(y8[[i]]$U[, , 2:5]), rep(0, length(tt)))
    expect_equal(rowSums(y8[[i]]$I[, , 2:5]), rep(0, length(tt)))
    expect_equal(rowSums(y8[[i]]$A[, , 2:5]), rep(0, length(tt)))
    expect_equal(rowSums(y8[[i]]$S[, , 2:5]), rep(0, length(tt)))
    expect_equal(rowSums(y8[[i]]$T[, , 2:5]), rep(0, length(tt)))

  }

  # restart_hes gives error if baseline model run provided already has hesitancy

  y9 <- run_onevax_xpvwrh(tt, gp, vea = 0, dur_v = 1e3, hes = 0.5)

  expect_error(lapply(y9, restart_hes, n_vax = 6, hes = 0.5, branching = TRUE))

  # restart_hes gives error if baseline model run provided contains vaccinated

  y10 <- run_onevax_xpvwrh(tt, gp, vea = 0, dur_v = 1e3, vbe = 1)

  expect_error(lapply(y10, restart_hes, n_vax = 6, hes = 0.5, branching = TRUE))

  # check booster_uptake defaults to r1r2 primary uptake
  y_r1r2_only <- run_onevax_xpvwrh(tt, gp, vea = 0, dur_v = 1,
                                  r1 = 0.8, r2 = 0.5,
                                  strategy = "VoD(L)+VoA(H)")

  y_r1r2_boost <- run_onevax_xpvwrh(tt, gp, vea = 0, dur_v = 1,
                                  r1 = 0.8, r2 = 0.5,
                                  booster_uptake = 0.4,
                                  strategy = "VoD(L)+VoA(H)")

  expect_equal(y_r1r2_only, y_r1r2_boost)

  # When everyone receives one dose (r2 = 0), V, W R, H all empty
  y11 <- run_onevax_xpvwrh(tt, gp, r1 = 1, r2 = 0, strategy = "VoD")

  for (i in seq_along(y11)) {

  expect_equal(sum(y11[[i]]$N[, , -c(1, 2)]), 0)
  expect_true(all(rowSums(y11[[i]]$N[, , 1]) > 0))
  expect_true(all(rowSums(y11[[i]]$N[-1, , 2]) > 0))

  }

  # When everyone receives two doses straight away (r1*r2 = 1),
  # P and H all empty

  y12 <- run_onevax_xpvwrh(tt, gp, r1 = 1, r2 = 1, strategy = "VoD")

  for (i in seq_along(y12)) {

    expect_equal(sum(y12[[i]]$N[, , c(2, 6)]), 0)
    expect_true(all(rowSums(y12[[i]]$N[, , 1]) > 0))
    expect_true(all(rowSums(y12[[i]]$N[-1, , 3:5]) > 0))

  }

  # if half the population gets 1 dose and half the population gets 2 doses
  # and durations are the same
  # then N should be the same for P and V

  y13 <- run_onevax_xpvwrh(tt, gp, r1 = 1, r2 = 0.5, strategy = "VoD")

  for (i in seq_along(y13)) {
  expect_equal(y13[[i]]$N[, , 2], y13[[i]]$N[, , 3])
  }
    # incidence should also be the same as all vex_p default to vex

  for (i in seq_along(y13)) {
  expect_equal(y13[[i]]$cum_incid[, , 2], y13[[i]]$cum_incid[, , 3])
  }

  # test vex_p work independently of vex
  y14 <- run_onevax_xpvwrh(tt, gp, r1 = 1, r2 = 0.5, vea_p = 1,
                           strategy = "VoD")

    # incidence in V but no incidence in P as vea_p = 1
    # individuals present in both

  for (i in seq_along(y14)) {
    expect_true(all(rowSums(y14[[i]]$N[-1, , 2]) > 0))
    expect_true(all(rowSums(y14[[i]]$N[-1, , 3]) > 0))
    expect_true(sum(y14[[i]]$cum_incid[-1, , 2]) == 0)
    expect_true(all(rowSums(y14[[i]]$cum_incid[-1, , 3]) > 0))
  }

  # test vea_p defaults to vea

  y_vea_only <- run_onevax_xpvwrh(tt, gp, r1 = 1, r2 = 0.5, vea = 0.5,
                           strategy = "VoD")
  y_vea_veap <- run_onevax_xpvwrh(tt, gp, r1 = 1, r2 = 0.5, vea = 0.5,
                                  vea_p = 0.5, strategy = "VoD")

  expect_equal(y_vea_only, y_vea_veap)

  # when vea, vea_p and vea_revax = 1 there is no incidence

  y_no_incid <- run_onevax_xpvwrh(tt, gp, r1 = 1, r2 = 0.5, vea = 1)

  for (i in seq_along(y_no_incid)) {
  expect_equal(sum(y_no_incid[[i]]$cum_incid[, , -1]), 0)
  }

  # check initial params generated correctly when coverage is specified
  pars <- lapply(gp[1], model_params)
  init_1 <- initial_params_xpvwrh(pars[[1]],
                                  coverage_p = 0.25, coverage_v = 0.25)

  # expect same number starting in P as in V
  expect_equal(init_1$U0[, 2], init_1$U0[, 3])

  init_2 <- initial_params_xpvwrh(pars[[1]],
                                 coverage_p = 0.5, coverage_v = 0.25)

  # expect double starting in P as in V
  expect_equal((init_2$U0[, 2]), init_2$U0[, 3] * 2)

  # expect error if sum of coverages and hes is greater than 1

  expect_error(initial_params_xpvwrh(pars[[1]], coverage_p = 0.5,
                                     coverage_v = 0.5, hes = 0.5))

  # expect error if sum of coverages is 100%
  # (if total coverage is 100%, no asymptomatic infection is seeded)

  expect_error(initial_params_xpvwrh(pars[[1]], coverage_p = 0.5,
                                     coverage_v = 0.5))

  # tests correct number of individuals are moving to V and to P
  # duration is large for v and p so assuming no waning

  y15 <- run_onevax_xpvwrh(tt, gp, r1 = 0.5, r2 = 0.5, dur_v = 1e90,
                           strategy = "VoD")

    # number getting vaccinated from X (timepoint 2 as all is 0 in timepoint 1)
    cum_vac <- sum(y15[[1]]$cum_vaccinated[2, , 1])
    cum_off <- sum(y15[[1]]$cum_offered[2, , 1])

    # 50% of individuals offered will accept and be vaccinated (0.25 + 0.25)
    prop_accept <- cum_vac / cum_off
    expect_equal(prop_accept, 0.5)

    # 50% of those treated, are vaccinated
    sum(y15[[1]]$cum_vaccinated[2, , 1])  / sum(y15[[1]]$cum_treated[2, , 1])

    # number total who received vaccination
    # i.e number in P and V at timepoint 2 given no waning
    p <- sum(y15[[1]]$N[2, , 2])
    v <- sum(y15[[1]]$N[2, , 3])

    # should be same number in P as in V
    expect_equal(p, v)

  # test individuals in P are actually waning back into X
    # move all individuals into P for starting conditions
    # no vaccination, duration of vaccine short so waning is fast
    # should see individuals accumulating in X, decreasing in P

  pars <- lapply(gp[1], model_params)
  ip <- lapply(pars, initial_params_xpvwrh, coverage_p = 0.99999999999999,
               t = 5)

  y16 <- run_onevax_xpvwrh(tt, gp, init_params = ip, dur_p = 1e-90)

    # entire population starts in P then wanes to X and stays there
  for (i in seq_along(y16)) {
      expect_equal(sum(y16[[i]]$N[1, , 2]), 6e+05)
      expect_equal(rowSums(y16[[i]]$N[2:(max(tt) + 1), , 1]),
                   rep(6e+05, max(tt)))

    # V W R are empty for all time
      expect_equal(rowSums(y16[[i]]$N[, , 3:6]), rep(0, 6))
}

  # test individuals still wane to W
    ip <- lapply(pars, initial_params_xpvwrh, coverage_v = 0.99999999999999
                 , t = 5)
    y17 <- run_onevax_xpvwrh(tt, gp, init_params = ip, dur_v = 1e-90)

    # entire population starts in V then wanes to W and stays there
    # note, W won't be equal to 6e+05 as 1. individuals in W die but 2. entrants
    # enter into X 3. no vaccination so replenishment of V and subsequently W

    for (i in seq_along(y17)) {
      expect_equal(sum(y17[[i]]$N[1, , 3]), 6e+05)
      expect_equal(rowSums(y17[[i]]$N[2:(max(tt) + 1), , 3]), rep(0, max(tt)))

      for (j in 2:max(tt) + 1) {
      expect_true(sum(y17[[i]]$N[j, , 4]) > 0)
      }

                                                 ######## difference in W
                                              ######### isnt quite 1200 though?
    # R, W are empty for all time
      expect_equal(rowSums(y17[[i]]$N[, , 5:6]), rep(0, 6))
}

})
