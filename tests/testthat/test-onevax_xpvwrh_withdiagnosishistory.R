#testing xpvwrh model

test_that("run_onevax_xpvwrh works correctly", {

  for (j in 1:5) {

    n_diag_rec <- j

    tt <- seq(0, 5)
    gp <- gono_params(1:2)
    y1 <- run_onevax_xpvwrh(tt, gp, vea = 0, dur_v = 1e3, n_erlang = 1,
                            n_diag_rec = n_diag_rec)[[1]]

    # check no-one is vaccinated with v switched off
    expect_true(all(y1$cum_vaccinated == 0))

    # check 100% vbe vaccinates all new entrants
    y2 <- run_onevax_xpvwrh(tt, gp, vea = 0, dur_v = 1e10, vbe = 1,
                            n_diag_rec = n_diag_rec)

    # cum_vaccinated = 12000 each year = number of entrants
    expect_equal(diff(rowSums(y2[[1]]$cum_vaccinated[, , 1])), rep(12e3,
                                                                   max(tt)))

    # this also means there should be people moving into 'V'
    expect_true(sum(y2[[1]]$N[, , 2 * n_diag_rec + 1]) > 0)

    # but not P
    expect_equal(sum(y2[[1]]$N[, , n_diag_rec + 1]), 0)

    # other compartments empty
    expect_equal(sum(y2[[1]]$N[, , -c(1:n_diag_rec,
                                      (2 * n_diag_rec + 1):(3 * n_diag_rec),
                                      (3 * n_diag_rec + 1):(4 * n_diag_rec))]),
                 0)

    # and no-one else
    expect_equal(sum(y2[[1]]$cum_vaccinated[, , 2:(6 * n_diag_rec)]), 0)
    expect_equal(sum(y2[[1]]$N[, , -c(1:n_diag_rec,
                                      (2 * n_diag_rec + 1):(3 * n_diag_rec),
                                      (3 * n_diag_rec + 1):(4 * n_diag_rec))]),
                 0)


    # check population is constant in size over time

    for (i in seq_along(y2))  {
      expect_equal(rowSums(y2[[i]]$N[, , ]), rep(6e+05, max(tt) + 1))
    }

    # check this is still the case when hesitancy > 0
    y2.1 <- run_onevax_xpvwrh(tt, gp, vea = 0, dur_v = 1e3, vbe = 1, hes = 0.5,
                              n_diag_rec = n_diag_rec)

    for (i in seq_along(y2.1)) {
      expect_equal(diff(rowSums(y2.1[[i]]$cum_vaccinated[, , 1])),
                   rep((12e3) / 2, max(tt)))

      expect_true(sum(y2.1[[i]]$N[, , (2 * n_diag_rec + 1)]) > 0)
      # check no vaccination in hesitant entrants for vbe = 100%

      expect_equal((rowSums(y2.1[[i]]$cum_vaccinated[, , 5 * n_diag_rec + 1])),
                   rep(0, max(tt) + 1))
    }
    # expect error if inputs are not of length 1 or equal to length of params

    uptake <- c(0, 2.5, 0.5, 0.75, 1)
    expect_error(y <- run_onevax_xpvwrh(tt, gp, vea = 0, dur_v = 1e3, vbe = 1,
                                        primary_uptake = uptake,
                                        n_diag_rec = n_diag_rec))

    # check can restart
    init_params <- lapply(y2, restart_params)
    y3 <- run_onevax_xpvwrh(seq(max(tt), length.out = 2, by = 1),
                            gp, init_params, vea = 0, dur_v = 1e3, vbe = 1,
                            n_diag_rec = n_diag_rec)

    for (i in seq_along(y2)) {
      expect_equal(y2[[i]]$U[length(tt), , ], y3[[i]]$U[1, , ])
      expect_equal(y2[[i]]$I[length(tt), , ], y3[[i]]$I[1, , ])
      expect_equal(y2[[i]]$A[length(tt), , ], y3[[i]]$A[1, , ])
      expect_equal(y2[[i]]$S[length(tt), , ], y3[[i]]$S[1, , ])
      expect_equal(y2[[i]]$T[length(tt), , ], y3[[i]]$T[1, , ])
    }

    y_h <- run_onevax_xpvwrh(tt, gp, vea = 0, dur_v = 1e3, vbe = 1, hes = 0.3,
                             n_diag_rec = n_diag_rec)

    # initial population split between non-vaccinated and hesitant only
    # other stratum empty
    init_pop <- 600000 * c(0.85, 0.15)
    for (i in seq_along(y_h)) {
      expect_equivalent(y_h[[i]]$N[1, , 1], init_pop * 0.7)
      expect_equivalent(y_h[[i]]$N[1, , 5 * n_diag_rec + 1], init_pop * 0.3)
      expect_true(all(y_h[[i]]$N[1, , (2:(5 * n_diag_rec))] == 0))
      if (n_diag_rec > 1) {
        expect_true(all(y_h[[i]]$N[1, , ((5 * n_diag_rec + 2):(6 * n_diag_rec))]
                        == 0))
      }
    }

    # incidence in hesitant population increasing
    for (i in seq_along(y_h)) {
      expect_true(all(y_h[[i]]$cum_incid[-1, , (5 * n_diag_rec + 1):
                                           (6 * n_diag_rec)] > 0))
    }

    y_h2 <- run_onevax_xpvwrh(tt, gp, vea = 0, dur_v = 1e3, vbe = 0, hes = 0.3,
                              n_diag_rec = n_diag_rec)

    for (i in seq_along(y_h2)) {
      # yearly population entrants enter X and H strata
      # in accordance with assigned proportion of hesitancy 'hes'
      # H and X stratum sizes remain constant in time

      if (n_diag_rec == 1) {
        expect_equal(y_h2[[i]]$N[1, , ], y_h2[[i]]$N[6, , ])
      } else if (n_diag_rec >= 2) {

        temp1 <- c(rowSums(y_h2[[i]]$N[1, , 1:n_diag_rec]),
                   rowSums(y_h2[[i]]$N[1, , 5 * n_diag_rec + (1:n_diag_rec)]))
        temp2 <- c(rowSums(y_h2[[i]]$N[6, , 1:n_diag_rec]),
                   rowSums(y_h2[[i]]$N[6, , 5 * n_diag_rec + (1:n_diag_rec)]))

        expect_equal(temp1, temp2)
      }
    }

    # H and X stratum equal when no vaccination and hes = 0.5
    y_h3 <- run_onevax_xpvwrh(tt, gp, vea = 0, dur_v = 1e3, hes = 0.5,
                              n_diag_rec = n_diag_rec)

    for (i in seq_along(y_h3)) {

      # Number of infections in X and H equal for no vaccination and hes = 0.5
      for (jj in 1:n_diag_rec) {
        expect_equal(y_h3[[i]]$N[, , jj], y_h3[[i]]$N[, , 5 * n_diag_rec + jj])
        expect_equal(y_h3[[i]]$cum_incid[, , jj],
                     y_h3[[i]]$cum_incid[, , 5 * n_diag_rec + jj])
      }

      expect_equal(y_h3[[i]]$cum_incid[, , 1],
                   y_h3[[i]]$cum_incid[, , 5 * n_diag_rec + 1])

      if (n_diag_rec == 2) {
        expect_equal(y_h3[[i]]$cum_incid[, , 2],
                     y_h3[[i]]$cum_incid[, , 5 * n_diag_rec + 2])
      }
    }
  }
})


test_that("The vaccination maps are being generated as expected", {

  gp <- gono_params(1:2)
  pars <- lapply(gp[1], model_params)
  vbe <- 1
  p <- set_strategy(strategy = "VoD(L)+VoA(H)", vbe > 0)
  n_erlang <- 1

  for (i in 1:5) {

    n_diag_rec <- i

    idx <- stratum_index_xpvwrh(n_erlang = n_erlang, n_diag_rec)

    if (n_diag_rec == 1) {
      i_eligible <- c(1, 1, 2, 4)
    } else if (n_diag_rec == 2) {
      i_eligible <- c(1, 2, 1, 2, 3, 4, 7, 8)
    } else {
      i_eligible <- c(1:n_diag_rec, 1:n_diag_rec, (n_diag_rec + (1:n_diag_rec)),
                      3 * n_diag_rec + (1:n_diag_rec))
    }

    expect_equal(i_eligible,  c(idx$X, idx$X, idx$P, idx$W))

    if (n_diag_rec == 1) {
      i_p <- c(2, 3, 3, 5)
    }else if (n_diag_rec == 2) {
      i_p <- c(3, 4, 5, 6, 5, 6, 9, 10)
    } else {
      i_p <- c((n_diag_rec + (1:n_diag_rec)), 2 * n_diag_rec + (1:n_diag_rec),
               2 * n_diag_rec + (1:n_diag_rec), 4 * n_diag_rec + (1:n_diag_rec))
    }

    i_p_temp_2 <- i_p
    i_p_temp_2[i_p %% n_diag_rec != 0] <- i_p[i_p %% n_diag_rec != 0] + 1

    expect_equal(i_p, c(idx$P[(1:n_diag_rec)], rep(idx$V[(1:n_diag_rec)],
                                                   n_erlang + 1),
                        idx$R[1:n_diag_rec]))

    vod_map <- create_vax_map_branching(n_vax = idx$n_vax, p$vod, i_eligible,
                                        i_p_temp_2, idx = idx)
    vos_map <- create_vax_map_branching(n_vax = idx$n_vax, p$vos, i_eligible,
                                        i_p, idx = idx)
    vbe_map <- create_vax_map_branching(n_vax = idx$n_vax, p$vbe, i_eligible,
                                        i_p, set_vbe = TRUE, idx = idx)

    for (d in (1:n_diag_rec)) {

      if (d < n_diag_rec) {

        expect_true(unique(vod_map[, d, d] == c(1, 1)))
        expect_true(unique(vod_map[, n_diag_rec + d + 1, d] == c(-1, -1)))
        expect_true(unique(vod_map[, 2 * n_diag_rec + d + 1, d] == c(-1, -1)))
        expect_true(unique(vod_map[, n_diag_rec + d, n_diag_rec + d] ==
                             c(1, 1)))
        expect_true(unique(vod_map[, 2 * n_diag_rec + d + 1, n_diag_rec + d] ==
                             c(-1, -1)))
        expect_true(unique(vod_map[, 3 * n_diag_rec + d, 3 * n_diag_rec + d] ==
                             c(1, 1)))
        expect_true(unique(vod_map[, 4 * n_diag_rec + d + 1, 3 * n_diag_rec + d]
                           == c(-1, -1)))
      } else {

        expect_true(unique(vod_map[, d, d] == c(1, 1)))
        expect_true(unique(vod_map[, n_diag_rec + d, d] == c(-1, -1)))
        expect_true(unique(vod_map[, 2 * n_diag_rec + d, d] == c(-1, -1)))
        expect_true(unique(vod_map[, n_diag_rec + d, n_diag_rec + d] ==
                             c(1, 1)))
        expect_true(unique(vod_map[, 2 * n_diag_rec + d, n_diag_rec + d] ==
                             c(-1, -1)))
        expect_true(unique(vod_map[, 3 * n_diag_rec + d, 3 * n_diag_rec + d] ==
                             c(1, 1)))
        expect_true(unique(vod_map[, 4 * n_diag_rec + d, 3 * n_diag_rec + d] ==
                             c(-1, -1)))
      }
    }

    for (d in (1:n_diag_rec)) {
      if (d < n_diag_rec) {
        expect_equal(sum(vod_map[, -c(d, n_diag_rec + d + 1,
                                      2 * n_diag_rec + d + 1), d]), 0)
        expect_equal(sum(vod_map[, -c(n_diag_rec + d, 2 * n_diag_rec + d + 1),
                                 n_diag_rec + d]), 0)
        expect_equal(sum(vod_map[, -c(3 * n_diag_rec + d,
                                      4 * n_diag_rec + d + 1),
                                 3 * n_diag_rec + d]), 0)
        expect_equal(sum(vod_map[, , c(2 * n_diag_rec + d + 1,
                                       5 * n_diag_rec + d + 1,
                                       4 * n_diag_rec + d)]), 0)
      } else {
        expect_equal(sum(vod_map[, -c(d, n_diag_rec + d, 2 * n_diag_rec + d),
                                 d]), 0)
        expect_equal(sum(vod_map[, -c(n_diag_rec + d, 2 * n_diag_rec + d),
                                 n_diag_rec + d]), 0)
        expect_equal(sum(vod_map[, -c(3 * n_diag_rec + d, 4 * n_diag_rec + d),
                                 3 * n_diag_rec + d]), 0)
        expect_equal(sum(vod_map[, , c(2 * n_diag_rec + d, 5 * n_diag_rec + d,
                                       4 * n_diag_rec + d)]), 0)
      }
    }

    # for vos, expect:
    for (d in (1:n_diag_rec)) {
      expect_true(unique(vos_map[, d, d] == c(0, 1)))
      expect_true(unique(vos_map[, n_diag_rec + d, d] == c(0, -1)))
      expect_true(unique(vos_map[, 2 * n_diag_rec + d, d] == c(0, -1)))
      expect_true(unique(vos_map[, n_diag_rec + d, n_diag_rec + d] == c(0, 1)))
      expect_true(unique(vos_map[, 2 * n_diag_rec + d, n_diag_rec + d] ==
                           c(0, -1)))
      expect_true(unique(vos_map[, 3 * n_diag_rec + d, 3 * n_diag_rec + d] ==
                           c(0, 1)))
      expect_true(unique(vos_map[, 4 * n_diag_rec + d, 3 * n_diag_rec + d] ==
                           c(0, -1)))
    }

    for (d in (1:n_diag_rec)) {
      expect_equal(sum(vos_map[, -c(d, n_diag_rec + d, 2 * n_diag_rec + d), d]),
                   0)
      expect_equal(sum(vos_map[, -c(n_diag_rec + d, 2 * n_diag_rec + d),
                               n_diag_rec + d]), 0)
      expect_equal(sum(vos_map[, -c(3 * n_diag_rec + d, 4 * n_diag_rec + d),
                               3 * n_diag_rec + d]), 0)
      expect_equal(sum(vos_map[, , c(2 * n_diag_rec + d, 5 * n_diag_rec + d,
                                     4 * n_diag_rec + d)]), 0)
    }

    # for vbe, expect:
    expect_true(unique(vbe_map[, 1, 1] == c(1, 1)))
    expect_true(unique(vbe_map[, 2 * n_diag_rec + 1, 1] == c(-1, -1)))
    expect_equal(sum(vbe_map[, -c(1, 2 * n_diag_rec + 1), 1]), 0)
    expect_equal(sum(vbe_map[, , 2:(6 * n_diag_rec)]), 0)

    # test uptake maps are generated as expected
    r1 <- c(0.25, 0.5)
    r2 <- c(0.5, 0.75)
    r2_p <- c(0.4, 0.8)
    booster_uptake <- c(0.3, 0.75)

    for (i in seq_along(r1)) {
      u <- create_uptake_map_xpvwrh(vod_map, r1[i], r2[i], r2_p[i],
                                    booster_uptake[i], idx,
                                    n_diag_rec = n_diag_rec, n_erlang = 1,
                                    screening_or_diagnosis = "diagnosis")

      acc_vax <- u * vod_map

      for (d in (1:n_diag_rec)) {

        if (d < n_diag_rec) {

          expect_true(unique(acc_vax[, d, d] == c(r1[i], r1[i])))
          expect_true(unique(acc_vax[, n_diag_rec + d + 1, d] ==
                               c(- (r1[i] * (1 - r2[i])),
                                 - (r1[i] * (1 - r2[i])))))

          expect_true(unique(acc_vax[, 2 * n_diag_rec + d + 1, d] ==
                               c(-r1[i] * r2[i], -r1[i] * r2[i])))

          expect_true(unique(acc_vax[, n_diag_rec + d, n_diag_rec + d] ==
                               c(r2_p[i], r2_p[i])))

          expect_true(unique(acc_vax[, 2 * n_diag_rec + d + 1, n_diag_rec + d]
                             == -c(r2_p[i], r2_p[i])))

          expect_true(unique(acc_vax[, 3 * n_diag_rec + d, 3 * n_diag_rec + d]
                             == c(booster_uptake[i], booster_uptake[i])))

          expect_true(unique(acc_vax[, 4 * n_diag_rec + d + 1,
                                     3 * n_diag_rec + d] ==
                               -c(booster_uptake[i], booster_uptake[i])))

        } else {
          expect_true(unique(acc_vax[, d, d] == c(r1[i], r1[i])))
          expect_true(unique(acc_vax[, n_diag_rec + d, d] ==
                               c(- (r1[i] * (1 - r2[i])),
                                 - (r1[i] * (1 - r2[i])))))
          expect_true(unique(acc_vax[, 2 * n_diag_rec + d, d] ==
                               c(-r1[i] * r2[i], -r1[i] * r2[i])))

          expect_true(unique(acc_vax[, n_diag_rec + d, n_diag_rec + d] ==
                               c(r2_p[i], r2_p[i])))
          expect_true(unique(acc_vax[, 2 * n_diag_rec + d, n_diag_rec + d] ==
                               -c(r2_p[i], r2_p[i])))

          expect_true(unique(acc_vax[, 3 * n_diag_rec + d, 3 * n_diag_rec + d]
                             == c(booster_uptake[i], booster_uptake[i])))

          expect_true(unique(acc_vax[, 4 * n_diag_rec + d, 3 * n_diag_rec + d]
                             == -c(booster_uptake[i], booster_uptake[i])))
        }

      }
    }
  }
})

test_that("Test VoD is working correctly", {
  # check VoD is working correctly --> producing the correct outputs
  # + the vaxmaps*uptakemaps generated are doing what they think they are
  # supposed to be doing

  tt <- seq(0, 5)
  gp <- gono_params(1:2)
  r1 <- c(0.25, 0.5)
  r2 <- c(0.5, 0.75)
  r2_p <- c(0.4, 0.8)
  booster_uptake <- c(0.3, 0.75)

  for (j in 1:5) {
    n_diag_rec <- j

    y3e <- run_onevax_xpvwrh(tt, gp, vea = 0.5, dur_v = 1, strategy = "VoD",
                             r1 = r1, r2 = r2, r2_p = r2_p,
                             booster_uptake = booster_uptake,
                             n_diag_rec = n_diag_rec)

    for (i in seq_along(y3e)) {

      ## all treated in X, p or W are offered vaccination
      expect_equal(y3e[[i]]$cum_offered[, ,
                                        c(1:n_diag_rec,
                                          n_diag_rec + 1:n_diag_rec,
                                          3 * n_diag_rec + 1:n_diag_rec)],
                   y3e[[i]]$cum_treated[, , c(1:n_diag_rec,
                                              n_diag_rec + 1:n_diag_rec,
                                              3 * n_diag_rec + 1:n_diag_rec)])

      # uptake % of offered are vaccinated
      expect_equal(rowSums(y3e[[i]]$cum_offered[, , 1:n_diag_rec] * r1[i]),
                   rowSums(y3e[[i]]$cum_vaccinated[, , 1:n_diag_rec]))

      # uptake % of 1st dose offered 2nd dose are vaccinated
      expect_equal(rowSums(y3e[[i]]$cum_offered[, , n_diag_rec + 1:n_diag_rec]
                           * r2_p[i]),
                   rowSums(y3e[[i]]$cum_vaccinated[, ,
                                                   n_diag_rec + 1:n_diag_rec]))

      # uptake % of offered booster are vaccinated
      expect_equal(rowSums(y3e[[i]]$cum_offered[, ,
                                                3 * n_diag_rec + 1:n_diag_rec]
                           * booster_uptake[i]),
                   rowSums(y3e[[i]]$cum_vaccinated[, ,
                                                   3 * n_diag_rec +
                                                     1:n_diag_rec]))

      # and no-one else
      expect_equal(sum(y3e[[i]]$cum_vaccinated[, , -c(1:n_diag_rec,
                                                      n_diag_rec + 1:n_diag_rec,
                                                      3 * n_diag_rec +
                                                        1:n_diag_rec)]), 0)

      ## all treated in X, p or W are offered vaccination
      # and no-one else
      expect_equal(sum(y3e[[i]]$cum_offered[, , -c(1:n_diag_rec,
                                                   n_diag_rec + 1:n_diag_rec,
                                                   3 * n_diag_rec +
                                                     1:n_diag_rec)]), 0)

      # uptake % of offered booster are vaccinated
      # and no-one else
      expect_equal(sum(y3e[[i]]$cum_vaccinated[, , -c(1:n_diag_rec,
                                                      n_diag_rec + 1:n_diag_rec,
                                                      3 * n_diag_rec +
                                                        1:n_diag_rec)]), 0)

      # no-one is lost
      expect_equal(apply(y3e[[i]]$N, 1, sum), rep(6e5, 6), tolerance = 1e-5)
    }
  }
})


test_that("Test VoA is working correctly", {
  # check VoA is working correctly

  tt <- seq(0, 5)
  gp <- gono_params(1:2)
  r1 <- c(0.25, 0.5)
  r2 <- c(0.5, 0.75)
  r2_p <- c(0.4, 0.8)
  booster_uptake <- c(0.3, 0.75)

  for (j in 1:5) {

    n_diag_rec <- j

    y4e <- run_onevax_xpvwrh(tt, gp, vea = 0.5, dur_v = 1, strategy = "VoA",
                             r1 = r1, r2 = r2, r2_p = r2_p,
                             booster_uptake = booster_uptake,
                             n_diag_rec = n_diag_rec)

    for (i in seq_along(y4e)) {

      ## all treated or screened in X, P, or W are offered vaccination
      expect_equal(y4e[[i]]$cum_offered[, , c(1:n_diag_rec,
                                              n_diag_rec + 1:n_diag_rec,
                                              3 * n_diag_rec +
                                                (1:n_diag_rec))],
                   y4e[[i]]$cum_treated[, , c(1:n_diag_rec,
                                              n_diag_rec + 1:n_diag_rec,
                                              3 * n_diag_rec +
                                                (1:n_diag_rec))] +
                     y4e[[i]]$cum_screened[, , c(1:n_diag_rec,
                                                 n_diag_rec + 1:n_diag_rec,
                                                 3 * n_diag_rec +
                                                   (1:n_diag_rec))])

      # and no-one else
      expect_equal(sum(y4e[[i]]$cum_offered[, , -c(1:n_diag_rec,
                                                   n_diag_rec +
                                                     (1:n_diag_rec),
                                                   3 * n_diag_rec +
                                                     (1:n_diag_rec))]), 0)

      # uptake % of offered are vaccinated
      expect_equal(rowSums(y4e[[i]]$cum_offered[, , 1:n_diag_rec] * r1[i]),
                   rowSums(y4e[[i]]$cum_vaccinated[, , 1:n_diag_rec]))

      # uptake % of 1st dose offered 2nd dose are vaccinated
      expect_equal(rowSums(y4e[[i]]$cum_offered[, , n_diag_rec + 1:n_diag_rec]
                           * r2_p[i]),
                   rowSums(y4e[[i]]$cum_vaccinated[, , n_diag_rec
                                                   + 1:n_diag_rec]))

      # uptake % of offered booster are vaccinated
      expect_equal(rowSums(y4e[[i]]$cum_offered[, , 3 * n_diag_rec
                                                + 1:n_diag_rec]
                           * booster_uptake[i]),
                   rowSums(y4e[[i]]$cum_vaccinated[, , 3 * n_diag_rec
                                                   + 1:n_diag_rec]))
      # and no-one else
      expect_equal(sum(y4e[[i]]$cum_vaccinated[, , -c(1:n_diag_rec,
                                                      n_diag_rec +
                                                        1:n_diag_rec,
                                                      3 * n_diag_rec +
                                                        1:n_diag_rec)]), 0)
      # no-one is lost
      expect_equal(apply(y4e[[i]]$N, 1, sum), rep(6e5, 6), tolerance = 1e-5)
    }
  }
})


test_that("Test vaccination targeting", {

  tt <- seq(0, 5)
  gp <- gono_params(1:2)
  r1 <- c(0.25, 0.5)
  r2 <- c(0.5, 0.75)
  r2_p <- c(0.4, 0.8)
  booster_uptake <- c(0.3, 0.75)

  for (j in 1:5) {

    n_diag_rec <- j

    # check vaccination targeting
    y5e <- run_onevax_xpvwrh(tt, gp, vea = 0.5, dur_v = 1,
                             strategy = "VoD(L)+VoA(H)",
                             r1 = r1, r2 = r2, r2_p = r2_p,
                             booster_uptake = booster_uptake,
                             n_diag_rec = n_diag_rec)

    for (i in seq_along(y5e)) {
      ## L who are treated in X, P, or W are offered vaccination
      expect_equal(y5e[[i]]$cum_offered[, 1,
                                        c(1:n_diag_rec, n_diag_rec +
                                            1:n_diag_rec,
                                          3 * n_diag_rec + (1:n_diag_rec))],
                   y5e[[i]]$cum_treated[, 1,
                                        c(1:n_diag_rec, n_diag_rec +
                                            1:n_diag_rec, 3 * n_diag_rec +
                                            (1:n_diag_rec))])
      ## H who are treated or screened in X, P, or W are offered vaccination
      expect_equal(y5e[[i]]$cum_offered[, 2, c(1:n_diag_rec, n_diag_rec +
                                                 1:n_diag_rec, 3 * n_diag_rec +
                                                 (1:n_diag_rec))],
                   y5e[[i]]$cum_treated[, 2,  c(1:n_diag_rec, n_diag_rec +
                                                  1:n_diag_rec, 3 * n_diag_rec +
                                                  (1:n_diag_rec))] +
                     y5e[[i]]$cum_screened[, 2, c(1:n_diag_rec,
                                                  n_diag_rec + 1:n_diag_rec,
                                                  3 * n_diag_rec +
                                                    (1:n_diag_rec))])
      # and no-one else
      expect_equal(sum(y5e[[i]]$cum_offered[, , -c(1:n_diag_rec,
                                                   n_diag_rec + 1:n_diag_rec,
                                                   3 * n_diag_rec +
                                                     (1:n_diag_rec))]), 0)

      # uptake % of offered are vaccinated
      expect_equal(y5e[[i]]$cum_offered[, , 1:n_diag_rec] * r1[i],
                   y5e[[i]]$cum_vaccinated[, , 1:n_diag_rec])

      # uptake % of 1st dose offered 2nd dose are vaccinated
      expect_equal(rowSums(y5e[[i]]$cum_offered[, , n_diag_rec + 1:n_diag_rec]
                           * r2_p[i]),
                   rowSums(y5e[[i]]$cum_vaccinated[, , n_diag_rec +
                                                     1:n_diag_rec]))

      # uptake % of offered booster are vaccinated
      expect_equal(y5e[[i]]$cum_offered[, , 3 * n_diag_rec + 1:n_diag_rec]
                   * booster_uptake[i],
                   y5e[[i]]$cum_vaccinated[, , 3 * n_diag_rec + 1:n_diag_rec])

      # and no-one else
      expect_equal(sum(y5e[[i]]$cum_vaccinated[, , -c(1:n_diag_rec, n_diag_rec +
                                                        1:n_diag_rec,
                                                      3 * n_diag_rec +
                                                        (1:n_diag_rec))]), 0)

      # no-one is lost
      expect_equal(apply(y5e[[i]]$N, 1, sum), rep(6e5, 6), tolerance = 1e-5)
    }
  }
})

test_that("Test vaccination according to history is working correctly", {

  tt <- seq(0, 5)
  gp <- gono_params(1:2)
  r1 <- c(0.25, 0.5)
  r2 <- c(0.5, 0.75)
  r2_p <- c(0.4, 0.8)
  booster_uptake <- c(0.3, 0.75)

  for (j in 1:5) {

    n_diag_rec <- j

    # check vaccination targeting
    y5e <- run_onevax_xpvwrh(tt, gp, vea = 0.5, dur_v = 1, strategy = "VaH",
                             r1 = r1, r2 = r2, r2_p = r2_p,
                             booster_uptake = booster_uptake,
                             n_diag_rec = n_diag_rec)

    for (i in seq_along(y5e)) {
      ## Individuals with no diagnosis history who are treated in X, P, or W
      ## are offered vaccination

      expect_equal(y5e[[i]]$cum_offered[, , c(1, n_diag_rec + 1,
                                              3 * n_diag_rec +  1)],
                   y5e[[i]]$cum_treated[, , c(1, n_diag_rec + 1,
                                              3 * n_diag_rec +  1)])


      ## Individuals with diagnosis history who are are treated or screened in
      ## X, P, or W are offered vaccination
      if (n_diag_rec > 1) {

        expect_equal(y5e[[i]]$cum_offered[, , c(2:n_diag_rec,
                                                n_diag_rec + 2:n_diag_rec,
                                                3 * n_diag_rec +
                                                  (2:n_diag_rec))],
                     y5e[[i]]$cum_treated[, ,  c(2:n_diag_rec,
                                                 n_diag_rec + 2:n_diag_rec,
                                                 3 * n_diag_rec +
                                                   (2:n_diag_rec))] +
                       y5e[[i]]$cum_screened[, , c(2:n_diag_rec,
                                                   n_diag_rec + 2:n_diag_rec,
                                                   3 * n_diag_rec +
                                                     (2:n_diag_rec))])
      }

      # and no-one else
      expect_equal(sum(y5e[[i]]$cum_offered[, , -c(1:n_diag_rec,
                                                   n_diag_rec + 1:n_diag_rec,
                                                   3 * n_diag_rec +
                                                     (1:n_diag_rec))]), 0)

      # uptake % of offered are vaccinated
      expect_equal(y5e[[i]]$cum_offered[, , 1:n_diag_rec] * r1[i],
                   y5e[[i]]$cum_vaccinated[, , 1:n_diag_rec])

      # uptake % of 1st dose offered 2nd dose are vaccinated
      expect_equal(rowSums(y5e[[i]]$cum_offered[, , n_diag_rec + 1:n_diag_rec]
                           * r2_p[i]),
                   rowSums(y5e[[i]]$cum_vaccinated[, , n_diag_rec +
                                                     1:n_diag_rec]))

      # uptake % of offered booster are vaccinated
      expect_equal(y5e[[i]]$cum_offered[, , 3 * n_diag_rec + 1:n_diag_rec]
                   * booster_uptake[i],
                   y5e[[i]]$cum_vaccinated[, , 3 * n_diag_rec + 1:n_diag_rec])

      # and no-one else
      expect_equal(sum(y5e[[i]]$cum_vaccinated[, , -c(1:n_diag_rec,
                                                      n_diag_rec + 1:n_diag_rec,
                                                      3 * n_diag_rec +
                                                        (1:n_diag_rec))]), 0)
      # no-one is lost
      expect_equal(apply(y5e[[i]]$N, 1, sum), rep(6e5, 6), tolerance = 1e-5)
    }
  }
})


test_that("Test length of uptake vector must be 1 or length gp", {

  tt <- seq(0, 5)
  gp <- gono_params(1:2)
  r1 <- c(0.25, 0.5)
  r2 <- c(0.5, 0.75)
  r2_p <- c(0.4, 0.8)
  booster_uptake <- c(0.3, 0.75)

  for (j in 1:5) {

    n_diag_rec <- j

    # check length of uptake vector must be 1 or length(gp)
    expect_error(run_onevax_xpvwrh(tt, gp, vea = 1, dur_v = 1e3,
                                   strategy = "VbE", r1 = c(0, 0.5, 1),
                                   n_diag_rec = n_diag_rec))
    expect_error(run_onevax_xpvwrh(tt, gp, vea = 1, dur_v = 1e3,
                                   strategy = "VbE", r1 = 1,
                                   r2 = c(0, 0.5, 1),
                                   n_diag_rec = n_diag_rec))
    expect_error(run_onevax_xpvwrh(tt, gp, vea = 1, dur_v = 1e3,
                                   strategy = "VbE", r1 = 1, r2 = 1,
                                   r2_p = c(0, 0.5, 1),
                                   n_diag_rec = n_diag_rec))
    expect_error(run_onevax_xpvwrh(tt, gp, vea = 1, dur_v = 1e3,
                                   strategy = "VbE", r1 = 1, r2 = 1,
                                   booster_uptake = c(0, 0.5, 1),
                                   n_diag_rec = n_diag_rec))

    # check length of vea must be 1
    expect_error(run_onevax_xpvwrh(tt, gp, vea = c(0, 1, 1), dur_v = 1e3,
                                   strategy = "VbE", r1 = 1,
                                   n_diag_rec = n_diag_rec))

    # check length of dur_v must be 1
    expect_error(run_onevax_xpvwrh(tt, gp, vea = 1, dur_v = c(0, 1e2, 1e3),
                                   strategy = "VbE", r1 = 1,
                                   n_diag_rec = n_diag_rec))

    # check length of dur_p must be 1
    expect_error(run_onevax_xpvwrh(tt, gp, vea = 1, dur_v = 1e3,
                                   dur_p = c(0, 1e2, 1e3),
                                   strategy = "VbE", r1 = 1,
                                   n_diag_rec = n_diag_rec))
  }
})


test_that("Test revaccination is working", {

  tt <- seq(0, 5)
  gp <- gono_params(1:2)
  r1 <- c(0.25, 0.5)
  r2 <- c(0.5, 0.75)
  r2_p <- c(0.4, 0.8)
  booster_uptake <- c(0.3, 0.75)

  for (j in 1:5) {

    n_diag_rec <- j

    ## test revax is working
    y6e <- run_onevax_xpvwrh(tt, gp, vea = 1, dur_v = 1e-10, dur_revax = 1e10,
                             strategy = "VoD(L)+VoA(H)", r1 = r1, r2 = r2,
                             r2_p = r2_p, booster_uptake = booster_uptake,
                             n_diag_rec = n_diag_rec)

    for (i in seq_along(y6e)) {
      ## L who are treated in X, P, or W are offered vaccination
      expect_equal(y6e[[i]]$cum_offered[, 1, c(1:n_diag_rec,
                                               n_diag_rec + 1:n_diag_rec,
                                               3 * n_diag_rec +
                                                 (1:n_diag_rec))],
                   y6e[[i]]$cum_treated[, 1, c(1:n_diag_rec,
                                               n_diag_rec + 1:n_diag_rec,
                                               3 * n_diag_rec +
                                                 (1:n_diag_rec))])
      ## H who are treated or screened in X, P, or W are offered vaccination
      expect_equal(y6e[[i]]$cum_offered[, 2, c(1:n_diag_rec,
                                               n_diag_rec + 1:n_diag_rec,
                                               3 * n_diag_rec +
                                                 (1:n_diag_rec))],
                   y6e[[i]]$cum_treated[, 2, c(1:n_diag_rec,
                                               n_diag_rec + 1:n_diag_rec,
                                               3 * n_diag_rec +
                                                 (1:n_diag_rec))] +
                     y6e[[i]]$cum_screened[, 2, c(1:n_diag_rec,
                                                  n_diag_rec + 1:n_diag_rec,
                                                  3 * n_diag_rec +
                                                    (1:n_diag_rec))])
      # and no-one else
      expect_equal(sum(y6e[[i]]$cum_offered[, , -c(1:n_diag_rec,
                                                   n_diag_rec + 1:n_diag_rec,
                                                   3 * n_diag_rec +
                                                     (1:n_diag_rec))]), 0)

      # uptake % of offered are vaccinated
      expect_equal(y6e[[i]]$cum_offered[, , 1:n_diag_rec] * r1[i],
                   y6e[[i]]$cum_vaccinated[, , 1:n_diag_rec])
      expect_equal(y6e[[i]]$cum_offered[, , n_diag_rec + 1:n_diag_rec]
                   * r2_p[i],
                   y6e[[i]]$cum_vaccinated[, , n_diag_rec + 1:n_diag_rec])
      expect_equal(y6e[[i]]$cum_offered[, , 3 * n_diag_rec + 1:n_diag_rec]
                   * booster_uptake[i],
                   y6e[[i]]$cum_vaccinated[, , 3 * n_diag_rec + 1:n_diag_rec])

      # stratum V empties immediately
      expect_equal(sum(y6e[[i]]$N[, , 2 * n_diag_rec + 1:n_diag_rec]), 0,
                   tolerance = 1e-5)

      # efficacy is perfect
      expect_equal(sum(y6e[[i]]$cum_incid[, , c(n_diag_rec + 1:n_diag_rec,
                                                2 * n_diag_rec + 1:n_diag_rec,
                                                4 * n_diag_rec +
                                                  1:n_diag_rec)]), 0)
      expect_true(all(y6e[[i]]$cum_incid[-1, , c(1:n_diag_rec,
                                                 3 * n_diag_rec +
                                                   1:n_diag_rec)] > 0))
    }

    y7e <- run_onevax_xpvwrh(tt, gp, vea = 0, vea_revax = 1, dur_v = 1,
                             strategy = "VoD(L)+VoA(H)", r1 = r1, r2 = r2,
                             r2_p = r2_p, booster_uptake = booster_uptake,
                             hes = 0.3, n_diag_rec = n_diag_rec)

    for (i in seq_along(y7e)) {
      ## L who are treated in X, P, or W are offered vaccination
      expect_equal(y7e[[i]]$cum_offered[, 1, c(1:n_diag_rec,
                                               n_diag_rec + 1:n_diag_rec,
                                               3 * n_diag_rec +
                                                 (1:n_diag_rec))],
                   y7e[[i]]$cum_treated[, 1, c(1:n_diag_rec,
                                               n_diag_rec + 1:n_diag_rec,
                                               3 * n_diag_rec +
                                                 (1:n_diag_rec))])
      ## H who are treated or screened in X, P or W are offered vaccination
      expect_equal(y7e[[i]]$cum_offered[, 2, c(1:n_diag_rec,
                                               n_diag_rec + 1:n_diag_rec,
                                               3 * n_diag_rec
                                               + (1:n_diag_rec))],
                   y7e[[i]]$cum_treated[, 2,  c(1:n_diag_rec,
                                                n_diag_rec + 1:n_diag_rec,
                                                3 * n_diag_rec +
                                                  (1:n_diag_rec))] +
                     y7e[[i]]$cum_screened[, 2,  c(1:n_diag_rec,
                                                   n_diag_rec + 1:n_diag_rec,
                                                   3 * n_diag_rec +
                                                     (1:n_diag_rec))])
      # and no-one else
      expect_equal(sum(y7e[[i]]$cum_offered[, , -c(1:n_diag_rec,
                                                   n_diag_rec + 1:n_diag_rec,
                                                   3 * n_diag_rec +
                                                     (1:n_diag_rec))]), 0)

      # uptake % of offered are vaccinated
      expect_equal(y7e[[i]]$cum_offered[, , 1:n_diag_rec] * r1[i],
                   y7e[[i]]$cum_vaccinated[, , 1:n_diag_rec])
      expect_equal(y7e[[i]]$cum_offered[, , n_diag_rec + 1:n_diag_rec]
                   * r2_p[i],
                   y7e[[i]]$cum_vaccinated[, , n_diag_rec + 1:n_diag_rec])
      expect_equal(y7e[[i]]$cum_offered[, , 3 * n_diag_rec + 1:n_diag_rec]
                   * booster_uptake[i],
                   y7e[[i]]$cum_vaccinated[, , 3 * n_diag_rec + 1:n_diag_rec])

      # efficacy is perfect in R
      expect_equal(sum(y7e[[i]]$cum_incid[, , 4 * n_diag_rec + 1:n_diag_rec]),
                   0)
      # but not in XPWVH
      expect_true(all(y7e[[i]]$cum_incid[-1, , -c(n_diag_rec + 1:n_diag_rec,
                                                  4 * n_diag_rec
                                                  + 1:n_diag_rec)] > 0))
    }
  }
})


test_that("Test restart with hesistancy is working", {

  ## test restart with hesitancy is working
  tt <- seq(0, 5)
  gp <- gono_params(1:2)

  for (j in 1:5) {

    n_diag_rec <- j
    y8 <- run_onevax_xpvwrh(tt, gp, vea = 0, dur_v = 1e3,
                            n_diag_rec = n_diag_rec)

    i_p <- lapply(y8, restart_hes, n_vax = 6 * n_diag_rec, hes = 0.5,
                  branching = TRUE, n_erlang = 1, n_diag_rec = n_diag_rec)
    y_hesres <- run_onevax_xpvwrh(tt, gp, init_params = i_p, vea = 0,
                                  dur_v = 1e3, hes = 0.5,
                                  n_diag_rec = n_diag_rec)

    n_erlang <- 2
    y_erlang <- run_onevax_xpvwrh(tt, gp, vea = 0, dur_v = 1e3,
                                  n_erlang = n_erlang, n_diag_rec = n_diag_rec)
    i_p <- lapply(y_erlang, restart_hes,
                  n_vax = 3 * n_diag_rec + 3 * n_diag_rec * n_erlang,
                  hes = 0.5, branching = TRUE, n_erlang = n_erlang,
                  n_diag_rec = n_diag_rec)
    y_hesres_erlang <- run_onevax_xpvwrh(tt, gp, init_params = i_p, vea = 0,
                                         dur_v = 1e3, hes = 0.5,
                                         n_erlang = n_erlang,
                                         n_diag_rec = n_diag_rec)

    # final timepoint y8 run = 2 * first timepoint of y_hesres run(as hes = 0.5)
    # (for XVWR)
    for (i in seq_along(y8)) {
      expect_equal(y8[[i]]$U[length(tt), , 1:5], y_hesres[[i]]$U[1, , 1:5] * 2)
      expect_equal(y8[[i]]$I[length(tt), , 1:5], y_hesres[[i]]$I[1, , 1:5] * 2)
      expect_equal(y8[[i]]$A[length(tt), , 1:5], y_hesres[[i]]$A[1, , 1:5] * 2)
      expect_equal(y8[[i]]$S[length(tt), , 1:5], y_hesres[[i]]$S[1, , 1:5] * 2)
      expect_equal(y8[[i]]$T[length(tt), , 1:5], y_hesres[[i]]$T[1, , 1:5] * 2)
    }

    # check this works for n_erlang > 1
    # final timepoint y_erlang run = 2 * first timepoint of y_hesres_erlang
    # run (as hes = 0.5)
    idx <- stratum_index_xpvwrh(n_erlang = n_erlang, n_diag_rec = n_diag_rec)
    for (i in seq_along(y_erlang)) {
      expect_equal(y_erlang[[i]]$U[length(tt), , 1:(idx$n_vax - n_diag_rec)],
                   y_hesres_erlang[[i]]$U[1, , 1:(idx$n_vax - n_diag_rec)] * 2)
      expect_equal(y_erlang[[i]]$I[length(tt), , 1:(idx$n_vax - n_diag_rec)],
                   y_hesres_erlang[[i]]$I[1, , 1:(idx$n_vax - n_diag_rec)] * 2)
      expect_equal(y_erlang[[i]]$A[length(tt), , 1:(idx$n_vax - n_diag_rec)],
                   y_hesres_erlang[[i]]$A[1, , 1:(idx$n_vax - n_diag_rec)] * 2)
      expect_equal(y_erlang[[i]]$S[length(tt), , 1:(idx$n_vax - n_diag_rec)],
                   y_hesres_erlang[[i]]$S[1, , 1:(idx$n_vax - n_diag_rec)] * 2)
      expect_equal(y_erlang[[i]]$T[length(tt), , 1:(idx$n_vax - n_diag_rec)],
                   y_hesres_erlang[[i]]$T[1, , 1:(idx$n_vax - n_diag_rec)] * 2)
    }


    # restart_hes moves X -> H only, and correctly
    for (i in seq_along(y_hesres)) {
      for (d in 1:n_diag_rec) {
        expect_equal(y_hesres[[i]]$U[, , d],
                     y_hesres[[i]]$U[, , 5 * n_diag_rec + d])
        expect_equal(y_hesres[[i]]$I[, , d],
                     y_hesres[[i]]$I[, , 5 * n_diag_rec + d])
        expect_equal(y_hesres[[i]]$A[, , d],
                     y_hesres[[i]]$A[, , 5 * n_diag_rec + d])
        expect_equal(y_hesres[[i]]$S[, , d],
                     y_hesres[[i]]$S[, , 5 * n_diag_rec + d])
        expect_equal(y_hesres[[i]]$T[, , d],
                     y_hesres[[i]]$T[, , 5 * n_diag_rec + d])
      }

      expect_equal(rowSums(y8[[i]]$U[, , ((n_diag_rec + 1):5 * n_diag_rec)]),
                   rep(0, length(tt)))
      expect_equal(rowSums(y8[[i]]$I[, , ((n_diag_rec + 1):5 * n_diag_rec)]),
                   rep(0, length(tt)))
      expect_equal(rowSums(y8[[i]]$A[, , ((n_diag_rec + 1):5 * n_diag_rec)]),
                   rep(0, length(tt)))
      expect_equal(rowSums(y8[[i]]$S[, , ((n_diag_rec + 1):5 * n_diag_rec)]),
                   rep(0, length(tt)))
      expect_equal(rowSums(y8[[i]]$T[, , ((n_diag_rec + 1):5 * n_diag_rec)]),
                   rep(0, length(tt)))
    }

    # the same for n_erlang > 1
    for (i in seq_along(y_hesres)) {

      for (d in 1:n_diag_rec) {
        expect_equal(y_hesres_erlang[[i]]$U[, , d],
                     y_hesres_erlang[[i]]$U[, , (idx$n_vax - n_diag_rec)  + d])
        expect_equal(y_hesres_erlang[[i]]$I[, , d],
                     y_hesres_erlang[[i]]$I[, , (idx$n_vax - n_diag_rec) + d])
        expect_equal(y_hesres_erlang[[i]]$A[, , d],
                     y_hesres_erlang[[i]]$A[, , (idx$n_vax - n_diag_rec) + d])
        expect_equal(y_hesres_erlang[[i]]$S[, , d],
                     y_hesres_erlang[[i]]$S[, , (idx$n_vax - n_diag_rec) + d])
        expect_equal(y_hesres_erlang[[i]]$T[, , d],
                     y_hesres_erlang[[i]]$T[, , (idx$n_vax - n_diag_rec) + d])
      }

      expect_equal(rowSums(y_hesres_erlang[[i]]$U[, , (n_diag_rec + 1):
                                                    (idx$n_vax - n_diag_rec)]),
                   rep(0, length(tt)))
      expect_equal(rowSums(y_hesres_erlang[[i]]$I[, , (n_diag_rec + 1):
                                                    (idx$n_vax - n_diag_rec)]),
                   rep(0, length(tt)))
      expect_equal(rowSums(y_hesres_erlang[[i]]$A[, , (n_diag_rec + 1):
                                                    (idx$n_vax - n_diag_rec)]),
                   rep(0, length(tt)))
      expect_equal(rowSums(y_hesres_erlang[[i]]$S[, , (n_diag_rec + 1):
                                                    (idx$n_vax - n_diag_rec)]),
                   rep(0, length(tt)))
      expect_equal(rowSums(y_hesres_erlang[[i]]$T[, , (n_diag_rec + 1):
                                                    (idx$n_vax - n_diag_rec)]),
                   rep(0, length(tt)))
    }

    # restart_hes gives error if baseline model run  already has hesitancy

    y9 <- run_onevax_xpvwrh(tt, gp, vea = 0, dur_v = 1e3, hes = 0.5,
                            n_diag_rec = n_diag_rec)

    expect_error(lapply(y9, restart_hes, n_vax = 6, hes = 0.5, branching = TRUE,
                        n_erlang = n_erlang, n_diag_rec = n_diag_rec))

    # restart_hes gives error if baseline model run provided contains vaccinated

    y10 <- run_onevax_xpvwrh(tt, gp, vea = 0, dur_v = 1e3, vbe = 1,
                             n_diag_rec = n_diag_rec)

    expect_error(lapply(y10, restart_hes, n_vax = 6, hes = 0.5,
                        branching = TRUE))

    # check booster_uptake defaults to r1r2 primary uptake
    y_r1r2_only <- run_onevax_xpvwrh(tt, gp, vea = 0, dur_v = 1,
                                     r1 = 0.8, r2 = 0.5,
                                     strategy = "VoD(L)+VoA(H)",
                                     n_diag_rec = n_diag_rec)

    y_r1r2_boost <- run_onevax_xpvwrh(tt, gp, vea = 0, dur_v = 1,
                                      r1 = 0.8, r2 = 0.5,
                                      booster_uptake = 0.4,
                                      strategy = "VoD(L)+VoA(H)",
                                      n_diag_rec = n_diag_rec)

    expect_equal(y_r1r2_only, y_r1r2_boost)

    # When everyone receives one dose (r2 = 0), V, W R, H all empty
    y11 <- run_onevax_xpvwrh(tt, gp, r1 = 1, r2 = 0, strategy = "VoD",
                             n_diag_rec = n_diag_rec)

    for (i in seq_along(y11)) {

      expect_equal(sum(y11[[i]]$N[, , -c(1:n_diag_rec,
                                         (n_diag_rec + 1):(2 * n_diag_rec))]),
                   0)
      expect_true(all(rowSums(y11[[i]]$N[, , 1:n_diag_rec]) > 0))
      expect_true(all(rowSums(y11[[i]]$N[-1, , (n_diag_rec + 1):
                                           (2 * n_diag_rec)]) > 0))
    }

    # When everyone receives two doses straight away (r1*r2 = 1),
    # P and H all empty

    y12 <- run_onevax_xpvwrh(tt, gp, r1 = 1, r2 = 1, strategy = "VoD",
                             n_diag_rec = n_diag_rec)

    for (i in seq_along(y12)) {

      expect_equal(sum(y12[[i]]$N[, , c((n_diag_rec + 1):(2 * n_diag_rec),
                                        (5 * n_diag_rec + 1):
                                          (6 * n_diag_rec))]), 0)
      expect_true(all(rowSums(y12[[i]]$N[, , 1:n_diag_rec]) > 0))
      expect_true(all(rowSums(y12[[i]]$N[-1, , (2 * n_diag_rec + 1):
                                           (5 * n_diag_rec)]) > 0))
    }

    # if half the population gets 1 dose and half the population gets 2 doses
    # and durations are the same
    # then N should be the same for P and V

    y13 <- run_onevax_xpvwrh(tt, gp, r1 = 1, r2 = 0.5, strategy = "VoD",
                             n_diag_rec = n_diag_rec)

    for (i in seq_along(y13)) {
      for (d in 1:n_diag_rec) {
        expect_equal(y13[[i]]$N[, , (n_diag_rec + d)],
                     y13[[i]]$N[, , (2 * n_diag_rec + d)])
      }
    }

    # incidence should also be the same as all vex_p default to vex
    for (i in seq_along(y13)) {
      for (d in 1:n_diag_rec) {
        expect_equal(y13[[i]]$cum_incid[, , (n_diag_rec + d)],
                     y13[[i]]$cum_incid[, , (2 * n_diag_rec + d)])
      }
    }
  }
})

test_that("Test vea_p works independently of vex", {

  tt <- seq(0, 5)
  gp <- gono_params(1:2)

  for (j in 1:5) {

    n_diag_rec <- j

    # test vex_p work independently of vex
    y14 <- run_onevax_xpvwrh(tt, gp, r1 = 1, r2 = 0.5, vea_p = 1,
                             strategy = "VoD", n_diag_rec = n_diag_rec)

    # incidence in V but no incidence in P as vea_p = 1
    # individuals present in both
    for (i in seq_along(y14)) {
      expect_true(all(rowSums(y14[[i]]$N[-1, , (n_diag_rec + 1):
                                           (2 * n_diag_rec)]) > 0))
      expect_true(all(rowSums(y14[[i]]$N[-1, , (2 * n_diag_rec + 1):
                                           (3 * n_diag_rec)]) > 0))
      expect_true(sum(y14[[i]]$cum_incid[-1, , (n_diag_rec + 1):
                                           (2 * n_diag_rec)]) == 0)
      expect_true(all(rowSums(y14[[i]]$cum_incid[-1, , (2 * n_diag_rec + 1):
                                                   (3 * n_diag_rec)]) > 0))
    }

    # test vea_p defaults to vea
    y_vea_only <- run_onevax_xpvwrh(tt, gp, r1 = 1, r2 = 0.5, vea = 0.5,
                                    strategy = "VoD", n_diag_rec = n_diag_rec)
    y_vea_veap <- run_onevax_xpvwrh(tt, gp, r1 = 1, r2 = 0.5, vea = 0.5,
                                    vea_p = 0.5, strategy = "VoD",
                                    n_diag_rec = n_diag_rec)

    expect_equal(y_vea_only, y_vea_veap)

    # when vea, vea_p and vea_revax = 1, no incidence among vacc individuals

    y_no_incid <- run_onevax_xpvwrh(tt, gp, r1 = 1, r2 = 0.5, vea = 1,
                                    n_diag_rec = n_diag_rec)

    for (i in seq_along(y_no_incid)) {
      expect_equal(sum(y_no_incid[[i]]$cum_incid[, , -(1:n_diag_rec)]), 0)
    }

    # check initial params generated correctly when coverage is specified
    pars <- lapply(gp[1], model_params)
    init_1 <- initial_params_xpvwrh(pars[[1]], coverage_p = 0.25,
                                    coverage_v = 0.25, n_diag_rec = n_diag_rec)

    # expect same number starting in P as in V
    expect_equal(init_1$U0[, n_diag_rec + (1:n_diag_rec)],
                 init_1$U0[, 2 * n_diag_rec + (1:n_diag_rec)])

    init_2 <- initial_params_xpvwrh(pars[[1]], coverage_p = 0.5,
                                    coverage_v = 0.25, n_diag_rec = n_diag_rec)

    # expect double starting in P as in V
    expect_equal((init_2$U0[, n_diag_rec + (1:n_diag_rec)]),
                 init_2$U0[, 2 * n_diag_rec + (1:n_diag_rec)] * 2)

    # expect error if sum of coverages and hes is greater than 1
    expect_error(initial_params_xpvwrh(pars[[1]], coverage_p = 0.5,
                                       coverage_v = 0.5, hes = 0.5,
                                       n_diag_rec = n_diag_rec))

    # expect error if sum of coverages is 100%
    # (if total coverage is 100%, no asymptomatic infection is seeded)
    expect_error(initial_params_xpvwrh(pars[[1]], coverage_p = 0.5,
                                       coverage_v = 0.5,
                                       n_diag_rec = n_diag_rec))

    # tests correct number of individuals are moving to V and to P
    # duration is large for v and p so assuming no waning
    y15 <- run_onevax_xpvwrh(tt, gp, r1 = 0.5, r2 = 0.5, dur_v = 1e90,
                             strategy = "VoD", n_diag_rec = n_diag_rec)

    # number getting vaccinated from X (timepoint 2 as all is 0 in timepoint 1)
    cum_vac <- sum(y15[[1]]$cum_vaccinated[2, , 1:n_diag_rec])
    cum_off <- sum(y15[[1]]$cum_offered[2, , 1:n_diag_rec])

    # 50% of individuals offered will accept and be vaccinated (0.25 + 0.25)
    prop_accept <- cum_vac / cum_off
    expect_equal(prop_accept, 0.5)

    # 50% of those treated, are vaccinated
    expect_equal(sum(y15[[1]]$cum_vaccinated[2, , 1:n_diag_rec])
                 / sum(y15[[1]]$cum_treated[2, , 1:n_diag_rec]), 0.5)

    # number total who received vaccination
    # i.e number in P and V at timepoint 2 given no waning
    p <- sum(y15[[1]]$N[2, , n_diag_rec + (1:n_diag_rec)])
    v <- sum(y15[[1]]$N[2, , 2 * n_diag_rec + (1:n_diag_rec)])

    # should be same number in P as in V
    expect_equal(p, v)

    # test individuals in P are actually waning back into X
    # move all individuals into P for starting conditions
    # no vaccination, duration of vaccine short so waning is fast
    # should see individuals accumulating in X, decreasing in P

    pars <- lapply(gp[1], model_params, n_diag_rec = n_diag_rec)
    ip <- lapply(pars, initial_params_xpvwrh, coverage_p = 0.99999999999999,
                 t = 5, n_diag_rec = n_diag_rec)

    y16 <- run_onevax_xpvwrh(tt, gp, init_params = ip, dur_p = 1e-90,
                             n_diag_rec = n_diag_rec)

    # entire population starts in P then wanes to X and stays there
    for (i in seq_along(y16)) {
      expect_equal(sum(y16[[i]]$N[1, , n_diag_rec + (1:n_diag_rec)]), 6e+05)
      expect_equal(rowSums(y16[[i]]$N[2:(max(tt) + 1), , 1:n_diag_rec]),
                   rep(6e+05, max(tt)))

      # V W R are empty for all time
      expect_equal(rowSums(y16[[i]]$N[, , (2 * n_diag_rec + 1):
                                        (6 * n_diag_rec)]), rep(0, 6))
    }

    # test individuals in P move to V
    # dur_v and dur_p are huge so no waning
    # (note cannot start entire population in P as model set up to have asymp
    # infections to seed transmission in X or H only)

    ip <- lapply(pars, initial_params_xpvwrh, coverage_p = 0.9,
                 t = 5, n_diag_rec = n_diag_rec)
    y_p <- run_onevax_xpvwrh(tt, gp, init_params = ip, r2_p = 1, dur_v = 1e100,
                             strategy = "VoD(L)+VoA(H)",
                             n_diag_rec = n_diag_rec)

    # 0.9 * population starts in P
    # no vaccination from X, and others
    # only vaccination from P
    # vaccination into V so stratum > 0 for time > 1

    ## Note - waning being initialised with numerical error val? (~1e-99?)
    expect_equal(sum(y_p[[i]]$N[1, , n_diag_rec + (1:n_diag_rec)]),
                 6e+05 * 0.9)
    expect_equal(sum(y_p[[i]]$cum_vaccinated[, , -((n_diag_rec + 1):
                                                     (2 * n_diag_rec))]), 0)
    expect_true(unique(rowSums(y_p[[i]]$cum_vaccinated[-1, ,
                                                       ((n_diag_rec + 1):
                                                          (2 * n_diag_rec))])
                       > 0))

    expect_true(unique(rowSums(y_p[[i]]$N[-1, ,
                                          (2 * n_diag_rec + 1):
                                            (3 * n_diag_rec)]) > 0))

    # test individuals still wane to W
    ip <- lapply(pars, initial_params_xpvwrh, coverage_v = 0.99999999999999
                 , t = 5, n_diag_rec = n_diag_rec)
    y17 <- run_onevax_xpvwrh(tt, gp, init_params = ip, dur_v = 1e-90,
                             n_diag_rec = n_diag_rec)

    # entire population starts in V then wanes to W and stays there
    # note, W won't equal 6e+05 as 1. individuals in W die but 2. entrants
    # enter into X 3. no vaccination so replenishment of V and subsequently W

    for (i in seq_along(y17)) {
      expect_equal(sum(y17[[i]]$N[1, , (2 * n_diag_rec + 1):(3 * n_diag_rec)]),
                   6e+05)
      expect_equal(rowSums(y17[[i]]$N[2:(max(tt) + 1), ,
                                      (2 * n_diag_rec + 1):(3 * n_diag_rec)]),
                   rep(0, max(tt)))

      for (j in 2:max(tt) + 1) {
        expect_true(sum(y17[[i]]$N[j, , (3 * n_diag_rec + 1):
                                     (4 * n_diag_rec)]) > 0)
      }

      # R, H are empty for all time
      expect_equal(rowSums(y17[[i]]$N[, , (5 * n_diag_rec + 1):
                                        (6 * n_diag_rec)]), rep(0, 6))
    }
  }
})

test_that("run_onevax_xpvwrh works when n_erlang > 1", {

  for (j in 1:5) {

    n_diag_rec <- j
    tt <- seq(0, 5)
    pars <- lapply(gono_params()[1], model_params, n_diag_rec = n_diag_rec)

    n_erlang_2 <- initial_params_xpvwrh(pars[[1]], n_erlang = 2,
                                        n_diag_rec = n_diag_rec)
    expect_equal(dim(n_erlang_2[[1]]), c(2, 9 * n_diag_rec))

    n_erlang_3 <- initial_params_xpvwrh(pars[[1]], n_erlang = 3,
                                        n_diag_rec = n_diag_rec)
    expect_equal(dim(n_erlang_3[[1]]), c(2, 12 * n_diag_rec))

    n_erlang_10 <- initial_params_xpvwrh(pars[[1]], n_erlang = 10,
                                         n_diag_rec = n_diag_rec)
    expect_equal(dim(n_erlang_10[[1]]), c(2, 33 * n_diag_rec))

    # correct number of erlang compartments in run function
    y <- run_onevax_xpvwrh(tt, gono_params = gono_params()[1], vea = 0,
                           dur_v = 1e3, n_erlang = 3, r1 = 0.5,
                           r2 = 0.1, booster_uptake = 1,
                           strategy = "VoD(L)+VoA(H)",
                           n_diag_rec = n_diag_rec)

    expect_true(dim(y[[1]]$N)[3] == 12 * n_diag_rec)

    # everything defaults to the n_vax = xpvwrh format when n_erlang = 1 or
    # not supplied
    y <- run_onevax_xpvwrh(tt, gono_params = gono_params()[1], vea = 0,
                           dur_v = 1e3, r1 = 0.5, r2 = 0.1, booster_uptake = 1,
                           strategy = "VoD(L)+VoA(H)", n_diag_rec = n_diag_rec)
    expect_true(dim(y[[1]]$N)[3] == 6 * n_diag_rec)

    # expect error if n_erlang supplied to run function doesn't match initial
    # conditions
    expect_error(run_onevax_xpvwrh(tt, gono_params = gono_params()[1], vea = 0,
                                   dur_v = 0.1, booster_uptake = 1,
                                   strategy = "VoD(L)+VoA(H)",
                                   init_params = list(n_erlang_2),
                                   n_erlang = 3, n_diag_rec = n_diag_rec))

    # all erlang compartments are protected
    # vea = 1 no on gets infected
    gp <- gono_params()[1]

    n_erlang <- 2
    idx <- stratum_index_xpvwrh(n_erlang, n_diag_rec = n_diag_rec)
    y <- run_onevax_xpvwrh(tt, gono_params = gp, vea = 1, vea_p = 1,
                           vea_revax = 1,
                           dur_v = 1e3, r1 = 0.5,
                           r2 = 0.1, booster_uptake = 1,
                           strategy = "VoD(L)+VoA(H)",
                           n_erlang = n_erlang,
                           n_diag_rec = n_diag_rec)

    # people are vaccinated (indiviudals exist in P1-P2, V1-V2, R1-R2)
    for (i in 1: (2 * n_diag_rec + 3 * n_diag_rec * n_erlang)) {
      expect_true(all(y[[1]]$N[-1, , c(i)] > 0))
    }

    # and no one vaccinated is infected
    expect_true(all(y[[1]]$cum_incid[-1, , (n_diag_rec + 1):
                                       (n_diag_rec + (n_diag_rec * n_erlang))]
                    == 0))  #P

    expect_true(all(y[[1]]$cum_incid[-1, ,
                                     (n_diag_rec + (n_diag_rec * n_erlang) + 1)
                                     :(n_diag_rec + 2 *
                                         (n_diag_rec * n_erlang))] == 0))  #V

    expect_true(all(y[[1]]$cum_incid[-1, , (2 * n_diag_rec + 2 *
                                              (n_diag_rec * n_erlang) + 1)
                                     :(2 * n_diag_rec + 3 *
                                         (n_diag_rec * n_erlang))] == 0))  #R

    expect_true(all(y[[1]]$cum_incid[-1, , c(idx$P, idx$V, idx$R)] == 0))


    # and non-vaccinated do get infected
    expect_true(all(y[[1]]$cum_incid[-1, , 1:n_diag_rec] > 0)) #X
    expect_true(all(y[[1]]$cum_incid[-1, , (n_diag_rec + 2 *
                                              (n_diag_rec * n_erlang) + 1)
                                     :(2 * n_diag_rec + 2 *
                                         (n_diag_rec * n_erlang))] > 0)) #W

    expect_true(all(y[[1]]$cum_incid[-1, , c(idx$X, idx$W)] > 0))

    # expect vaccination maps to be correctly generated, moving people
    # to the correct strata

    pars <- lapply(gp[1], model_params, n_diag_rec = n_diag_rec)
    vbe <- 1
    p <- set_strategy(strategy = "VoD(L)+VoA(H)", vbe > 0)

    n_erlang <- 2
    idx <- stratum_index_xpvwrh(n_erlang, n_diag_rec = n_diag_rec)

    i_eligible <- c(1:n_diag_rec, 1:n_diag_rec,
                    (n_diag_rec + 1):(n_diag_rec + n_diag_rec * n_erlang),
                    (n_diag_rec + 2 * n_diag_rec * n_erlang + 1)
                    :(2 * n_diag_rec + 2 * n_diag_rec * n_erlang))
    expect_equal(i_eligible, c(idx$X, idx$X, idx$P, idx$W))

    i_p <- c((n_diag_rec + 1) : (n_diag_rec  + n_diag_rec),
             rep((n_diag_rec + n_diag_rec * n_erlang + 1):
                   (n_diag_rec + n_diag_rec * n_erlang + n_diag_rec),
                 n_erlang + 1), ((2 * n_diag_rec + 2 *
                                    n_diag_rec * n_erlang + 1)
                                 :(2 * n_diag_rec + 2 * n_diag_rec *
                                     n_erlang + n_diag_rec)))
    expect_equal(i_p, c(idx$P[1:n_diag_rec], rep(idx$V[1:n_diag_rec],
                                                 n_erlang + 1),
                        idx$R[1:n_diag_rec]))

    i_p_temp_2 <- i_p
    i_p_temp_2[i_p %% n_diag_rec != 0] <- i_p[i_p %% n_diag_rec != 0] + 1

    vod_map <- create_vax_map_branching(n_vax = n_vax, p$vod, i_eligible,
                                        i_p_temp_2, idx = idx)
    vos_map <- create_vax_map_branching(n_vax = n_vax, p$vos, i_eligible,
                                        i_p, idx = idx)
    vbe_map <- create_vax_map_branching(n_vax = n_vax, p$vbe, i_eligible,
                                        i_p, set_vbe = TRUE, idx = idx)

    # for vod, expect:
    # all 1's indicating vaccination, are in the right place
    # and all -1's indicating where people are being vaccinated to are
    # in the right place

    for (d in 1:n_diag_rec) {
      expect_true(unique(vod_map[, idx$X[d], idx$X[d]]  == c(1, 1)))
    }

    for (d in 1:n_diag_rec) {
      if (d < n_diag_rec) {
        expect_true(unique(vod_map[, idx$P[d + 1], idx$X[d]]  == c(-1, -1)))
      } else {
        expect_true(unique(vod_map[, idx$P[d], idx$X[d]]  == c(-1, -1)))
      }
    }

    for (e in 1:n_erlang) {
      for (d in 1:n_diag_rec) {
        expect_true(unique(vod_map[, idx$P[(e - 1) * n_diag_rec + d],
                                   idx$P[(e - 1) * n_diag_rec + d]] == c(1, 1)))
      }
    }

    for (e in 1:n_erlang) {
      for (d in 1:n_diag_rec) {
        if (d < n_diag_rec) {
          expect_true(unique(vod_map[, idx$V[d + 1], idx$P[(e - 1) *
                                                             n_diag_rec + d]]
                             == c(-1, -1)))
        } else {
          expect_true(unique(vod_map[, idx$V[d], idx$P[(e - 1) *
                                                         n_diag_rec + d]]
                             == c(-1, -1)))
        }
      }
    }

    for (d in 1:n_diag_rec) {
      expect_true(unique(vod_map[, idx$W[d], idx$W[d]] == c(1, 1)))
    }

    for (d in 1:n_diag_rec) {
      if (d < n_diag_rec) {
        expect_true(unique(vod_map[, idx$R[d + 1], idx$W[d]] == c(-1, -1)))
      } else {
        expect_true(unique(vod_map[, idx$R[d], idx$W[d]] == c(-1, -1)))
      }
    }

    # everything else is 0 in those i_eligible strata

    for (d in 1:n_diag_rec) {
      if (d < n_diag_rec) {
        expect_equal(sum(vod_map[, -c(idx$X[d], idx$P[d + 1], idx$V[d + 1]),
                                 idx$X[d]]), 0)
      } else {
        expect_equal(sum(vod_map[, -c(idx$X[d], idx$P[d], idx$V[d]),
                                 idx$X[d]]), 0)
      }
    }

    for (e in 1:n_erlang) {
      for (d in 1:n_diag_rec) {
        if (d < n_diag_rec) {
          expect_equal(sum(abs(vod_map[, -c(idx$P[(e - 1) * n_diag_rec + d],
                                            idx$V[d + 1]),
                                       idx$P[(e - 1) * n_diag_rec + d]])), 0)
        } else {
          expect_equal(sum(abs(vod_map[, -c(idx$P[(e - 1) * n_diag_rec + d],
                                            idx$V[d]),
                                       idx$P[(e - 1) * n_diag_rec + d]])), 0)
        }
      }
    }

    for (d in 1:n_diag_rec) {
      if (d < n_diag_rec) {
        expect_equal(sum(vod_map[, -c(idx$W[d],  idx$R[d + 1]), idx$W[d]]), 0)
      } else {
        expect_equal(sum(vod_map[, -c(idx$W[d],  idx$R[d]), idx$W[d]]), 0)
      }
    }

    expect_equal(sum(vod_map[, , -c(i_eligible)]), 0)

    # for vos, expect:
    # all 1's indicating vaccination, are in the right place
    # and all -1's indicating where people are being vaccinated to are
    # in the right place

    for (d in 1:n_diag_rec) {
      expect_true(unique(vos_map[, idx$X[d], idx$X[d]] == c(0, 1)))
      expect_true(unique(vos_map[, idx$P[d], idx$X[d]] == c(0, -1)))
      expect_true(unique(vos_map[, idx$V[d], idx$X[d]] == c(0, -1)))
      expect_true(unique(vos_map[, idx$W[d], idx$W[d]] == c(0, 1)))
      expect_true(unique(vos_map[, idx$R[d], idx$W[d]] == c(0, -1)))
    }

    for (e in 1:n_erlang) {
      for (d in 1:n_diag_rec) {
        expect_true(unique(vos_map[, idx$P[(e - 1) * n_diag_rec + d],
                                   idx$P[(e - 1) * n_diag_rec + d]]
                           == c(0, 1)))
        expect_true(unique(vos_map[, idx$V[d], idx$P[(e - 1) * n_diag_rec + d]]
                           == c(0, -1)))
      }
    }

    # everything else is 0 in those i_eligible strata
    for (d in 1:n_diag_rec) {
      expect_equal(sum(vos_map[, -c(idx$X[d], idx$P[d], idx$V[d]), idx$X[d]]),
                   0)
      expect_equal(sum(abs(vos_map[, -c(idx$W[d], idx$R[d]), idx$W[d]])), 0)

    }

    for (e in 1:n_erlang) {
      for (d in 1:n_diag_rec) {
        expect_equal(sum(abs(vos_map[, -c(idx$P[(e - 1) * n_diag_rec + d],
                                          idx$V[d]),
                                     idx$P[(e - 1) * n_diag_rec + d]])), 0)
      }
    }

    # everything is entirely 0 when not in i_eligible strata
    expect_equal(sum(vos_map[, , -c(i_eligible)]), 0)

    # for vbe, expect:
    expect_true(unique(vbe_map[, idx$X[1], idx$X[1]] == c(1, 1)))
    expect_true(unique(vbe_map[, idx$V[1], idx$X[1]] == c(-1, -1)))
    expect_equal(sum(vbe_map[, -c(idx$X[1], idx$V[1]), idx$X[1]]), 0)
    expect_equal(sum(vbe_map[, , 2:(3 * n_diag_rec + 3 * n_diag_rec *
                                      n_erlang)]), 0)

    # test uptake maps are generated as expected
    r1 <- c(0.25, 0.5)
    r2 <- c(0.5, 0.75)
    r2_p <- c(0.4, 0.8)
    booster_uptake <- c(0.3, 0.75)

    for (i in seq_along(r1)) {
      u <- create_uptake_map_xpvwrh(vod_map, r1[i], r2[i], r2_p[i],
                                    booster_uptake[i],
                                    idx = idx, n_diag_rec = n_diag_rec,
                                    n_erlang = n_erlang,
                                    screening_or_diagnosis = "diagnosis")

      # vaccination mapping multiplied by corresponding uptake % mapping
      acc_vax <- u * vod_map

      # uptakes from X are at the expected values
      # uptakes from P1 and P2 to V are at the expected values
      # uptakes from W to R1 are at the expected values


      for (d in 1:n_diag_rec) {
        expect_true(unique(acc_vax[, idx$X[d], idx$X[d]] == c(r1[i], r1[i])))

        expect_true(unique(acc_vax[, idx$W[d], idx$W[d]] ==
                             c(booster_uptake[i], booster_uptake[i])))

        if (d < n_diag_rec) {
          expect_true(unique(acc_vax[, idx$P[d +  1],
                                     idx$X[d]] == c(- (r1[i] * (1 - r2[i])),
                                                    - (r1[i] * (1 - r2[i])))))

          expect_true(unique(acc_vax[, idx$R[d + 1], idx$W[d]] ==
                               -c(booster_uptake[i],
                                  booster_uptake[i])))
        } else {
          expect_true(unique(acc_vax[, idx$P[d],
                                     idx$X[d]] == c(- (r1[i] * (1 - r2[i])),
                                                    - (r1[i] * (1 - r2[i])))))

          expect_true(unique(acc_vax[, idx$V[d],
                                     idx$X[d]] == c(-r1[i] * r2[i],
                                                    -r1[i] * r2[i])))

          expect_true(unique(acc_vax[, idx$R[d], idx$W[d]] ==
                               -c(booster_uptake[i], booster_uptake[i])))
        }
      }


      for (e in 1:n_erlang) {
        for (d in 1:n_diag_rec) {
          expect_true(unique(acc_vax[, idx$P[(e - 1) * n_diag_rec + d],
                                     idx$P[(e - 1) * n_diag_rec + d]] ==
                               c(r2_p[i], r2_p[i])))

          if (d < n_diag_rec) {
            expect_true(unique(acc_vax[, idx$V[d + 1], idx$P[(e - 1) *
                                                               n_diag_rec + d]]
                               == -c(r2_p[i], r2_p[i])))
          } else {
            expect_true(unique(acc_vax[, idx$V[d], idx$P[(e - 1) *
                                                           n_diag_rec + d]]
                               == -c(r2_p[i], r2_p[i])))
          }
        }
      }
    }

    # people wane as expected
    # entire population starts in V1 and wanes through V2 to W only
    ip_v <- lapply(pars, initial_params_xpvwrh, coverage_v = 0.99999999999999,
                   n_erlang = 2, n_diag_rec = n_diag_rec, t = FALSE)
    y_v <- run_onevax_xpvwrh(tt, gp, init_params = ip_v, n_erlang = 2,
                             n_diag_rec = n_diag_rec)

    # people in no diagnoses categories of V1, V2, W
    # (no infections so no diagnoses!)
    for (e in 1:n_erlang) {
      expect_true(all(y_v[[1]]$N[-1, , idx$V[(e - 1) * n_diag_rec + 1]] > 0))
    }

    expect_true(all(y_v[[1]]$N[-1, ,  idx$W[1]] > 0))

    # no one in P1, P2, R1, R2 (people born into X so not included here)
    expect_true(all(y_v[[1]]$N[, , c(idx$P, idx$R)] ==  0))

    # people are flowing from V1(4) -> V2(5) -> W(6)
    # would expect the number in V1 to be greater than V2 and W if rate
    # of waning is constant and people are moving from one to the next
    # to the next in time

    for (e in 1:n_erlang) {
      if (e < n_erlang) {
        expect_true(sum(rowSums(y_v[[1]]$N[, , idx$V[(e - 1) *
                                                       n_diag_rec + 1]])) >
                      sum(rowSums(y_v[[1]]$N[, , idx$V[e * n_diag_rec + 1]])))
      } else {
        expect_true(sum(rowSums(y_v[[1]]$N[, , idx$V[(e - 1) *
                                                       n_diag_rec + 1]])) >
                      sum(rowSums(y_v[[1]]$N[, , idx$W[1]])))
      }
    }

    # entire population starts in P1 and wanes through P2 to X only
    # ( + no further vaccination)
    ip_p <- lapply(pars, initial_params_xpvwrh,
                   coverage_p = 0.99999999999999,
                   n_erlang = 2, n_diag_rec = n_diag_rec, t = FALSE)
    y_p <- run_onevax_xpvwrh(tt, gp, init_params = ip_p, n_erlang = 2,
                             n_diag_rec = n_diag_rec)

    # people in no diagnoses categories of X, P1, P2
    # (no infections so no diagnoses!)

    expect_true(all(y_p[[1]]$N[-1, , idx$X[1]] > 0))

    for (e in 1:n_erlang) {
      expect_true(all(y_p[[1]]$N[-1, , idx$P[(e - 1) * n_diag_rec + 1]] > 0))
    }

    # no one in V1, V2, W, R1 or R2
    expect_true(all(y_p[[1]]$N[, , c(idx$V, idx$W, idx$R)] ==  0))

    # people are flowing from P1(2) -> P2(3)

    for (e in 1:(n_erlang - 1)) {
      expect_true(sum(rowSums(y_p[[1]]$N[, , idx$P[e]])) >
                    sum(rowSums(y_p[[1]]$N[, , idx$P[e + 1]])))
    }

    # people also wane from R1 and R2
    n_erlang <- 2
    idx <- stratum_index_xpvwrh(n_erlang, n_diag_rec = n_diag_rec)

    #run for ages with very short duration of primary vaccination,
    # 100% uptake and very long duration of booster vacination to get
    # people to converge in R1

    tt <- seq(0, 100)
    y_long <- run_onevax_xpvwrh(tt, gp, vea = 0, dur_v = 0.1,
                                n_erlang = n_erlang,
                                r1 = 1, r2 = 1, booster_uptake = 1,
                                dur_revax = 1e10, strategy = "VoD(L)+VoA(H)",
                                n_diag_rec = n_diag_rec)

    i_p <- lapply(y_long, restart_hes, n_vax = idx$n_vax, n_erlang = n_erlang,
                  n_diag_rec = n_diag_rec, branching = TRUE)

    tt <- seq(0, 5)
    # run model with most people in W + R1
    y_r <- run_onevax_xpvwrh(tt, gp, init_params = i_p, r1 = 0, r2 = 0,
                             dur_v = 0.1, dur_revax = 10, n_erlang = n_erlang,
                             n_diag_rec = n_diag_rec)

    # when lots of people start revaccinated and (almost) none start
    # primarily vaccinated, people flow from R1 -> R2
    for (e in 1:(n_erlang - 1)) {
      expect_true(sum(rowSums(y_r[[1]]$N[, , idx$R[e]])) >
                    sum(rowSums(y_r[[1]]$N[, , idx$R[e + 1]])))
    }

    # check protection being assigned correctly
    n_erlang_vec <- seq(1, 10)

    for (d in 1:n_diag_rec) {
      for (e in n_erlang_vec) {
        n_erlang <- e
        idx <- stratum_index_xpvwrh(n_erlang = n_erlang,
                                    n_diag_rec = n_diag_rec)
        i_v <- c(idx$P, idx$V, idx$R)

        p_vec <- set_protection(i_v = i_v, idx = idx, n_vax = idx$n_vax,
                                ve_p = 0.5, ve = 0.75, ve_revax = 1)
        expect_equal(p_vec[idx$P], rep(0.5, n_erlang * n_diag_rec))
        expect_equal(p_vec[idx$V], rep(0.75, n_erlang * n_diag_rec))
        expect_equal(p_vec[idx$R], rep(1, n_erlang * n_diag_rec))
        expect_equal(p_vec[c(idx$X, idx$W, idx$H)], rep(0, 3 * n_diag_rec))

      }
    }


    ### waning maps are generated correctly with the correct rates
    n_erlang <- 3
    idx <- stratum_index_xpvwrh(n_erlang = n_erlang, n_diag_rec = n_diag_rec)

    # i_v and i_w generated in the same way as in vax_params function
    i_v <- c(idx$P, idx$V, idx$R)
    i_w <- c(idx$P[-(1:n_diag_rec)], idx$X, idx$V[-(1:n_diag_rec)], idx$W,
             idx$R[-(1:n_diag_rec)], idx$W)

    # set durations
    dur_p <- 1e03
    dur_v <- 2e03
    dur_revax <- 3e03

    w <- create_waning_map_branching(idx$n_vax, i_v, i_w,
                                     n_erlang / c(dur_p, dur_v, dur_revax),
                                     n_erlang = n_erlang,
                                     n_diag_rec = n_diag_rec)
    # width and height of waning map are both equal to n_vax
    expect_true(dim(w)[1] == c(idx$n_vax))
    expect_true(dim(w)[2] == c(idx$n_vax))

    #people wane from P1 -> P2 -> P3 -> X at the same, expected rate
    #people wane from V1 -> V2 -> V3 -> W at the same, expected rate
    #people wane from R1 -> R2 -> R3 -> W at the same, expected rate
    for (e in 1:n_erlang) {
      for (d in 1:n_diag_rec) {

        expect_true(all(w[idx$P[(e - 1) * n_diag_rec + d],
                          idx$P[(e - 1) * n_diag_rec + d]]
                        == -n_erlang / dur_p))
        expect_true(all(w[idx$V[(e - 1) * n_diag_rec + d],
                          idx$V[(e - 1) * n_diag_rec + d]]
                        == -n_erlang / dur_v))
        expect_true(all(w[idx$R[(e - 1) * n_diag_rec + d],
                          idx$R[(e - 1) * n_diag_rec + d]]
                        == -n_erlang / dur_revax))

        if (e < n_erlang) {
          expect_true(all(w[idx$P[(e) * n_diag_rec + d],
                            idx$P[(e - 1) * n_diag_rec + d]]
                          == n_erlang / dur_p))
          expect_true(all(w[idx$V[(e) * n_diag_rec + d],
                            idx$V[(e - 1) * n_diag_rec + d]]
                          == n_erlang / dur_v))
          expect_true(all(w[idx$R[(e) * n_diag_rec + d],
                            idx$R[(e - 1) * n_diag_rec + d]]
                          == n_erlang / dur_revax))

        } else {
          expect_true(all(w[idx$X[d], idx$P[(e - 1) * n_diag_rec + d]] ==
                            n_erlang / dur_p))
          expect_true(all(w[idx$W[d], idx$V[(e - 1) * n_diag_rec + d]] ==
                            n_erlang / dur_v))
          expect_true(all(w[idx$W[d], idx$R[(e - 1) * n_diag_rec + d]] ==
                            n_erlang / dur_revax))
        }
      }
    }
  }
})
