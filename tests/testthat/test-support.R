test_that("aggregate works", {
  dt <- 1 / 12
  tt <- seq(0, 1, dt)
  y <- run_onevax_xvwv(tt, gono_params(1:2), vea = 1, dur = 4, uptake = 1,
                       strategy = "VoD", vbe = 1)
  z <- aggregate(y, "cum_incid", FALSE)
  expect_equal(dim(z), c(2, length(tt)))
  expect_equal(z[1, ], apply(y[[1]]$cum_incid, 1, sum))
  expect_equal(z[2, ], apply(y[[2]]$cum_incid, 1, sum))

  z1 <- aggregate(y, "cum_incid", TRUE)
  expect_equal(z1, t(apply(z, 1, diff)) / dt)

  z2 <- aggregate(y, "cum_incid", FALSE, f = mean)
  expect_equal(z2, colMeans(z))

  z3 <- aggregate(y, "cum_incid", TRUE, f = mean)
  expect_equal(z3, colMeans(z1))

  z4 <- aggregate(y, "cum_incid", TRUE, stratum = 1, f = mean)
  z5 <- aggregate(y, "cum_incid", TRUE, stratum = 3, f = mean)
  expect_equal(z4 + z5, z3)

})

test_that("extract_flows_xpvwrh works", {
  tt <- seq(0, 2)
  y <- run_onevax_xpvwrh(tt, gono_params(1:2), vea = 1, dur_v = 4, r1 = 1,
                         r2 = 1, strategy = "VoD", vbe = 1)
  z <- extract_flows_xpvwrh(y)

  expect_equal(z$cum_treated[1, ], z$treated[1, ])
  expect_equal(z$cum_treated[2, ] - z$cum_treated[1, ], z$treated[2, ])
  expect_equal(z$vaccinated, t(aggregate(y, "cum_vaccinated", as_incid = TRUE)))
  expect_equal(t(aggregate(y, "cum_vaccinated", as_incid = TRUE)),
               (z$primary_total + z$part_to_full +
                  z$revaccinated))

  # number revaccinated at t = 2, is the difference in the # cumulatively
  # vaccinated in strata W (4) between t = 2 and t = 3
  expect_equal(z$revaccinated[2, ],
               sapply(y, function(x) diff(rowSums(x$cum_vaccinated[-1, , 4]))))
  # asymptomatic and symptomatic diagnoses are calculated for vaccine protected
  # and non-vaccine protected strata correctly when adjusted baseline supplied
  # (extract_flows_xpvwrh() is within compare_baseline_xpvwrh() and baseline
  # comparison occurs within the latter after the former called )

  tt <- seq(0, 2)
  y  <- run_onevax_xpvwrh(tt, gono_params(1:2), vea = 1, dur_v = 4, r1 = 1,
                          r2 = 1, strategy = "VoD")
  y0 <- run_onevax_xpvwrh(tt, gono_params(1:2), vea = 0, dur_v = 4, r1 = 0,
                          r2 = 0, strategy = "VoD")

  baseline <- extract_flows_xpvwrh(adjust_baseline(y0, y))
  cp <- list(qaly_loss_per_diag_s = c(0.002, 0.001),
             unit_cost_manage_symptomatic = c(98, 99),
             unit_cost_manage_asymptomatic = c(93, 94),
             unit_cost_screen_uninfected = c(70, 71))

  z <- compare_baseline_xpvwrh(y, baseline, 1, 1, cp, 0, 1, 1)

  # inc_diag_a = inc_diag_a_xwh + inc_diag_a_pvr (same for s)
  expect_equal(z$inc_diag_a_xwh + z$inc_diag_a_pvr, z$inc_diag_a)
  expect_equal(z$inc_diag_s_xwh + z$inc_diag_s_pvr, z$inc_diag_s)

  # adjusting baseline hasn't broken asymp and symp diagoses totals, expect
  # sum of diagnoses to be similar to inc_treated
  # (won't be exact as diagnoses occurs before treatment in a timepoint)
  expect_equal(z$inc_diag_a + z$inc_diag_s, z$inc_treated, tolerance = 0.05)

  # for vaccine uptake > 1 and efficacy > 1, inc_diag_a_pvr and
  # inc_diag_s_pvr -'ve
  expect_true(all(z$inc_diag_a_pvr < 0))
  expect_true(all(z$inc_diag_s_pvr < 0))

  # sense check: adjusted baseline should give same resutls as vaccine with
  # 0 uptake
  # if vaccination uptake = 0 in vaccine run, inc_diag_a_pvr, inc_diag_s_pvr = 0
  # because no one is entering vaccine strata in the model run
  # and inc_diag_a_xwh, inc_diag_s_xwh = 0
  # because there is no difference in diagnoses between the adjusted baseline
  # and the model run with no uptake

  baseline <- extract_flows_xpvwrh(adjust_baseline(y0, y0))
  z <- compare_baseline_xpvwrh(y0, baseline, 1, 1, cp, 0, 1, 1)

  expect_true(all(round(z$inc_diag_a_xwh, 10) == 0))
  expect_true(all(round(z$inc_diag_s_xwh, 10) == 0))
  expect_true(all(round(z$inc_diag_a_pvr, 10) == 0))
  expect_true(all(round(z$inc_diag_s_pvr, 10) == 0))

  #extract_flows_xpvwrh() works with n_erlang > 1. Doesn't throw error
  #2. If people overall take longer to wane (with n_erlang > 1)
  # fewer diagnoses overall as more people protected at a given time

  y_erlang <- run_onevax_xpvwrh(tt, gono_params(1:2), vea = 1, dur_v = 4,
                                r1 = 1, r2 = 1, strategy = "VoD", n_erlang = 3)
  y0_erlang <- run_onevax_xpvwrh(tt, gono_params(1:2), vea = 0, dur_v = 4,
                                 r1 = 0, r2 = 0, strategy = "VoD", n_erlang = 3)
  baseline_erlang <- extract_flows_xpvwrh(adjust_baseline(y0_erlang, y_erlang))
  z_erlang <- compare_baseline_xpvwrh(y_erlang, baseline_erlang,
                                      1, 1, cp, 0, 1, 1)

  y <- run_onevax_xpvwrh(tt, gono_params(1:2), vea = 1, dur_v = 4,
                         r1 = 1, r2 = 1, strategy = "VoD", n_erlang = 1)
  y0 <- run_onevax_xpvwrh(tt, gono_params(1:2), vea = 0, dur_v = 4,
                          r1 = 0, r2 = 0, strategy = "VoD", n_erlang = 1)

  baseline <- extract_flows_xpvwrh(adjust_baseline(y0, y))
  z <- compare_baseline_xpvwrh(y, baseline, 1, 1, cp, 0, 1, 1)

  #overall
  expect_true(all(z_erlang$inc_diag_a < z$inc_diag_a))
  expect_true(all(z_erlang$inc_diag_s < z$inc_diag_s))
})

test_that("extract_flows works", {
  tt <- seq(0, 2)
  y <- run_onevax_xvwv(tt, gono_params(1:2), vea = 1, dur = 4, uptake = 1,
                       strategy = "VoD", vbe = 1)
  z <- extract_flows(y)

  expect_equal(z$cum_treated[1, ], z$treated[1, ])
  expect_equal(z$cum_treated[2, ] - z$cum_treated[1, ], z$treated[2, ])
  expect_equal(z$vaccinated, t(aggregate(y, "cum_vaccinated", as_incid = TRUE)))
  expect_equal(z$revaccinated[2, ],
               sapply(y, function(x) diff(rowSums(x$cum_vaccinated[-1, , 3]))))
  expect_equal(z$offered_primary, z$vaccinated - z$revaccinated - z$vbe)
})

test_that("extract_flows_trial works", {

  #n_erlang = 1, n_diag_rec = 1
  gp <- gono_params_trial(1:3)
  n_erlang <- 1
  n_diag_rec <- 1
  tt <- seq.int(0, 10)
  N <- 600
  set.seed(1)

  #stochastic
  y1s <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                              vea = 0.5, vei = 0, ved = 0, ves = 0,
                              n_erlang = n_erlang, n_diag_rec = n_diag_rec,
                              stochastic = TRUE, N = N)

  y1s_flows <- extract_flows_trial(y1s)

  #core code from extract_flows is unchanged , and works for stochastic
  z <- y1s_flows
  idx <- stratum_index_xvw_trial(n_erlang, n_diag_rec)
  expect_equal(z$cum_treated[1, ], z$treated[1, ])
  expect_equal(z$cum_treated[2, ] - z$cum_treated[1, ], z$treated[2, ])
  #same for new outputs
  expect_equal(z$diag_a_X_all_diaghist, t(aggregate(y1s, "cum_diag_a",
                                                    stratum = idx$X,
                                                    stochastic = TRUE,
                                                    as_incid = TRUE)))
  expect_equal(z$diag_a_VW_all_diaghist, t(aggregate(y1s, "cum_diag_a",
                                                     stratum = c(idx$V, idx$W),
                                                     stochastic = TRUE,
                                                     as_incid = TRUE)))
  expect_equal(z$diag_s_X_all_diaghist, t(aggregate(y1s, "cum_diag_s",
                                                    stratum = idx$X,
                                                    stochastic = TRUE,
                                                    as_incid = TRUE)))
  expect_equal(z$diag_s_VW_all_diaghist, t(aggregate(y1s, "cum_diag_s",
                                                     stratum = c(idx$V, idx$W),
                                                     stochastic = TRUE,
                                                     as_incid = TRUE)))

  #deterministic
  y1d <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                              vea = 0.5, vei = 0, ved = 0, ves = 0,
                              n_erlang = n_erlang, n_diag_rec = n_diag_rec,
                              stochastic = FALSE, N = N)

  y1d_flows <- extract_flows_trial(y1d)

  #core code from extract_flows is unchanged, and works for deterministic
  z <- y1d_flows
  idx <- stratum_index_xvw_trial(n_erlang, n_diag_rec)
  expect_equal(z$cum_treated[1, ], z$treated[1, ])
  expect_equal(z$cum_treated[2, ] - z$cum_treated[1, ], z$treated[2, ])
  #same for new outputs
  expect_equal(z$diag_a_X_all_diaghist, t(aggregate(y1d, "cum_diag_a",
                                                    stratum = idx$X,
                                                    as_incid = TRUE)))
  expect_equal(z$diag_a_VW_all_diaghist, t(aggregate(y1d, "cum_diag_a",
                                                     stratum = c(idx$V, idx$W),
                                                     as_incid = TRUE)))
  expect_equal(z$diag_s_X_all_diaghist, t(aggregate(y1d, "cum_diag_s",
                                                    stratum = idx$X,
                                                    as_incid = TRUE)))
  expect_equal(z$diag_s_VW_all_diaghist, t(aggregate(y1d, "cum_diag_s",
                                                     stratum = c(idx$V, idx$W),
                                                     as_incid = TRUE)))


  #n_erlang = 1, n_diag_rec > 1
  n_diag_rec <- 2
  idx <- stratum_index_xvw_trial(n_erlang, n_diag_rec)
  set.seed(1)

  #stochastic
  y2s <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                              vea = 0.5, vei = 0, ved = 0, ves = 0,
                              n_erlang = n_erlang, n_diag_rec = n_diag_rec,
                              stochastic = TRUE, N = N)

  y2s_flows <- extract_flows_trial(y2s)
  #deterministic
  y2d <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                              vea = 0.5, vei = 0, ved = 0, ves = 0,
                              n_erlang = n_erlang, n_diag_rec = n_diag_rec,
                              stochastic = FALSE, N = N)

  y2d_flows <- extract_flows_trial(y2d)

  #for n_diag_rec > 1, cumulative diagnoses in first diagnosis history strata
  #is < diagnoses across all diagnosis history strata
  #and the same for cumulative incidence

  #stochastic
  expect_true(all(y2s_flows$cum_diag_a_X_first_diag_hist
                  < y2s_flows$cum_diag_a_X_all_diaghist))
  expect_true(all(y2s_flows$cum_diag_a_VW_first_diag_hist
                  < y2s_flows$cum_diag_a_VW_all_diaghist))
  expect_true(all(y2s_flows$cum_diag_s_X_first_diag_hist
                  < y2s_flows$cum_diag_s_X_all_diaghist))
  expect_true(all(y2s_flows$cum_diag_s_VW_first_diag_hist
                  < y2s_flows$cum_diag_s_VW_all_diaghist))
  expect_true(all(y2s_flows$cum_incid_X_first_diag_hist
                  < y2s_flows$cum_incid_X_all_diaghist))
  expect_true(all(y2s_flows$cum_incid_VW_first_diag_hist
                  < y2s_flows$cum_incid_VW_all_diaghist))

  #deterministic - no migration to diagnosis history strata so first diagnosis
  #history strata should = number diagnoses across all diagnosis history
  expect_true(all(y2d_flows$cum_diag_a_X_first_diag_hist
                  < y2d_flows$cum_diag_a_X_all_diaghist))
  expect_true(all(y2d_flows$cum_diag_a_VW_first_diag_hist
                  < y2d_flows$cum_diag_a_VW_all_diaghist))
  expect_true(all(y2d_flows$cum_diag_s_X_first_diag_hist
                  < y2d_flows$cum_diag_s_X_all_diaghist))
  expect_true(all(y2d_flows$cum_diag_s_VW_first_diag_hist
                  < y2d_flows$cum_diag_s_VW_all_diaghist))
  expect_true(all(y2d_flows$cum_incid_X_first_diag_hist
                  < y2d_flows$cum_incid_X_all_diaghist))
  expect_true(all(y2d_flows$cum_incid_VW_first_diag_hist
                  < y2d_flows$cum_incid_VW_all_diaghist))

  #n_erlang > 1, n_diag_rec > 1
  n_diag_rec <- 2
  n_erlang <- 3
  idx <- stratum_index_xvw_trial(n_erlang, n_diag_rec)
  set.seed(1)

  #stochastic
  y3s <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                              vea = 0.5, vei = 0, ved = 0, ves = 0,
                              n_erlang = n_erlang, n_diag_rec = n_diag_rec,
                              stochastic = TRUE, N = N)
  y3s_flows <- extract_flows_trial(y3s)

  #deterministic
  y3d <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                              vea = 0.5, vei = 0, ved = 0, ves = 0,
                              n_erlang = n_erlang, n_diag_rec = n_diag_rec,
                              stochastic = FALSE, N = N)

  y3d_flows <- extract_flows_trial(y3d)

  #when n_erlang > 1  and n_diag_rec > 1 tests still pass
  #stochastic
  z <- y3s_flows
  expect_equal(z$diag_a_X_all_diaghist, t(aggregate(y3s, "cum_diag_a",
                                                    stratum = idx$X,
                                                    stochastic = TRUE,
                                                    as_incid = TRUE)))
  expect_equal(z$diag_a_VW_all_diaghist, t(aggregate(y3s, "cum_diag_a",
                                                     stratum = c(idx$V, idx$W),
                                                     stochastic = TRUE,
                                                     as_incid = TRUE)))
  expect_equal(z$diag_s_X_all_diaghist, t(aggregate(y3s, "cum_diag_s",
                                                    stratum = idx$X,
                                                    stochastic = TRUE,
                                                    as_incid = TRUE)))

  expect_equal(z$diag_s_VW_all_diaghist, t(aggregate(y3s, "cum_diag_s",
                                                     stratum = c(idx$V, idx$W),
                                                     stochastic = TRUE,
                                                     as_incid = TRUE)))

  expect_true(all(y3s_flows$cum_diag_a_X_first_diag_hist
                  < y3s_flows$cum_diag_a_X_all_diaghist))
  expect_true(all(y3s_flows$cum_diag_a_VW_first_diag_hist
                  < y3s_flows$cum_diag_a_VW_all_diaghist))
  expect_true(all(y3s_flows$cum_diag_s_X_first_diag_hist
                  < y3s_flows$cum_diag_s_X_all_diaghist))
  expect_true(all(y3s_flows$cum_diag_s_VW_first_diag_hist
                  < y3s_flows$cum_diag_s_VW_all_diaghist))

  #deterministic
  z <- y3d_flows
  expect_equal(z$diag_a_X_all_diaghist, t(aggregate(y3d, "cum_diag_a",
                                                    stratum = idx$X,
                                                    as_incid = TRUE)))
  expect_equal(z$diag_a_VW_all_diaghist, t(aggregate(y3d, "cum_diag_a",
                                                     stratum = c(idx$V, idx$W),
                                                     as_incid = TRUE)))
  expect_equal(z$diag_s_X_all_diaghist, t(aggregate(y3d, "cum_diag_s",
                                                    stratum = idx$X,
                                                    as_incid = TRUE)))
  expect_equal(z$diag_s_VW_all_diaghist, t(aggregate(y3d, "cum_diag_s",
                                                     stratum = c(idx$V, idx$W),
                                                     as_incid = TRUE)))

  expect_true(all(y3d_flows$cum_diag_a_X_first_diag_hist
                  < y3d_flows$cum_diag_a_X_all_diaghist))
  expect_true(all(y3d_flows$cum_diag_a_VW_first_diag_hist
                  < y3d_flows$cum_diag_a_VW_all_diaghist))
  expect_true(all(y3d_flows$cum_diag_s_X_first_diag_hist
                  < y3d_flows$cum_diag_s_X_all_diaghist))
  expect_true(all(y3d_flows$cum_diag_s_VW_first_diag_hist
                  < y3d_flows$cum_diag_s_VW_all_diaghist))

  #number of person years spent exposed is increasing over time
  #person years at final timepiont > person years at the start

  expect_true(all(y3d_flows$N_person_yrs_exp_X.I[1, ]
                  < y3d_flows$N_person_yrs_exp_X.I[10, ]))
  expect_true(all(y3d_flows$N_person_yrs_exp_VW.I[1, ]
                  < y3d_flows$N_person_yrs_exp_VW.I[10, ]))
  expect_true(all(y3s_flows$N_person_yrs_exp_X.I[1, ]
                  < y3s_flows$N_person_yrs_exp_X.I[10, ]))
  expect_true(all(y3s_flows$N_person_yrs_exp_VW.I[1, ]
                  < y3s_flows$N_person_yrs_exp_VW.I[10, ]))

  #if vea = 1, person years exposed in V, increases by the number of people in
  #the trial arm each year (N/2)

  N <- 600
  set.seed(1)
  y4s <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                              vea = 1, vei = 0, ved = 0, ves = 0,
                              n_erlang = n_erlang, n_diag_rec = n_diag_rec,
                              stochastic = TRUE, N = N)
  y4s_flows <- extract_flows_trial(y4s)

  expect_true(all(round(diff(y4s_flows$N_person_yrs_exp_VW.I), 0) == N / 2))

})

test_that("generating indices for never-diagnosed VW strata is correct", {

  #for n_erlang 1 but diagnosis history > 1
  n_erlang <- 1
  n_diag_rec <- 3
  gp <- gono_params_trial(1:3)
  tt <- seq.int(0, 10)
  N <- 600
  set.seed(1)

  #stochastic
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0.5, vei = 0, ved = 0, ves = 0,
                            n_erlang = n_erlang, n_diag_rec = n_diag_rec,
                            stochastic = TRUE, N = N)

  #code for strat, n_diag_hist and n_erlang as is in extract_flows_trial()
  strata <- dimnames(y[[1]]$N)[[3]]
  n_diag_rec <- sum(grepl("W", strata)) # count diag hist categories
  n_erlang <- sum((grepl("V", strata))) / n_diag_rec # count erlang
  idx <- stratum_index_xvw_trial(n_erlang, n_diag_rec)
  idx$never_diag <- seq(idx$V[1], by = n_diag_rec,
                        length.out = n_erlang + 1)

  expect_true(all(idx$never_diag == c(4, 7)))

  # same but n_erlang > 1
  n_erlang <- 4
  n_diag_rec <- 3
  gp <- gono_params_trial(1:3)
  tt <- seq.int(0, 10)
  N <- 600
  set.seed(1)

  #stochastic
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0.5, vei = 0, ved = 0, ves = 0,
                            n_erlang = n_erlang, n_diag_rec = n_diag_rec,
                            stochastic = TRUE, N = N)

  strata <- dimnames(y[[1]]$N)[[3]]
  n_diag_rec <- sum(grepl("W", strata)) # count diag hist categories
  n_erlang <- sum((grepl("V", strata))) / n_diag_rec # count erlang
  idx <- stratum_index_xvw_trial(n_erlang, n_diag_rec)
  idx$never_diag <- seq(idx$V[1], by = n_diag_rec,
                        length.out = n_erlang + 1)

  expect_true(all(idx$never_diag == c(4, 7, 10, 13, 16)))

  # same but diagnosis history is 1
  n_erlang <- 4
  n_diag_rec <- 1
  gp <- gono_params_trial(1:3)
  tt <- seq.int(0, 10)
  N <- 600
  set.seed(1)

  #stochastic
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0.5, vei = 0, ved = 0, ves = 0,
                            n_erlang = n_erlang, n_diag_rec = n_diag_rec,
                            stochastic = TRUE, N = N)

  strata <- dimnames(y[[1]]$N)[[3]]
  n_diag_rec <- sum(grepl("W", strata)) # count diag hist categories
  n_erlang <- sum((grepl("V", strata))) / n_diag_rec # count erlang
  idx <- stratum_index_xvw_trial(n_erlang, n_diag_rec)
  idx$never_diag <- seq(idx$V[1], by = n_diag_rec,
                        length.out = n_erlang + 1)

  expect_true(all(idx$never_diag == c(2, 3, 4, 5, 6)))

})

test_that("gonovax_year works as expected", {
  expect_equal(gonovax_year(2009), 0)
  expect_error(gonovax_year(2006),
               "Negative dates, gonovax_year likely applied twice")
  expect_equal(gonovax_year_as_year(2), 2011)
  expect_error(gonovax_year_as_year(-1), "'gonovax_year' must be at least 1")
})


test_that("dbetabinom", {

  ## when prob = 1/2, rho = 1/3 (equivalently a = b = 1),
  ## equivalent to discrete uniform from 0 to size
  expect_equal(dbetabinom(4, 35, 1 / 2, 1 / 3), 1 / 36)
  expect_equal(dbetabinom(4, 35, 1 / 2, 1 / 3, log = TRUE), log(1 / 36))


  ## compare with extraDistr::dbbinom (uses a and b as parameters - to match
  ## use prob = a / (a + b) and rho = 1 / (a + b + 1))
  f <- function(x, size, a, b, log = FALSE) {
    prob  <- a / (a + b)
    rho <- 1 / (a + b + 1)

    dbetabinom(x, size, prob, rho, log)
  }

  ## compare to extraDistr::dbbinom(15, 63, 2, 5) ~ 0.03613356
  expect_equal(f(15, 63, 2, 5), 0.03613356, tolerance = 1e-8)
  ## compare to extraDistr::dbbinom(15, 63, 2, 5, log = TRUE) ~ -3.320533
  expect_equal(f(15, 63, 2, 5, log = TRUE), -3.320533, tolerance = 1e-6)

  ## compare to extraDistr::dbbinom(672, 50454, 3, 2) ~ 4.18089e-08
  expect_equal(f(672, 50454, 3, 2), 4.18089e-08, tolerance = 1e-13)
  ## compare to extraDistr::dbbinom(15, 63, 2, 5, log = TRUE) ~ -16.99016
  expect_equal(f(672, 50454, 3, 2, log = TRUE), -16.99016, tolerance = 1e-5)
})

test_that("adjust_baseline is working as expected", {

  tt <- seq(0, 2)
  gp <- gono_params(1)

  # partially effective vaccine that everyone accepts when offered, all move to
  # V
  y <- run_onevax_xpvwrh(tt, gp, vea = 0.5, dur_v = 4, r1 = 1,
                         r2 = 0.5, strategy = "VoD")

  # no vaccine uptake = nobody in P, V, W, or R
  y0 <- run_onevax_xpvwrh(tt, gp, vea = 0, dur_v = 4, r1 = 0,
                          r2 = 0, strategy = "VoD")

  y0_adj <- adjust_baseline(y0, y)

  #expect y0 to have absence of people in P, V, W and R but y0_adjust to have
  #people in P, V, W, and R (for asymp and symp diagnoses)

  idx <- stratum_index_xpvwrh(n_erlang = 1)
  #a
  expect_true(all(y0[[1]]$cum_diag_a[, , c(idx$P, idx$V, idx$W, idx$R)] == 0))
  expect_true(all(y0_adj[[1]]$cum_diag_a[-1, , c(idx$P,
                                                 idx$V, idx$W, idx$R)] != 0))
  #s
  expect_true(all(y0[[1]]$cum_diag_s[, , c(idx$P, idx$V, idx$W, idx$R)] == 0))
  expect_true(all(y0_adj[[1]]$cum_diag_s[-1, , c(idx$P,
                                                 idx$V, idx$W, idx$R)] != 0))

  #expect the overall number of diagnoses in y0 and y0_adj to be the same
  #for each time point
  #(they should simply be spread across strata)
  expect_equal(rowSums(y0[[1]]$cum_diag_a[, , c(idx$X)]),
               rowSums(y0_adj[[1]]$cum_diag_a[, , c(idx$X, idx$P,
                                                    idx$V, idx$W, idx$R)]))
  expect_equal(rowSums(y0[[1]]$cum_diag_s[, , c(idx$X)]),
               rowSums(y0_adj[[1]]$cum_diag_s[, , c(idx$X, idx$P,
                                                    idx$V, idx$W, idx$R)]))

})
