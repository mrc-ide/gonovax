######### testing #######
context("trial model check")

test_that("no new infections if lambda is 0", {

  # no vaccination
  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])
  params$lambda <- 0
  mod <- model_trial$new(user = params, unused_user_action = "ignore")
  tt <- seq.int(0, 5) / 365
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)


  expect_true(all(y$I == 0))             #incubating remain at 0
  expect_true(all(y$S == 0))             #no symptomatic infections
  expect_true(all(y$cum_incid == 0))
  expect_true(all(unlist(y) >= 0))                     #no compartments negative
  expect_true(all(round(rowSums(y$N)) == params$N0))   #number of participants
  #doesn't change

  # presence of vaccination
  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]],
                               vax_params = vax_params_xvw_trial(), p_v = 0.5)
  tt <- seq.int(0, 5) / 365
  params$lambda <- 0
  mod2 <- model_trial$new(user = params, unused_user_action = "ignore")
  y <- mod2$run(t = tt)

  y <- mod2$transform_variables(y)

  expect_true(all(y$I == 0))
  expect_true(all(y$S == 0))
  expect_true(all(y$cum_incid == 0))

  expect_true(all(unlist(y) >= 0))
  expect_true(all(round(rowSums(y$N)) == params$N0))

})

test_that("number of individuals always >= 0 even at high lambda", {
  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])

  params$lambda <- 10
  mod <- model_trial$new(user = params, unused_user_action = "ignore")
  tt <- seq.int(0, 5) / 365
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)

  expect_true(all(unlist(y) >= 0))
  expect_true(all(round(rowSums(y$N)) == params$N0))

})

test_that("all individuals are in the high activity group", {

  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])
  mod <- model_trial$new(user = params, unused_user_action = "ignore")
  tt <- seq.int(0, 5) / 365
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)

  #low activity group empty
  expect_true(all(y$U[, 1, ] == 0))
  expect_true(all(y$I[, 1, ] == 0))
  expect_true(all(y$A[, 1, ] == 0))
  expect_true(all(y$S[, 1, ] == 0))
  expect_true(all(y$T[, 1, ] == 0))

  expect_true(all(unlist(y) >= 0))
  expect_true(all(round(rowSums(y$N)) == params$N0))

})


test_that("there are no symptomatic infections when psi = 0", {
  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])
  params$psi <- 0
  mod <- model_trial$new(user = params, unused_user_action = "ignore")
  tt <- seq.int(0, 5) / 365
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  expect_true(any(y$I > 0))
  expect_true(all(y$S == 0))
  expect_true(all(y$cum_diag_s == 0))
  expect_true(all(y$cum_diag_a[-1, , ] >= 0))
  expect_true(all(unlist(y) >= 0))

  expect_true(all(unlist(y) >= 0))
  expect_true(all(round(rowSums(y$N)) == params$N0))

})

test_that("there are no asymptomatic infections when psi = 1", {
  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])
  params$psi <- 1                   #proportion symptomatic
  params$S0[, ] <- params$A0[, ]    #moves starting individuals from A to S
  params$A0[, ] <- 0                #so still infections seeded
  mod <- model_trial$new(user = params, unused_user_action = "ignore")
  tt <- seq.int(0, 5) / 365
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  expect_true(any(y$I > 0))
  expect_true(all(y$A == 0))
  expect_true(all(y$cum_diag_a == 0))
  expect_true(all(y$cum_diag_s[-1, , ] >= 0))

  expect_true(all(unlist(y) >= 0))
  expect_true(all(round(rowSums(y$N)) == params$N0))

})

test_that("all individuals are uninfected at t = 0", {
  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])
  mod <- model_trial$new(user = params, unused_user_action = "ignore")
  tt <- seq.int(0, 5) / 365
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)

  expect_true(all(y$I[1, , ] == 0))
  expect_true(all(y$A[1, , ] == 0))
  expect_true(all(y$S[1, , ] == 0))
  expect_true(all(y$T[1, , ] == 0))
  expect_true(all(y$U[1, 2, ] > 0))

  expect_true(all(unlist(y) >= 0))
  expect_true(all(round(rowSums(y$N)) == params$N0))

})

test_that("no-one is treated when mu and eta = 0", {
  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])
  params$mu <- params$eta <-  0
  mod <- model_trial$new(user = params, unused_user_action = "ignore")

  tt <- seq.int(0, 5) / 365
  y <- mod$run(tt)
  y <- mod$transform_variables(y)
  expect_true(all(y$T == 0))
  expect_true(all(y$cum_treated == 0))

  expect_true(all(unlist(y) >= 0))
  expect_true(all(round(rowSums(y$N)) == params$N0))

})


test_that("Model works with vaccination and waning", {
  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, 5) / 365
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0, vei = 0, ved = 0, ves = 0)

  #Total N in high group in X and  V is equal and non zero at t = 0
  expect_true(all(y[[1]]$N[1, 2, 1] == y[[1]]$N[1, 2, 2]))
  expect_true(all(y[[1]]$N[1, 2, 1] > 0))
  expect_true(all(y[[1]]$N[1, 2, 2] > 0))

  #Total N in waning is 0 at t =0 and increases over time
  expect_true(all(y[[1]]$N[1, 2, 3] == 0))
  length_tt <- length(y[[1]]$N[, 2, 3])
  expect_true(all(y[[1]]$N[2:length_tt, 2, 3] > 0))

  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])
  expect_true(all(unlist(y) >= 0))
  expect_true(all(round(rowSums(y[[1]]$N)) == params$N0))

})


test_that("VEa behaves as expected ", {
  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, 5) / 365
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 1, vei = 0, ved = 0, ves = 0)

  # VEa = 1, no infections in V, X >0, W >0 (after t=0)

  expect_true(all(y[[1]]$cum_incid[, 2, 2] == 0))
  expect_true(all(y[[1]]$cum_incid[2:6, 2, 1] > 0))
  expect_true(all(y[[1]]$cum_incid[2:6, 2, 3] > 0))

  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])

  expect_true(all(unlist(y) >= 0))
  expect_true(all(round(rowSums(y[[1]]$N)) == params$N0))

  # VEa = 0, infections in X = V + W
  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, 5) / 365
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0, vei = 0, ved = 0, ves = 0)


  x <- as.numeric(y[[1]]$cum_incid[6, 2, 1])
  vw <- as.numeric(y[[1]]$cum_incid[6, 2, 2]) +
    as.numeric(y[[1]]$cum_incid[6, 2, 3])
  expect_equal(x, vw)

  expect_true(all(unlist(y) >= 0))
  expect_true(all(round(rowSums(y[[1]]$N)) == params$N0))

})


test_that("VEs behaves as expected ", {
  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, 5) / 365
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0, vei = 0, ved = 0, ves = 1)


  # VEs = 1 symptomatic diagnoses in V = 0, but for X & W > 0
  expect_true(all(y[[1]]$cum_diag_s[, 2, 2] == 0))
  expect_true(all(y[[1]]$cum_diag_s[2:6, 2, 1] > 0))
  expect_true(all(y[[1]]$cum_diag_s[2:6, 2, 3] > 0))

  # VEs = 1 asymptomatic diagnoses in V > 0
  expect_true(all(y[[1]]$cum_diag_a[2:6, 2, 2] > 0))

  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])

  expect_true(all(unlist(y) >= 0))
  expect_true(all(round(rowSums(y[[1]]$N)) == params$N0))

})


test_that("VEd behaves as expected ", {
  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, 5) / 365
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0, vei = 0, ved = 1, ves = 0)

  # VEd = 1 cumulative incid in V > X + W

  v <- y[[1]]$cum_incid[6, 2, 2]
  xw <- y[[1]]$cum_incid[6, 2, 1] + y[[1]]$cum_incid[6, 2, 3]

  expect_true(v > xw)

  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])

  expect_true(all(unlist(y) >= 0))
  expect_true(all(round(rowSums(y[[1]]$N)) == params$N0))

})


test_that("VEi behaves as expected ", {
  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, 5) / 365
  y0 <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                             vea = 0, vei = 0, ved = 0, ves = 0)

  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0, vei = 1, ved = 0, ves = 0)

  # Incidence in V is the same for VEi = 0 and VEi = 1

  expect_equal(y0[[1]]$cum_incid[6, 2, 2],
               y[[1]]$cum_incid[6, 2, 2])

  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])

  expect_true(all(unlist(y) >= 0))
  expect_true(all(round(rowSums(y[[1]]$N)) == params$N0))

  expect_true(all(unlist(y0) >= 0))
  expect_true(all(round(rowSums(y0[[1]]$N)) == params$N0))

})


test_that("aggregated time series output correctly", {
  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])
  mod <- model_trial$new(user = params, unused_user_action = "ignore")
  tt <- seq.int(0, 5) / 365
  y <- mod$run(tt)
  y <- mod$transform_variables(y)
  expect_equal(y$tot_treated, apply(y$cum_treated, 1, sum))
  expect_equal(y$tot_attended, apply(y$cum_screened, 1, sum) + y$tot_treated)

  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])

  expect_true(all(unlist(y) >= 0))
  expect_true(all(round(rowSums(y$N)) == params$N0))

})

test_that("outputs for vaccination run as expected", {

  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, 5) / 365
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0.5, vei = 0.5, ved = 0.5, ves = 0.5)

  expect_equal(y[[1]]$U, array(c(0, 0, 0, 0, 0, 0, 3e+05, 298770.012957928,
                                 297547.433626743, 296334.908051234,
                                 295135.061477033, 293950.297077376, 0, 0, 0,
                                 0, 0, 0, 3e+05, 299383.713797674,
                                 298770.717502468, 298162.935800213,
                                 297562.116361515, 296969.7952273, 0, 0, 0, 0,
                                 0, 0, 0, 0.819388141002698, 1.63374216215416,
                                 2.4431142416253, 3.24757521226149,
                                 4.04721268804341),
                               dim = c(6, 2, 3), dimnames = list(NULL,
                                                                 c("L", "H"),
                                                                 c("X.I",
                                                                   "V1.I",
                                                                   "W.I"))))
})


test_that("n_erlang = n is working as expected", {

  #when n_erlang = n, nvax = n + 2, as generated in stratum_index_xvw_trial
  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, 5) / 365
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0.5, vei = 0.5, ved = 0.5, ves = 0.5)

  expect_equal(dim(y[[1]]$N)[3], 3)

  n_erlang <- 2
  y2 <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                             vea = 0.5, vei = 0.5, ved = 0.5, ves = 0.5,
                             n_erlang = n_erlang)

  expect_equal(dim(y2[[1]]$N)[3], n_erlang + 2)
  idx <- stratum_index_xvw_trial(n_erlang)
  expect_equal(dim(y2[[1]]$N)[3], idx$n_vax)

  n_erlang <- 5
  y3 <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                             vea = 0.5, vei = 0.5, ved = 0.5, ves = 0.5,
                             n_erlang = n_erlang)
  expect_equal(dim(y3[[1]]$N)[3], n_erlang + 2)
  idx <- stratum_index_xvw_trial(n_erlang)
  expect_equal(dim(y3[[1]]$N)[3], idx$n_vax)

  #for n_diag_rec > 2
  n_erlang <- 2
  n_diag_rec <- 2
  tt <- seq.int(0, 1)
  y_2_2 <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                                vea = 0.5, vei = 0.5, ved = 0.5, ves = 0.5,
                                n_erlang = n_erlang,
                                stochastic = FALSE,
                                n_diag_rec = n_diag_rec)

  idx <- stratum_index_xvw_trial(n_erlang, n_diag_rec)
  expect_equal(idx$n_vax, 8)
  expect_equal(dim(y_2_2[[1]]$N)[3], idx$n_vax)

  #vax_params_xvw_trial generates the correct waning maps (in for n_erlang > 1)
  dur <- 1e03
  n_erlang <- 1 #n_vax is 3
  idx <- stratum_index_xvw_trial(n_erlang)
  erlang_1 <- vax_params_xvw_trial(vea = 1, dur = dur, n_erlang = n_erlang)$w
  matrix_1 <- matrix(data = (c(rep(0, 4), -n_erlang / dur, rep(0, 2),
                               n_erlang / dur, 0)), ncol = idx$n_vax,
                     byrow = TRUE)
  expect_equal(erlang_1, matrix_1)

  n_erlang <- 2 #n_vax is 4
  idx <- stratum_index_xvw_trial(n_erlang)
  erlang_2 <- vax_params_xvw_trial(vea = 1, dur = dur, n_erlang = n_erlang)$w
  matrix_2 <- matrix(data = (c(rep(0, 5), -n_erlang / dur, rep(0, 3),
                               n_erlang / dur, -n_erlang / dur, rep(0, 3),
                               n_erlang / dur, 0)), ncol = idx$n_vax,
                     byrow = TRUE)
  expect_equal(erlang_2, matrix_2)

  # when vea is perfect, there are no infections in any of the
  # erlang strata

  n_erlang <- 3   #n_vax is 5
  idx <- stratum_index_xvw_trial(n_erlang)
  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, 5) / 365
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 1, vei = 0, ved = 0, ves = 0,
                            n_erlang = n_erlang)

  expect_true(all(y[[1]]$cum_incid[, , 2] == 0)) #V1
  expect_true(all(y[[1]]$cum_incid[, , 3] == 0)) #V2
  expect_true(all(y[[1]]$cum_incid[, , 4] == 0)) #V3
  expect_true(all(y[[1]]$cum_incid[, , c(idx$V)] == 0)) #all V

  expect_true(all(y[[1]]$cum_incid[-1, 2, 1] != 0)) #X still has infections
  expect_true(all(y[[1]]$cum_incid[-1, 2, idx$X] != 0))
  expect_true(all(y[[1]]$cum_incid[-1, 2, 5] != 0)) #W not under protection, has
  expect_true(all(y[[1]]$cum_incid[-1, 2, idx$W] != 0)) # infections

  n_diag_rec <- 2
  n_erlang <- 3   #n_vax is 10
  idx <- stratum_index_xvw_trial(n_erlang, n_diag_rec)
  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, 1)
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1,
                            vea = 1, vei = 0, ved = 0, ves = 0,
                            n_erlang = n_erlang,
                            stochastic = FALSE,
                            n_diag_rec = n_diag_rec)

  expect_true(all(y[[1]]$cum_incid[, , c(idx$V)] == 0)) #all of V no infections
  expect_true(all(y[[1]]$cum_incid[-1, 2, idx$X] != 0)) #X still has infections
  expect_true(sum(y[[1]]$cum_incid[-1, 2, idx$W]) != 0) #W not under protection,
  #has infections

})

test_that("expected cumulative outputs are cumulative", {

  years <- 40
  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, 40)

  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0.5, vei = 0.5, ved = 0.5, ves = 0.5,
                            stochastic = FALSE)

  expect_true(all(diff(y[[1]]$cum_diag_a[c(years - 10, years + 1), 2, ]) > 0))
  expect_true(all(diff(y[[1]]$cum_incid[c(years - 10, years + 1), 2, ]) > 0))
  expect_true(all(diff(y[[1]]$cum_diag_s[c(years - 10, years + 1), 2, ]) > 0))
  expect_true(all(diff(y[[1]]$cum_treated[c(years - 10, years + 1), 2, ]) > 0))
  expect_true(all(diff(y[[1]]$cum_screened[c(years - 10, years + 1), 2, ]) > 0))
  expect_true(all(diff(y[[1]]$cum_pye_trial_pov[c(years - 10, years + 1), 2, ]) > 0))
  expect_true(all(diff(y[[1]]$cum_pye_true[c(years - 10, years + 1), 2, ]) > 0))

  n_diag_rec <- 2
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0.5, vei = 0.5, ved = 0.5, ves = 0.5,
                            stochastic = FALSE, n_diag_rec = n_diag_rec)

  expect_true(all(diff(y[[1]]$cum_diag_a[c(years - 10, years + 1), 2,
                                         c(2, 4, 6)]) > 0))
  expect_true(all(diff(y[[1]]$cum_incid[c(years - 10, years + 1),
                                        2, c(2, 4, 6)]) > 0))
  expect_true(all(diff(y[[1]]$cum_diag_s[c(years - 10, years + 1),
                                         2, c(2, 4, 6)]) > 0))
  expect_true(all(diff(y[[1]]$cum_treated[c(years - 10, years + 1),
                                          2, c(2, 4, 6)]) > 0))
  expect_true(all(diff(y[[1]]$cum_screened[c(years - 10, years + 1),
                                           2, c(2, 4, 6)]) > 0))
  expect_true(all(diff(y[[1]]$cum_pye_trial_pov[c(years - 10, years + 1),
                                                2, c(2, 4, 6)]) > 0))
  expect_true(all(diff(y[[1]]$cum_pye_true[c(years - 10, years + 1),
                                           2, c(2, 4, 6)]) > 0))

})


test_that("model outputs for a basic run haven't changed between updates", {

  gp <- gono_params_trial(1)[1]
  n_erlang <- 1
  tt <- seq.int(0, 1)
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0, vei = 0, ved = 0, ves = 0,
                            n_erlang = n_erlang,
                            stochastic = FALSE)

  expect_equal(y[[1]]$U,
               array(c(rep(0, 2), 300000, 250162.7, rep(0, 2), 300000, 249912.6,
                       rep(0, 3), 250.0376), dim = c(2, 2, 3),
                     dimnames = list(NULL, c("L", "H"),
                                     c("X.I", "V1.I", "W.I"))),
               tolerance = 1e-2)

  #and with n_erlang and diagnosis history:
  gp <- gono_params_trial(1)[1]
  n_erlang <- 2
  n_diag_rec <- 2
  tt <- seq.int(0, 1)
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0, vei = 0, ved = 0, ves = 0,
                            n_erlang = n_erlang,
                            stochastic = FALSE, n_diag_rec = n_diag_rec)

  expect_equal(y[[1]]$U,
               array(c(rep(0, 2), 300000, 136724.6, rep(0, 3),
                       113438.1, rep(0, 2), 300000, 136451.4,
                       rep(0, 3), 113211.5, rep(0, 3), 272.9028,
                       rep(0, 3), 226.4229, rep(0, 3),
                       0.2730848, rep(0, 3), 0.2265739),
                     dim = c(2, 2, 8),
                     dimnames = list(NULL, c("L", "H"),
                                     c("X.I", "X.II", "V1.I", "V1.II", "V2.I",
                                       "V2.II", "W.I", "W.II"))),
               tolerance = 1e-2)

})

test_that("correct number of individuals are set up in each trial arm", {

  gp <- gono_params_trial(1)[1]
  n_erlang <- 1
  tt <- seq.int(0, 1)

  N <- 500

  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0, vei = 0, ved = 0, ves = 0,
                            n_erlang = n_erlang,
                            stochastic = FALSE,
                            N = N)

  expect_equal(y[[1]]$N[1, 2, 1], N / 2)
  expect_equal(y[[1]]$N[1, 2, 2], N / 2)
  expect_equal(y[[1]]$N[1, 2, 3], 0)

  # for n_diag_rec > 1

  n_diag_rec <- 3
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0, vei = 0, ved = 0, ves = 0,
                            n_erlang = n_erlang,
                            stochastic = FALSE,
                            N = N, n_diag_rec = n_diag_rec)

  # N/2 in X.I and V1.I for t = 0
  expect_true(all(y[[1]]$N[1, 2, c(1, 4)] == N / 2))

  # 0 elsewhere
  expect_true(all(y[[1]]$N[1, 2, c(2, 3, 5:9)] == 0))

})

test_that("for n_diag_rec > 1, total N summed over X or V+W is
          the same and correct", {

            gp <- gono_params_trial(1)[1]
            n_erlang <- 1
            tt <- seq.int(0, 5)
            n_diag_rec <- 3
            N <- 6e+05

            y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e03,
                                      vea = 0, vei = 0, ved = 0, ves = 0,
                                      n_erlang = n_erlang,
                                      stochastic = FALSE,
                                      n_diag_rec = n_diag_rec, N = N)

            #X.I, X.II, X.III
            expect_equal(rowSums(y[[1]]$N[, 2, 1:3]), rep(N / 2, length(tt)))

            #V1.I, V1.II, V1.III, W.I, W.II, WIII
            expect_equal(rowSums(y[[1]]$N[, 2, 4:9]), rep(N / 2, length(tt)))

            # and for n_erlang > 1
            n_erlang <- 2
            y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e03,
                                      vea = 0, vei = 0, ved = 0, ves = 0,
                                      n_erlang = n_erlang,
                                      stochastic = FALSE,
                                      n_diag_rec = n_diag_rec, N = N)

            #X.I, X.II, X.III
            expect_equal(rowSums(y[[1]]$N[, 2, 1:3]), rep(N / 2, length(tt)))
            expect_equal(rowSums(y[[1]]$N[, 2, 4:12]), rep(N / 2, length(tt)))

          })

test_that("for n_diag_rec > 1, the number treated = the number recorded
          as diagnosed", {

            gp <- gono_params_trial(1)[1]
            n_erlang <- 1
            tt <- seq.int(0, 5)
            n_diag_rec <- 2
            N <- 6e+05

            #no waning for simplicity
            y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e100000000,
                                      vea = 0, vei = 0, ved = 0, ves = 0,
                                      n_erlang = n_erlang,
                                      stochastic = FALSE,
                                      n_diag_rec = n_diag_rec, N = N)

            # number treated in '.I' becomes the number in '.II'
            expect_true(all(round(y[[1]]$cum_treated[, 2, 1], 3) ==
                              round(y[[1]]$N[, 2, 2], 3)))
            expect_true(all(round(y[[1]]$cum_treated[, 2, 3], 3) ==
                              round(y[[1]]$N[, 2, 4], 3)))
            expect_true(all(round(y[[1]]$cum_treated[, 2, 5], 3) ==
                              round(y[[1]]$N[, 2, 6], 3)))

            n_diag_rec <- 3

            y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e100000000,
                                      vea = 0, vei = 0, ved = 0, ves = 0,
                                      n_erlang = n_erlang,
                                      stochastic = FALSE,
                                      n_diag_rec = n_diag_rec, N = N)

            # number treated each year in '.I' becomes the N gained in '.II' AND
            # '.III' diagnosis history strata for that year
            # (because some people could get diagnosed twice in one year)
            # and number treated in '.II' becomes the number in '.III' diagnosis
            # history

            expect_true(all(round(diff(y[[1]]$cum_treated[, 2, 1]), 3) ==
                              round(diff(rowSums(y[[1]]$N[, 2, 2:3])), 3)))
            expect_true(all(round(diff(y[[1]]$cum_treated[, 2, 2]), 3) ==
                              round(diff(y[[1]]$N[, 2, 3]), 3)))

            diff(y[[1]]$cum_treated[, 2, 1])
            diff(rowSums(y[[1]]$N[, 2, 2:3]))

          })

test_that("the number of person-years exposed is as expected", {
 
  # Previous method of aggregating over UIAS or U to get time individuals
  # spent exposed, underestimated pyes as values depended on how often
  # the model was output. Increasing frequency of output improved estimates.
  # The new ODE in odin model code should get around this by summing
  # the pye spent in UIAS or U over continuous time, and can be output directly
  # from the model result as 'cum_pye_trial_pov' or 'cum_pye_true'
  
  gp <- gono_params_trial(1)[1]
  n_erlang <- 1
  n_diag_rec <- 2
  N <- 600
  
  idx <- stratum_index_xvw_trial(n_erlang = n_erlang, n_diag_rec = n_diag_rec)
  idx$never_diag <- seq(idx$V[1], by = n_diag_rec, length.out = n_erlang + 1)

  # run
  # output every 1 year 
  tt <- seq(0, 5, 1)
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e100000000,
                            vea = 0, vei = 0, ved = 0, ves = 0,
                            n_erlang = n_erlang,
                            stochastic = FALSE,
                            n_diag_rec = n_diag_rec, N = N)
  
  # output every 20th of a year 
  inc <- 0.05
  tt_2 <- seq(0, 5, inc)
  y_2 <- run_onevax_xvw_trial(tt = tt_2, gp, dur = 1e100000000,
                              vea = 0, vei = 0, ved = 0, ves = 0,
                              n_erlang = n_erlang,
                              stochastic = FALSE,
                              n_diag_rec = n_diag_rec, N = N)
  
  # expect the cumulative person years calculated by odin to be the same
  # regardless of how often we output 
  
  # for U I A S pye
  cum_pye_trial_pov_odin <- y[[1]]$cum_pye_trial_pov[length(tt), 2, idx$X[1]]
  cum_pye_trial_pov_odin_2 <- y_2[[1]]$cum_pye_trial_pov[length(tt_2),
                                                         2, idx$X[1]]
  
  # for U only pye
  cum_pye_true_odin <- y[[1]]$cum_pye_true[length(tt), 2, idx$X[1]]
  cum_pye_true_odin_2 <- y_2[[1]]$cum_pye_true[length(tt_2), 2, idx$X[1]]
  
  # compare less frequent with more frequent outputting:
  expect_equal(cum_pye_trial_pov_odin,
               cum_pye_trial_pov_odin_2, tolerance = 1e-5)
  expect_equal(cum_pye_true_odin,
               cum_pye_true_odin_2, tolerance = 1e-5)
  
  # expect cumulative aggregated pyes (old mehtod) to be smaller than those
  # calculated by odin as they will have missed some person-time between outputs 
  
  # for UIAS pye
  pye_trial_pov_agg <- t(aggregate(y[1], "N", stratum = idx$X[1]))[, -1]
  cum_pye_trial_pov_agg <- cumsum(pye_trial_pov_agg)[length(tt) - 1]
  
  # for U only pye
  pye_true_agg <- t(aggregate(y[1], "U", stratum = idx$X[1]))[, -1]
  cum_pye_true_agg <- cumsum(pye_true_agg)[length(tt) - 1]
  
  expect_true(cum_pye_trial_pov_agg < cum_pye_trial_pov_odin)
  expect_true(cum_pye_true_agg < cum_pye_true_odin)
  
  # expect aggregated pyes to be greater when outputs are more frequent 
  # note: outputs here divided by 1/inc because individuals are only
  # contributing 1/inc'th of person-time (not a year)
  
  # for UIAS pye
  pye_trial_pov_agg_2 <- (t(aggregate(y_2[1],
                                      "N", stratum = idx$X[1]))[, -1])/(1/inc)
  cum_pye_trial_pov_agg_2 <- cumsum(pye_trial_pov_agg_2)[length(tt_2) - 1]
  
  # for U only pye
  pye_true_agg_2 <- (t(aggregate(y_2[1],
                                 "U", stratum = idx$X[1]))[, -1])/(1/inc)
  cum_pye_true_agg_2 <- cumsum(pye_true_agg_2)[length(tt_2) - 1]
  
  expect_true(cum_pye_trial_pov_agg_2 > cum_pye_trial_pov_agg)
  expect_true(cum_pye_true_agg_2 > cum_pye_true_agg)
  
  # expect pye calculated through aggregation over UIAS/U when output more
  # frequently to tend to pyes calculated through odin 
  # as a check that odin output is closer to the 'actual' amount of pyes exposed
  
  expect_true(cum_pye_trial_pov_odin - cum_pye_trial_pov_agg > 
                cum_pye_trial_pov_odin - cum_pye_trial_pov_agg_2)
  
  expect_true(cum_pye_true_odin - cum_pye_true_agg > 
                cum_pye_true_odin - cum_pye_true_agg_2) 
  
  # expect pye to be greater in the vaccinated > unvaccinated AND
  # expect pye to increase by N/2 each timepoint when vaccination perfect

  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e100000000,
                            vea = 1, vei = 0, ved = 0, ves = 0,
                            n_erlang = n_erlang,
                            stochastic = FALSE,
                            n_diag_rec = n_diag_rec, N = N)

   x  <- (aggregate(y, "cum_pye_trial_pov", stratum = idx$X[1]))[-1]
   vw <- (aggregate(y, "cum_pye_trial_pov", stratum = idx$never_diag))[-1]

  expect_true(all(vw > x))
  expect_equal(diff(vw), rep(300, 4))


})
