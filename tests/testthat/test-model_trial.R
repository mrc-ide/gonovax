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

  # VEd = 1 rate of natural clearance (nu) is equal to the rate of care seeking
  # (mu)
  # if same number of people start in A as in S and parameters
  # lambda and eta = 0 (i.e other routes into A & S and out of A respectively)
  # the rate of movement of people A -> U and S -> T in 'V' should be the same

  # replace lambda, eta = 0
  elements_to_replace <- c("eta", "lambda")
  for (element in elements_to_replace) {
    gp[[1]][[element]] <- 0
  }

  # set starting conditions
  pars <- c(demographic_params_trial(N = 600), gp)
  n_vax <- 3
  U0 <- I0 <- A0 <- S0 <- T0 <- array(0, c(2, n_vax))

  #half of population in A0 and half in S0
  #everyone begins vaccinated
  N0 <- pars$N0 * outer(pars$q, c(0, 1, 0))
  A0 <- pars$N0 * outer(pars$q, c(0, 0.5, 0))
  S0 <- pars$N0 * outer(pars$q, c(0, 0.5, 0))

  init_params <- list(U0 = U0, I0 = I0, A0 = A0, S0 = S0, T0 = T0)

  y2 <- run_onevax_xvw_trial(tt = tt, gono_params = gp,
                             initial_params_trial = list(init_params),
                             vea = 0, vei = 0, ved = 1, ves = 0,
                             dur = 1e99)

  expect_equal(y2[[1]]$A[, , 2], y2[[1]]$S[, , 2]) # 2 = V stratum

  # if ved = 0.5, the effective rate of clearance is about double the rate of
  # natural clearance
  y_0 <- run_onevax_xvw_trial(tt = tt, gono_params = gp,
                              initial_params_trial = list(init_params),
                              vea = 0, vei = 0, ved = 0, ves = 0,
                              dur = 1e99)
  y_0.5 <- run_onevax_xvw_trial(tt = tt, gono_params = gp,
                                initial_params_trial = list(init_params),
                                vea = 0, vei = 0, ved = 0.5, ves = 0,
                                dur = 1e99)

  expect_equal(tolerance = 1, 2 * abs(round(diff(y_0[[1]]$A[, 2, 2]), 0)),
               abs(round(diff(y_0.5[[1]]$A[, 2, 2]), 0)))

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

  expect_equal(tolerance = 0.1,
               y[[1]]$U, array(c(0, 0, 0, 0, 0, 0, 300000.0, 298770.0,
                                 297547.4, 296334.9, 295135.1, 293950.3,
                                 0, 0, 0,
                                 0, 0, 0, 300000.0, 299386.7, 298791.9,
                                 298227.1, 297698.8, 297210.3, 0, 0, 0, 0,
                                 0, 0, 0, 0.8193902, 1.6337723, 2.4432536,
                                 3.2479778, 4.0481123),
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

test_that("can set n_AU conditional statement works as expected", {

#when mu = 0, defaults to nu
  #and when ved = 0, and mu not = 0, mu's cancel out, and rate is nu
  #so in both cases movement A->U should be the same
  tt <- seq.int(0, 5) / 365

  gp_mu_og <- gono_params_trial(1)[1]
  gp_mu_zero <- gono_params_trial(1)[1]

  # replace parameters with zero
  elements_to_replace <- c("eta", "lambda", "mu")
  for (element in elements_to_replace) {
    gp_mu_zero[[1]][[element]] <- 0
  }
  elements_to_replace <- c("eta", "lambda")
  for (element in elements_to_replace) {
    gp_mu_og[[1]][[element]] <- 0
  }

  # set starting conditions
  pars <- c(demographic_params_trial(N = 600), gp)
  n_vax <- 3
  U0 <- I0 <- A0 <- S0 <- T0 <- array(0, c(2, n_vax))

  #all of population in A0 because we are only interested in movement out of A
  N0 <- pars$N0 * outer(pars$q, c(0, 1, 0))
  A0 <- pars$N0 * outer(pars$q, c(0, 1, 0))

  init_params <- list(U0 = U0, I0 = I0, A0 = A0, S0 = S0, T0 = T0)

  y1 <- run_onevax_xvw_trial(tt = tt, gono_params = gp_mu_zero,
                             initial_params_trial = list(init_params),
                             vea = 0, vei = 0, ved = 0, ves = 0,
                             dur = 1e99)

  y2 <- run_onevax_xvw_trial(tt = tt, gono_params = gp_mu_og,
                            initial_params_trial = list(init_params),
                            vea = 0, vei = 0, ved = 0, ves = 0,
                            dur = 1e99)

  expect_equal(y1, y2)

  #adding in ved when mu = 0, makes no difference
  y3 <- run_onevax_xvw_trial(tt = tt, gono_params = gp_mu_zero,
                             initial_params_trial = list(init_params),
                             vea = 0, vei = 0, ved = 0.5, ves = 0,
                             dur = 1e99)
  expect_equal(y3, y2)

  #adding in ved when mu is not 0, does make a difference
  y4 <- run_onevax_xvw_trial(tt = tt, gono_params = gp_mu_og,
                             initial_params_trial = list(init_params),
                             vea = 0, vei = 0, ved = 0.5, ves = 0,
                             dur = 1e99)

  expect_error(expect_equal(y3, y4))

})
