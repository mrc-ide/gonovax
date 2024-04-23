######### testing #######
context("stochastic trial model check")

# the same as test-model_trial_stochastic but for stochastic = TRUE, plus some
# extra tests

test_that("no new infections if lambda is 0", {

  # no vaccination
  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])
  params$lambda <- 0
  mod <- model_trial_stochastic$new(user = params,
                                    unused_user_action = "ignore")
  tt <- seq.int(0, 1)
  y <- mod$run(step = tt)
  y <- mod$transform_variables(y)


  expect_true(all(y$I == 0))             #incubating remain at 0
  expect_true(all(y$S == 0))             #no symptomatic infections
  expect_true(all(y$cum_incid == 0))
  expect_true(all(unlist(y) >= 0))                     #no compartments negative
  expect_true(all(round(rowSums(y$N)) == params$N0))   #number of participants
  #doesn't change
  # presence of vaccination

  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]],
                               vax_params = vax_params_xvw_trial(
                                 stochastic = TRUE
                               ), p_v = 0.5)
  tt <- seq.int(0, 1)
  params$lambda <- 0
  mod2 <- model_trial_stochastic$new(user = params,
                                     unused_user_action = "ignore")
  y <- mod2$run(step = tt)

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
  mod <- model_trial_stochastic$new(user = params,
                                    unused_user_action = "ignore")
  tt <- seq.int(0, 1)
  y <- mod$run(step = tt)
  y <- mod$transform_variables(y)

  expect_true(all(unlist(y) >= 0))
  expect_true(all(round(rowSums(y$N)) == params$N0))

})

test_that("all individuals are in the high activity group", {

  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])
  mod <- model_trial_stochastic$new(user = params,
                                    unused_user_action = "ignore")
  tt <- seq.int(0, 1)
  y <- mod$run(step = tt)
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
  mod <- model_trial_stochastic$new(user = params,
                                    unused_user_action = "ignore")
  tt <- seq.int(0, 1)
  y <- mod$run(step = tt)
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
  mod <- model_trial_stochastic$new(user = params,
                                    unused_user_action = "ignore")
  tt <- seq.int(0, 1)
  y <- mod$run(step = tt)
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
  mod <- model_trial_stochastic$new(user = params,
                                    unused_user_action = "ignore")
  tt <- seq.int(0, 1)
  y <- mod$run(step = tt)
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
  params$mu <- params$eta <- 0
  mod <- model_trial_stochastic$new(user = params,
                                    unused_user_action = "ignore")

  tt <- seq.int(0, 1)
  y <- mod$run(tt)
  y <- mod$transform_variables(y)
  expect_true(all(y$T == 0))
  expect_true(all(y$cum_treated == 0))

  expect_true(all(unlist(y) >= 0))
  expect_true(all(round(rowSums(y$N)) == params$N0))

})


test_that("Model works with vaccination and waning", {
  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, 1)
  set.seed(1)
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0, vei = 0, ved = 0, ves = 0,
                            stochastic = TRUE)

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
  tt <- seq.int(0, 1)
  set.seed(1)
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 1, vei = 0, ved = 0, ves = 0,
                            stochastic = TRUE)

  # VEa = 1, no infections in V, X >0, W >0 (after t=0)

  expect_true(all(y[[1]]$cum_incid[, 2, 2] == 0))
  expect_true(all(y[[1]]$cum_incid[2:2, 2, 1] > 0))
  expect_true(y[[1]]$cum_incid[2, 2, 3] > 0)

  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])

  expect_true(all(unlist(y) >= 0))
  expect_true(all(round(rowSums(y[[1]]$N)) == params$N0))

  # VEa = 0, infections in X = V + W
  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, 1)


  set.seed(1)
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0, vei = 0, ved = 0, ves = 0,
                            stochastic = TRUE)

  x <- round(as.numeric(y[[1]]$cum_incid[2, 2, 1]), -2)
  vw <- round(as.numeric(y[[1]]$cum_incid[2, 2, 2]) +
                as.numeric(y[[1]]$cum_incid[2, 2, 3]), -2)


  expect_true(all(unlist(y) >= 0))
  expect_true(all(round(rowSums(y[[1]]$N)) == params$N0))

})


test_that("VEs behaves as expected ", {
  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, 1)
  set.seed(1)
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0, vei = 0, ved = 0, ves = 1,
                            stochastic = TRUE)


  # VEs = 1 symptomatic diagnoses in V = 0, but for X & W > 0
  expect_true(all(y[[1]]$cum_diag_s[, 2, 2] == 0))
  expect_true(y[[1]]$cum_diag_s[2, 2, 1] > 0)
  expect_true(y[[1]]$cum_diag_s[2, 2, 3] > 0)

  # VEs = 1 asymptomatic diagnoses in V > 0
  expect_true(y[[1]]$cum_diag_a[2, 2, 2] > 0)

  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])

  expect_true(all(unlist(y) >= 0))
  expect_true(all(round(rowSums(y[[1]]$N)) == params$N0))

})


test_that("VEd behaves as expected ", {
  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, 1)
  set.seed(1)
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e1000000000,
                            vea = 0, vei = 0, ved = 1, ves = 0,
                            stochastic = TRUE)

  # VEd = 1 cumulative incid in V > X (when waning doesn't occur!)
  # running for a year rather than 5 days so enough timepoints
  # to see if there is a difference

  vw <- y[[1]]$cum_incid[2, 2, 2] + y[[1]]$cum_incid[2, 2, 3]

  x <- y[[1]]$cum_incid[2, 2, 1]


  expect_true(vw > x)

  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])

  expect_true(all(unlist(y) >= 0))
  expect_true(all(round(rowSums(y[[1]]$N)) == params$N0))

  ###########
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

  y2 <- run_onevax_xvw_trial(tt = seq(0, 5), gono_params = gp,
                             initial_params_trial = list(init_params),
                             vea = 0, vei = 0, ved = 1, ves = 0,
                             dur = 1e99, stochastic = TRUE)

  expect_equal(y2[[1]]$A[, , 2], y2[[1]]$S[, , 2]) # 2 = V stratum

})


test_that("VEi behaves as expected ", {
  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, 1)
  set.seed(1)
  y0 <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                             vea = 0, vei = 0, ved = 0, ves = 0,
                             stochastic = TRUE)

  set.seed(1)
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0, vei = 1, ved = 0, ves = 0,
                            stochastic = TRUE)

  # Incidence in V is the same for VEi = 0 and VEi = 1

  expect_equal(y0[[1]]$cum_incid[2, 2, 2],
               y[[1]]$cum_incid[2, 2, 2])

  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])

  expect_true(all(unlist(y) >= 0))
  expect_true(all(round(rowSums(y[[1]]$N)) == params$N0))

  expect_true(all(unlist(y0) >= 0))
  expect_true(all(round(rowSums(y0[[1]]$N)) == params$N0))

})


test_that("aggregated time series output correctly", {
  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])
  mod <- model_trial_stochastic$new(user = params,
                                    unused_user_action = "ignore")
  tt <- seq.int(0, 1)
  y <- mod$run(tt)
  y <- mod$transform_variables(y)
  expect_equal(y$tot_treated, apply(y$cum_treated, 1, sum))
  expect_equal(y$tot_attended, apply(y$cum_screened, 1, sum) + y$tot_treated)

  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])

  expect_true(all(unlist(y) >= 0))
  expect_true(all(round(rowSums(y$N)) == params$N0))

})

test_that("expected cumulative outputs are cumulative", {

  years <- 5 * 365

  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, years)
  set.seed(1)
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0.5, vei = 0.5, ved = 0.5, ves = 0.5,
                            stochastic = TRUE)

  expect_true(all(diff(y[[1]]$cum_diag_a[c(years - 10, years + 1), 2, ]) > 0))
  expect_true(all(diff(y[[1]]$cum_incid[c(years - 10, years + 1), 2, ]) > 0))
  expect_true(all(diff(y[[1]]$cum_diag_s[c(years - 10, years + 1), 2, ]) > 0))
  expect_true(all(diff(y[[1]]$cum_treated[c(years - 10, years + 1), 2, ]) > 0))
  expect_true(all(diff(y[[1]]$cum_screened[c(years - 10, years + 1), 2, ]) > 0))

  set.seed(1)
  n_diag_rec <- 2
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0.5, vei = 0.5, ved = 0.5, ves = 0.5,
                            stochastic = TRUE, n_diag_rec = n_diag_rec)

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


test_that("n_erlang = n is working as expected", {
  #when n_erlang = n, nvax = n + 2, as generated in stratum_index_xvw_trial
  #(for n_diag_rec = 1 only)
  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, 1)
  set.seed(1)
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0.5, vei = 0.5, ved = 0.5, ves = 0.5,
                            stochastic = TRUE)

  expect_equal(dim(y[[1]]$N)[3], 3)

  n_erlang <- 2
  set.seed(1)
  y2 <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                             vea = 0.5, vei = 0.5, ved = 0.5, ves = 0.5,
                             n_erlang = n_erlang,
                             stochastic = TRUE)

  expect_equal(dim(y2[[1]]$N)[3], n_erlang + 2)
  idx <- stratum_index_xvw_trial(n_erlang)
  expect_equal(dim(y2[[1]]$N)[3], idx$n_vax)

  n_erlang <- 5
  set.seed(1)
  y3 <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                             vea = 0.5, vei = 0.5, ved = 0.5, ves = 0.5,
                             n_erlang = n_erlang,
                             stochastic = TRUE)
  expect_equal(dim(y3[[1]]$N)[3], n_erlang + 2)
  idx <- stratum_index_xvw_trial(n_erlang)
  expect_equal(dim(y3[[1]]$N)[3], idx$n_vax)

  #for n_diag_rec > 2
  n_erlang <- 2
  n_diag_rec <- 2
  tt <- seq.int(0, 1)
  set.seed(1)
  y2 <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                             vea = 0.5, vei = 0.5, ved = 0.5, ves = 0.5,
                             n_erlang = n_erlang,
                             stochastic = TRUE,
                             n_diag_rec = n_diag_rec)

  idx <- stratum_index_xvw_trial(n_erlang, n_diag_rec)
  expect_equal(idx$n_vax, 8)
  expect_equal(dim(y2[[1]]$N)[3], idx$n_vax)

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

  n_diag_rec <- 2
  n_erlang <- 2 #n_vax is 8
  idx <- stratum_index_xvw_trial(n_erlang, n_diag_rec)
  erlang_2_n_diag_rec_2 <- vax_params_xvw_trial(vea = 1, dur = dur,
                                                n_erlang = n_erlang,
                                                n_diag_rec = n_diag_rec)$w
  matrix_2_n_diag_rec_2 <- matrix(data = (c(rep(0, 18),
                                            -n_erlang / dur,
                                            rep(0, 8),
                                            -n_erlang / dur,
                                            rep(0, 6),
                                            n_erlang / dur, 0,
                                            -n_erlang / dur, rep(0, 6),
                                            n_erlang / dur, 0,
                                            -n_erlang / dur,
                                            rep(0, 6),
                                            n_erlang / dur,
                                            rep(0, 8),
                                            n_erlang / dur, 0, 0)),
                                  ncol = idx$n_vax,
                                  byrow = TRUE)
  expect_equal(erlang_2_n_diag_rec_2, matrix_2_n_diag_rec_2)

  # when vea is perfect, there are no infections in any of the
  # erlang strata

  n_erlang <- 3   #n_vax is 5
  idx <- stratum_index_xvw_trial(n_erlang)
  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, 1)
  set.seed(1)
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1,
                            vea = 1, vei = 0, ved = 0, ves = 0,
                            n_erlang = n_erlang,
                            stochastic = TRUE)

  expect_true(all(y[[1]]$cum_incid[, , 2] == 0)) #V1
  expect_true(all(y[[1]]$cum_incid[, , 3] == 0)) #V2
  expect_true(all(y[[1]]$cum_incid[, , 4] == 0)) #V3
  expect_true(all(y[[1]]$cum_incid[, , c(idx$V)] == 0)) #all V

  expect_true(all(y[[1]]$cum_incid[-1, 2, 1] != 0)) #X still has infections
  expect_true(all(y[[1]]$cum_incid[-1, 2, idx$X] != 0))
  expect_true(sum(y[[1]]$cum_incid[-1, 2, 5]) != 0) #W not under protection, has
  expect_true(sum(y[[1]]$cum_incid[-1, 2, idx$W]) != 0) # infections

  n_diag_rec <- 2
  n_erlang <- 3   #n_vax is 10
  idx <- stratum_index_xvw_trial(n_erlang, n_diag_rec)
  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, 1)
  set.seed(1)
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1,
                            vea = 1, vei = 0, ved = 0, ves = 0,
                            n_erlang = n_erlang,
                            stochastic = TRUE,
                            n_diag_rec = n_diag_rec)

  expect_true(all(y[[1]]$cum_incid[, , c(idx$V)] == 0)) #all of V no infections
  expect_true(all(y[[1]]$cum_incid[-1, 2, idx$X] != 0)) #X still has infections
  expect_true(sum(y[[1]]$cum_incid[-1, 2, idx$W]) != 0) #W not under protection,
  #has infections

})

test_that("create_vax_map_branching is appropriately generating the
          diagnosis history mapping (diag_rec)", {
            #default is n_erlang = 1, n_diag_rec = 1
            #expect mapping to not move anyone upon diagnosis
            n_erlang <- 1
            n_diag_rec <- 1
            #copied from vax_params_xvw_trial code:
            idx <- stratum_index_xvw_trial(n_erlang)

            # diagnosed from
            i_eligible <- seq_len(idx$n_vax)[seq_len(idx$n_vax) %%
                                               n_diag_rec != 0]

            # diagnosed to
            i_p <- seq_len(idx$n_vax)[seq_len(idx$n_vax) %% n_diag_rec != 1]

            # create diagnosis history mapping
            diag_rec <- create_vax_map_branching(idx$n_vax, c(0, 1), i_eligible,
                                                 i_p, set_vbe = FALSE, idx)

            expect_equal(idx$n_vax, 3)
            expect_equal(diag_rec, array(0, dim = c(2, idx$n_vax, idx$n_vax)))

            #n_erlang = 2, n_diag_rec = 2

            n_erlang <- 2
            n_diag_rec <- 2

            #copied from vax_params_xvw_trial code:
            idx <- stratum_index_xvw_trial(n_erlang, n_diag_rec)

            # diagnosed from
            i_eligible <- seq_len(idx$n_vax)[seq_len(idx$n_vax) %%
                                               n_diag_rec != 0]

            # diagnosed to
            i_p <- seq_len(idx$n_vax)[seq_len(idx$n_vax) %% n_diag_rec != 1]

            # create diagnosis history mapping
            diag_rec <- create_vax_map_branching(idx$n_vax, c(0, 1), i_eligible,
                                                 i_p, set_vbe = FALSE, idx)

            #
            expect_equal(diag_rec, array(c(0, 1, 0, -1, rep(0, idx$n_vax * 4),
                                           0, 1, 0, -1, rep(0, idx$n_vax * 4),
                                           0, 1, 0, -1, rep(0, idx$n_vax * 4),
                                           0, 1, 0, -1, rep(0, idx$n_vax * 2)),
                                         dim = c(2, idx$n_vax, idx$n_vax)))

          })

test_that("stochasticity has been incorporated", {

  gp <- gono_params_trial(1)[1]
  n_erlang <- 1
  tt <- seq.int(0, 10)
  set.seed(1)
  y1 <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                             vea = 0, vei = 0, ved = 0, ves = 0,
                             n_erlang = n_erlang,
                             stochastic = TRUE)
  set.seed(2)
  y2 <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                             vea = 0, vei = 0, ved = 0, ves = 0,
                             n_erlang = n_erlang,
                             stochastic = TRUE)


  #model runs with same parameter sets and conditions give different numbers
  # of diagnoses due to stochasticity
  expect_true(y1[[1]]$cum_diag_s[11, 2, 1] != y2[[1]]$cum_diag_s[11, 2, 1])
  expect_true(y1[[1]]$cum_diag_s[11, 2, 2] != y2[[1]]$cum_diag_s[11, 2, 2])

  expect_true(y1[[1]]$cum_diag_a[11, 2, 1] != y2[[1]]$cum_diag_a[11, 2, 1])
  expect_true(y1[[1]]$cum_diag_a[11, 2, 2] != y2[[1]]$cum_diag_a[11, 2, 2])

  # same with total treated and total attending care
  expect_true(y1[[1]]$cum_treated[11, 2, 1] != y2[[1]]$cum_treated[11, 2, 1])
  expect_true(y1[[1]]$cum_treated[11, 2, 2] != y2[[1]]$cum_treated[11, 2, 2])

  expect_true(y1[[1]]$tot_attended[11] != y2[[1]]$tot_attended[11])

  # same with incidence
  expect_true(y1[[1]]$cum_incid[11, 2, 1] != y2[[1]]$cum_incid[11, 2, 1])
  expect_true(y1[[1]]$cum_incid[11, 2, 2] != y2[[1]]$cum_incid[11, 2, 2])

  #when FOI = 0, screening is still different between stochastic runs even
  #when U[i, j] is the same
  gp[[1]]$lambda <- 0
  set.seed(1)
  y1 <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e10000,
                             vea = 0, vei = 0, ved = 0, ves = 0,
                             n_erlang = n_erlang,
                             stochastic = TRUE)
  set.seed(2)
  y2 <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e10000,
                             vea = 0, vei = 0, ved = 0, ves = 0,
                             n_erlang = n_erlang,
                             stochastic = TRUE)

  expect_true(y1[[1]]$cum_screened[11, 2, 1] != y2[[1]]$cum_screened[11, 2, 1])
  expect_true(y1[[1]]$cum_screened[11, 2, 2] != y2[[1]]$cum_screened[11, 2, 2])


})

test_that("waning parameters are being generated correctly", {

  vax <- vax_params_xvw_trial(vea = 0.5, vei = 0.5, ved = 0.5, ves = 0.5,
                              dur = 20, n_erlang = 1, stochastic = TRUE)

  # When stochastic = TRUE ...
  # the waning matrix should only be made up of 1, -1 and 0 (even though rate
  # is a draw made with rate 1/20)
  expect_equal(vax$w, matrix(c(rep(0, 4), -1, rep(0, 2), 1, 0),
                             nrow = 3, byrow = TRUE))

  #D should be a vector of length 3, the  diagonal of the waning matrix
  #but with -n_erlang/dur as the value for waning rate (at position 2)
  expect_equal(vax$D, c(0, -0.05, 0))

  #Edits made to vax_params_xvw_trial should not affect the expected waning maps
  # generated for the deterministic trial when stochastic = FALSE

  vax <- vax_params_xvw_trial(vea = 0.5, vei = 0.5, ved = 0.5, ves = 0.5,
                              dur = 20, n_erlang = 1, stochastic = FALSE)
  expect_equal(vax$w, matrix(c(rep(0, 4), -0.05, rep(0, 2), 0.05, 0),
                             nrow = 3, byrow = TRUE))

  #D is present as the diagonal of w (will not be used)
  expect_equal(vax$D, c(0, -0.05, 0))


})



test_that("correct number of individuals are set up in each trial arm", {

  gp <- gono_params_trial(1)[1]
  n_erlang <- 1
  tt <- seq.int(0, 1)
  set.seed(1)

  N <- 500

  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0, vei = 0, ved = 0, ves = 0,
                            n_erlang = n_erlang,
                            stochastic = TRUE,
                            N = N)

  expect_equal(y[[1]]$N[1, 2, 1], N / 2)
  expect_equal(y[[1]]$N[1, 2, 2], N / 2)
  expect_equal(y[[1]]$N[1, 2, 3], 0)

  # for n_diag_rec > 1

  n_diag_rec <- 3
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0, vei = 0, ved = 0, ves = 0,
                            n_erlang = n_erlang,
                            stochastic = TRUE,
                            N = N, n_diag_rec = n_diag_rec)

  # N/2 in X.I and V1.I for t = 0
  expect_true(all(y[[1]]$N[1, 2, c(1, 4)] == N / 2))

  # 0 elsewhere
  expect_true(all(y[[1]]$N[1, 2, c(2, 3, 5:9)] == 0))

})


test_that("model can be output more frequently than each year", {

  #NB stochastic model always runs in days but tt is supplied in year form
  # e.g c(0, 1, 2) same as in the deterministic trial, but gets converted to day
  # form under the hood e.g c(0, 365, 730)
  # testing this and that we can choose to output every 6 months or every
  # quarter

  gp <- gono_params_trial(1)[1]
  n_erlang <- 1
  tt <- seq.int(0, 10)
  set.seed(1)

  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0, vei = 0, ved = 0, ves = 0,
                            n_erlang = n_erlang,
                            stochastic = TRUE)

  #values of 'time' are t/365 to give the years
  #values of 't' are increments of 365
  #time should be whole numbers from 0 to 10

  expect_equal(y[[1]]$t, seq(0, 3650, 365))
  expect_equal(y[[1]]$time, seq(0, 10))


  #every 6 months
  tt <- seq(0, 2, 1 / 2)
  set.seed(1)

  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0, vei = 0, ved = 0, ves = 0,
                            n_erlang = n_erlang,
                            stochastic = TRUE)

  expect_equal(y[[1]]$t, round(seq(0, 730, 365 / 2)))
  expect_equal(round(y[[1]]$time, 2), seq(0, 2, 1 / 2))

  #every 3 months
  tt <- seq(0, 2, 1 / 4)
  set.seed(1)
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0, vei = 0, ved = 0, ves = 0,
                            n_erlang = n_erlang,
                            stochastic = TRUE)

  expect_equal(y[[1]]$t, round(seq(0, 730, 365 / 4)))
  expect_equal(round(y[[1]]$time, 2), seq(0, 2, 1 / 4))

})

test_that("for n_diag_rec > 1, total N summed over X or V+W is
          the same and correct", {

            gp <- gono_params_trial(1)[1]
            n_erlang <- 1
            tt <- seq.int(0, 5)
            set.seed(1)
            n_diag_rec <- 3
            N <- 6e+05

            y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e03,
                                      vea = 0, vei = 0, ved = 0, ves = 0,
                                      n_erlang = n_erlang,
                                      stochastic = TRUE,
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
                                      stochastic = TRUE,
                                      n_diag_rec = n_diag_rec, N = N)

            #X.I, X.II, X.III
            expect_equal(rowSums(y[[1]]$N[, 2, 1:3]), rep(N / 2, length(tt)))
            expect_equal(rowSums(y[[1]]$N[, 2, 4:12]), rep(N / 2, length(tt)))

          })

test_that("for n_diag_rec > 1,  the number diagnosed = the number recorded", {

  gp <- gono_params_trial(1)[1]
  n_erlang <- 1
  tt <- seq.int(0, 5)
  set.seed(1)
  n_diag_rec <- 2
  N <- 6e+05

  #no waning for simplicity
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e100000000,
                            vea = 0, vei = 0, ved = 0, ves = 0,
                            n_erlang = n_erlang,
                            stochastic = TRUE,
                            n_diag_rec = n_diag_rec, N = N)

  # number diagnosed in '.I' becomes the number in '.II'
  expect_true(all(round(y[[1]]$cum_diag_a[, 2, 1] +
                          y[[1]]$cum_diag_s[, 2, 1]) ==
                    round(y[[1]]$N[, 2, 2])))
  expect_true(all(round(y[[1]]$cum_diag_a[, 2, 3] +
                          y[[1]]$cum_diag_s[, 2, 3]) ==
                    round(y[[1]]$N[, 2, 4])))
  expect_true(all(round(y[[1]]$cum_diag_a[, 2, 5] +
                          y[[1]]$cum_diag_s[, 2, 5]) ==
                    round(y[[1]]$N[, 2, 6])))

  n_diag_rec <- 3

  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e100000000,
                            vea = 0, vei = 0, ved = 0, ves = 0,
                            n_erlang = n_erlang,
                            stochastic = TRUE,
                            n_diag_rec = n_diag_rec, N = N)

  # number diagnosed each year in '.I' becomes the N gained in '.II' &
  # '.III' diagnosis history strata for that year
  # (because some people could get diagnosed twice in one year)
  # and number treated in '.II' becomes the number in '.III' diagnosis
  # history
  expect_true(all(round(diff(y[[1]]$cum_diag_a[, 2, 1] +
                               y[[1]]$cum_diag_s[, 2, 1])) ==
                    round(diff(rowSums(y[[1]]$N[, 2, 2:3])))))
  expect_true(all(round(diff(y[[1]]$cum_diag_a[, 2, 2] +
                               y[[1]]$cum_diag_s[, 2, 2])) ==
                    round(diff(y[[1]]$N[, 2, 3]))))

  #as movement into next diagnosis history stratum occurs upon those
  #entering T, expect T in .I to be empty for all time
  #(recorded move to T in .II and so on so we don't expect these to be
  #empty)

  expect_true(all(y[[1]]$T[, , c(1, 4, 7)] == 0))

})

test_that("for n_diag_rec > 1, when lambda = 0,
          diagnosis history strata > '.I' are empty", {

            gp <- gono_params_trial(1)[1]
            gp[[1]]$lambda <- 0
            n_erlang <- 1
            tt <- seq.int(0, 5)
            set.seed(1)
            n_diag_rec <- 3

            y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e1000000,
                                      vea = 0, vei = 0, ved = 0, ves = 0,
                                      n_erlang = n_erlang,
                                      stochastic = TRUE,
                                      n_diag_rec = n_diag_rec)

            # no one has been diagnosed = no diagnosis history movement
            expect_true(all(rowSums(y[[1]]$N[, 2, c(2, 3, 5, 6, 8, 9)]) == 0))
            expect_true(all(rowSums(y[[1]]$N[, 2, c(1, 4)]) == 6e+05))

            # and no one in the low activity group
            expect_true(all(rowSums(y[[1]]$N[, 1, ]) == 0))

          })


test_that("no asymptomatic diagnoses recorded when asymp_recorded = FALSE", {

  gp <- gono_params_trial(1)[1]
  n_erlang <- 1
  tt <- seq.int(0, 5)
  n_diag_rec <- 2
  N <- 6e+05

  #no waning for simplicity
  y_symp_only <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e100000000,
                                      vea = 0, vei = 0, ved = 0, ves = 0,
                                      n_erlang = n_erlang,
                                      stochastic = TRUE,
                                      n_diag_rec = n_diag_rec, N = N,
                                      asymp_recorded = FALSE)

  y_symp_asymp <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e100000000,
                                       vea = 0, vei = 0, ved = 0, ves = 0,
                                       n_erlang = n_erlang,
                                       stochastic = TRUE,
                                       n_diag_rec = n_diag_rec, N = N)
  #(asymp_recorded defaults to TRUE)

  #when asymptomatic diagnoses not recorded, fewer individuals move to the
  #next diagnoses history stratum (.II) so stratum X.I & V.I will have larger
  #N for asymp_recorded = FALSE
  expect_true(all(y_symp_only[[1]]$N[-1, 2, c(1, 3)] >
                    y_symp_asymp[[1]]$N[-1, 2, c(1, 3)]))

  #asymptomatic individuals move to treatment but are not recorded in the trial
  #i.e number diagnosed in .I is not the same as the number gained in .II
  #when asymp_recorded = FALSE
  expect_true(all(round(y_symp_only[[1]]$cum_diag_a[-1, 2, 1] +
                          y_symp_only[[1]]$cum_diag_s[-1, 2, 1]) !=
                    round(y_symp_only[[1]]$N[-1, 2, 2])))
  expect_true(all(round(y_symp_only[[1]]$cum_diag_a[-1, 2, 3] +
                          y_symp_only[[1]]$cum_diag_s[-1, 2, 3]) !=
                    round(y_symp_only[[1]]$N[-1, 2, 4])))

  #T.I is no longer empty when asymp_recorded  = FALSE
  expect_true(all(y_symp_asymp[[1]]$T[-1, 2, c(1, 3)] == 0))
  expect_true(all(y_symp_only[[1]]$T[-1, 2, c(1, 3)] != 0))

  #number asymptomatic diagnoses = number entering T.I
  #(no recovery to U so we can keep track of entrance into T)
  gp[[1]][["rho"]] <- 0
  y_symp_only_ze_rho <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e100000000,
                                             vea = 0, vei = 0, ved = 0, ves = 0,
                                             n_erlang = n_erlang,
                                             stochastic = TRUE,
                                             n_diag_rec = n_diag_rec, N = N,
                                             asymp_recorded = FALSE)

  expect_true(all(round(diff(y_symp_only_ze_rho[[1]]$cum_diag_a[, 2, 1])) ==
                    round(diff(y_symp_only_ze_rho[[1]]$T[, 2, 1]))))

  #number symptomatic diagnoses = number entering .II
  expect_true(all(round(y_symp_only_ze_rho[[1]]$cum_diag_s[, 2, 1])
                  == round(y_symp_only_ze_rho[[1]]$N[, 2, 2])))

})
