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
})


test_that("n_erlang = n is working as expected", {

  #when n_erlang = n, nvax = n + 2, as generated in stratum_index_xvw_trial
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

 # Edits made to vax_params_xvw_trial should not affect the expected waning maps
 # generated for the deterministic trial when stochastic = FALSE

  vax <- vax_params_xvw_trial(vea = 0.5, vei = 0.5, ved = 0.5, ves = 0.5,
                               dur = 20, n_erlang = 1, stochastic = FALSE)
   expect_equal(vax$w, matrix(c(rep(0, 4), -0.05, rep(0, 2), 0.05, 0),
                                nrow = 3, byrow = TRUE))

   #D is present as the diagonal of w (will not be used)
   expect_equal(vax$D, c(0, -0.05, 0))


})

test_that("model outputs for a basic run haven't changed between updates", {

gp <- gono_params_trial(1)[1]
n_erlang <- 1
tt <- seq.int(0, 1)
set.seed(1)
y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                           vea = 0, vei = 0, ved = 0, ves = 0,
                           n_erlang = n_erlang,
                           stochastic = TRUE)

expect_equal(y[[1]]$U,
    array(c(rep(0, 2), 300000, 249077, rep(0, 2), 300000, 248844,
            rep(0, 3), 237), dim = c(2, 2, 3)))

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

  expect_true(all(y[[1]]$N[, 2, 1] == N / 2))
  expect_true(all(y[[1]]$N[, 2, 2] == N / 2))
  expect_true(all(y[[1]]$N[, 2, 3] == 0))

})
