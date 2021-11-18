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

  # presence of vaccination

  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]],
                               vax_params = vax_params_xvw_trial(),
                         coverage = 0.5)
  tt <- seq.int(0, 5) / 365
  params$lambda <- 0
  mod2 <- model_trial$new(user = params, unused_user_action = "ignore")
  y <- mod2$run(t = tt)

  y <- mod2$transform_variables(y)

  expect_true(all(y$I == 0))
  expect_true(all(y$S == 0))
  expect_true(all(y$cum_incid == 0))

})

test_that("number of individuals always >= 0 even at high lambda", {
  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])

  params$lambda <- 10
  mod <- model_trial$new(user = params, unused_user_action = "ignore")
  tt <- seq.int(0, 5) / 365
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)

  expect_true(all(unlist(y) >= 0))

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

})


test_that("there are no symptomatic infections when psi = 0", {
  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])
  params$psi <- 0
  mod <- model_trial$new(user = params, unused_user_action = "ignore")
  tt <- seq.int(0, 5) / 365
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  expect_true(any(y$I == 0))
  expect_true(all(y$S == 0))
  expect_true(all(y$cum_diag_s == 0))
  expect_true(all(y$cum_diag_a[-1, , ] >= 0))
  expect_true(all(unlist(y) >= 0))
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
  expect_true(any(y$I == 0))
  expect_true(all(y$A == 0))
  expect_true(all(y$cum_diag_a == 0))
  expect_true(all(y$cum_diag_s[-1, , ] >= 0))
  expect_true(all(unlist(y) >= 0))
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

})

test_that("no-one is treated when mu and eta = 0", {
  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])
  params$mu <- params$eta_h <- params$eta_l <-  0
  mod <- model_trial$new(user = params, unused_user_action = "ignore")

  tt <- seq.int(0, 5) / 365
  y <- mod$run(tt)
  y <- mod$transform_variables(y)
  expect_true(all(y$T == 0))
  expect_true(all(y$cum_treated == 0))
  expect_true(all(unlist(y) >= 0))
})


test_that("Model works with vaccination and waning", {
  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, 5) / 365
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                       vea = 0, vei = 0, ved = 0, ves = 0,
                       coverage = 0.5)

  #Total N in high group in X and  V is equal and non zero at t = 0
  expect_true(all(y[[1]]$N[1, 2, 1] == y[[1]]$N[1, 2, 2]))
  expect_true(all(y[[1]]$N[1, 2, 1] > 0))
  expect_true(all(y[[1]]$N[1, 2, 2] > 0))

  #Total N in waning is 0 at t =0 and increases over time
  expect_true(all(y[[1]]$N[1, 2, 3] == 0))
  length_tt <- length(y[[1]]$N[, 2, 3])
  expect_true(all(y[[1]]$N[2:length_tt, 2, 3] > 0))


})


test_that("VEa behaves as expected ", {
  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, 5) / 365
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 1, vei = 0, ved = 0, ves = 0,
                            coverage = 0.5)

  # VEa = 1, no infections in V, X >0, W >0 (after t=0)

  expect_true(all(y[[1]]$cum_incid[, 2, 2] == 0))
  expect_true(all(y[[1]]$cum_incid[2:6, 2, 1] > 0))
  expect_true(all(y[[1]]$cum_incid[2:6, 2, 3] > 0))

  # VEa = 0, infections in X = V + W
  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, 5) / 365
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0, vei = 0, ved = 0, ves = 0,
                            coverage = 0.5)


  x <- as.numeric(y[[1]]$cum_incid[6, 2, 1])
  vw <- as.numeric(y[[1]]$cum_incid[6, 2, 2]) +
    as.numeric(y[[1]]$cum_incid[6, 2, 3])
  expect_equal(x, vw)

})


test_that("VEs behaves as expected ", {
  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, 5) / 365
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0, vei = 0, ved = 0, ves = 1,
                            coverage = 0.5)


  # VEs = 1 symptomatic diagnoses in V = 0, but for X & W > 0
  expect_true(all(y[[1]]$cum_diag_s[, 2, 2] == 0))
  expect_true(all(y[[1]]$cum_diag_s[2:6, 2, 1] > 0))
  expect_true(all(y[[1]]$cum_diag_s[2:6, 2, 3] > 0))

  # VEs = 1 asymptomatic diagnoses in V > 0
  expect_true(all(y[[1]]$cum_diag_a[2:6, 2, 2] > 0))

})


test_that("VEd behaves as expected ", {
  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, 5) / 365
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0, vei = 0, ved = 1, ves = 0,
                            coverage = 0.5)

  # VEd = 1 cumulative incid in V > X + W

  v <- y[[1]]$cum_incid[6, 2, 2]
  xw <- y[[1]]$cum_incid[6, 2, 1] + y[[1]]$cum_incid[6, 2, 3]

  expect_true(v > xw)

})


test_that("VEi behaves as expected ", {
  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, 5) / 365
  y0 <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                             vea = 0, vei = 0, ved = 0, ves = 0,
                             coverage = 0.5)

  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                           vea = 0, vei = 1, ved = 0, ves = 0,
                           coverage = 0.5)

  # Incidence in V is the same for VEi = 0 and VEi = 1

  expect_equal(y0[[1]]$cum_incid[6, 2, 2],
               y[[1]]$cum_incid[6, 2, 2])

})


test_that("aggregated time series output correctly", {
  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])
  mod <- model_trial$new(user = params, unused_user_action = "ignore")
  tt <- seq.int(0, 5) / 365
  y <- mod$run(tt)
  y <- mod$transform_variables(y)
  expect_equal(y$tot_treated, apply(y$cum_treated, 1, sum))
  expect_equal(y$tot_attended, apply(y$cum_screened, 1, sum) + y$tot_treated)

})

