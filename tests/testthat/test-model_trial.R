######### testing #######
context("trial model check")

test_that("no new infections if lambda is 0", {

  # no vaccination
  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])   #no vaccination
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
  params$mu <- params$eta_h <- params$eta_l <-  0
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
                          297547.433626743, 296334.908051234, 295135.061477033,
                          293950.297077376, 0, 0, 0, 0, 0, 0, 3e+05,
                          299383.713797674, 298770.717502468, 298162.935800213,
                          297562.116361515, 296969.7952273, 0, 0, 0, 0, 0, 0, 0,
                          0.819388141002698, 1.63374216215416, 2.4431142416253,
                          3.24757521226149, 4.04721268804341),
                        dim = c(6, 2, 3)))
})


test_that("n_erlang = n is working as expected"){
  
  #when n_erlang = n, nvax = n + 2
  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, 5) / 365
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0.5, vei = 0.5, ved = 0.5, ves = 0.5)
  
  expect_equal(dim(y[[1]]$N)[3], 3)
  
  n_erlang = 2
  y2 <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0.5, vei = 0.5, ved = 0.5, ves = 0.5,
                            n_erlang = n_erlang)
  
  expect_equal(dim(y2[[1]]$N)[3], n_erlang + 2)
  
  n_erlang = 5
  y3 <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 0.5, vei = 0.5, ved = 0.5, ves = 0.5,
                            n_erlang = n_erlang)
  expect_equal(dim(y3[[1]]$N)[3], n_erlang + 2)
  
  #vax_params_xvw_trial generates the correct waning maps (in for n_erlang > 1)
  dur <- 1e03
  n_erlang <- 1
  erlang_1 <- vax_params_xvw_trial(vea = 1, dur = dur , n_erlang = n_erlang)$w
  matrix_1 <- matrix(data = (c(rep(0,4), -n_erlang/dur, rep(0,2),
                               n_erlang/dur, 0)), ncol = n_erlang + 2, byrow = TRUE)
  expect_equal(erlang_1, matrix_1)
  
  n_erlang <- 2
  erlang_2 <- vax_params_xvw_trial(vea = 1, dur = dur , n_erlang = n_erlang)$w
  matrix_2 <- matrix(data = (c(rep(0,5), -n_erlang/dur, rep(0,3),
                               n_erlang/dur, -n_erlang/dur, rep(0,3),
                               n_erlang/dur, 0)), ncol = n_erlang + 2, byrow = TRUE)
  expect_equal(erlang_2, matrix_2)
  
  # when vea is perfect, there are no infections in any of the
  # erlang strata
  
  n_erlang <- 3   #n_vax = 5
  gp <- gono_params_trial(1)[1]
  tt <- seq.int(0, 5) / 365
  y <- run_onevax_xvw_trial(tt = tt, gp, dur = 1e3,
                            vea = 1, vei = 0, ved = 0, ves = 0,
                            n_erlang = n_erlang)
  
  expect_true(all(y[[1]]$cum_incid[, , 2] == 0)) #V1
  expect_true(all(y[[1]]$cum_incid[, , 3] == 0)) #V2
  expect_true(all(y[[1]]$cum_incid[, , 4] == 0)) #V3
  
  expect_true(all(y[[1]]$cum_incid[-1, 2, 1] != 0)) #X still has infections  
  expect_true(all(y[[1]]$cum_incid[-1, 2, 5] != 0)) #W not under protection, has
                                                     # infections
  
}