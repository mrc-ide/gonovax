######### testing #######
context("trial model check")

test_that("no new infections if lambda is 0", {
  params <- gonovax::model_params_trial(gono_params_trial = gonovax::gono_params_trial(1)[[1]])                #gono_params() and model_params() generate params for testing
  params$lambda <- 0                                                             #set lambda to 0
  mod <- model_trial$new(user = params, unused_user_action = "ignore")     
  tt <- seq.int(0, 5) / 365                                             
  y <- mod$run(t = tt)                                                 
  y <- mod$transform_variables(y)                                                #turns into list
  
  
  expect_true(all(y$I ==0))             #incubating remain at 0
  expect_true(all(y$S ==0))             #no symptomatic infections
  expect_true(all(y$cum_incid == 0))
})

test_that("number of individuals always >= 0", {
  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]]) 
  
  params$lambda <- 2
  mod <- model_trial$new(user = params, unused_user_action = "ignore")     
  tt <- seq.int(0, 5) / 365                                             
  y <- mod$run(t = tt)                                                 
  y <- mod$transform_variables(y)
  
  expect_true(all(unlist(y) >= 0))
  
  params$lambda <- 3
  mod <- model_trial$new(user = params, unused_user_action = "ignore")     
  tt <- seq.int(0, 5) / 365                                             
  y <- mod$run(t = tt)                                                 
  y <- mod$transform_variables(y)
  
  expect_true(all(unlist(y) >= 0))
  
  params$lambda <- 4
  mod <- model_trial$new(user = params, unused_user_action = "ignore")     
  tt <- seq.int(0, 5) / 365                                             
  y <- mod$run(t = tt)                                                 
  y <- mod$transform_variables(y)
  
  expect_true(all(unlist(y) >= 0))
  
})

test_that("all individuals are in the high activity group", {

params <- gonovax::model_params_trial(gono_params_trial = gonovax::gono_params_trial(1)[[1]])
params$lambda <- 1.5                                                           
mod <- model_trial$new(user = params, unused_user_action = "ignore")     
tt <- seq.int(0, 5) / 365                                             
y <- mod$run(t = tt)                                                 
y <- mod$transform_variables(y)                                                
  
#low activity group empty
expect_true(all(y$U[ , 1, 1] ==0))
expect_true(all(y$I[ , 1, 1] ==0))
expect_true(all(y$A[ , 1, 1] ==0))
expect_true(all(y$S[ , 1, 1] ==0))
expect_true(all(y$T[ , 1, 1] ==0))
  
})

############ variatons of Lilith's tests:

test_that("there are no symptomatic infections when psi = 0", {
  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])
  params$lambda <- 1.5
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
  params$lambda <- 1.5
  params$psi <- 1                   #proportion symptomatic 
  params$S0[, ] <- params$A0[, ]    #moves starting individuals from A to S so still infections seeded
  params$A0[, ] <- 0
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
  params$lambda <- 1.5
  mod <- model_trial$new(user = params, unused_user_action = "ignore")
  tt <- seq.int(0, 5) / 365
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  
  expect_true(all(y$I[, , ][1,] == 0))
  expect_true(all(y$A[, , ][1,] == 0))
  expect_true(all(y$S[, , ][1,] == 0))
  expect_true(all(y$T[, , ][1,] == 0))
  expect_true(all(y$U[, , ][1,2] > 0))

})

test_that("no-one is treated when mu and eta = 0", {
  params <- model_params_trial(gono_params_trial = gono_params_trial(1)[[1]])
  params$lambda <- 1.5
  params$mu <- params$eta_h[] <- params$eta_l[] <-  0
  mod <- model_trial$new(user = params, unused_user_action = "ignore")
  
  tt <- seq.int(0, 5) / 365
  y <- mod$run(tt)
  y <- mod$transform_variables(y)
  expect_true(all(y$T == 0))
  expect_true(all(y$cum_treated == 0))
  expect_true(all(unlist(y) >= 0))
})


test_that("Bex model runs with no vaccination", {                                              #come back to this later 
  tt <- seq.int(0, 5) / 365
  params0 <- model_params(gono_params = gono_params(1)[[1]])
  mod0 <- model(user = params0, unused_user_action = "ignore")
  y0 <- mod0$run(t = tt)
  y0 <- mod0$transform_variables(y0)
  
  params1 <- model_params(gono_params = gono_params(1)[[1]],
                          vax_params = vax_params_xvwv(vbe = 0))
  mod1 <- model(user = params1, unused_user_action = "ignore")
  y1 <- mod1$run(t = tt)
  y1 <- mod1$transform_variables(y1)
  
  # check that nil vaccination gives same results as before
  expect_true(all(y1$U[, , 1, drop = FALSE] == y0$U))
  expect_true(all(y1$I[, , 1, drop = FALSE] == y0$I))
  expect_true(all(y1$A[, , 1, drop = FALSE] == y0$A))
  expect_true(all(y1$S[, , 1, drop = FALSE] == y0$S))
  expect_true(all(y1$T[, , 1, drop = FALSE] == y0$T))
  
  expect_true(all(y1$N[, , 2] == 0))
  expect_true(all(apply(y1$N, c(1, 2), sum) - 6e5 < 1e-6))
  
})


test_that("can initialise after time 0", {                                                    #come back to this 
  
  ## check with single parameter set
  params <- model_params(gono_params = gono_params(1)[[1]])
  mod <- model(user = params, unused_user_action = "ignore")
  
  tt <- seq.int(0, 5)
  y <- mod$run(tt)
  y <- mod$transform_variables(y)
  
  inits <- restart_params(y, n_vax = 1)
  
  expect_true(all(y$U[length(tt), , ] == inits$U0[, 1]))
  expect_true(all(y$I[length(tt), , ] == inits$I0[, 1]))
  expect_true(all(y$A[length(tt), , ] == inits$A0[, 1]))
  expect_true(all(y$S[length(tt), , ] == inits$S0[, 1]))
  expect_true(all(y$T[length(tt), , ] == inits$T0[, 1]))
  expect_true(5 == inits$t)
  
  ## check that restarting works properly
  tt1 <- seq.int(0, 10)
  y1 <- mod$run(tt1)
  y1 <- mod$transform_variables(y1)
  
  params2 <- model_params(gono_params = gono_params(1)[[1]],
                          init_params = inits)
  mod2 <- model(user = params2, unused_user_action = "ignore")
  y2 <- mod2$run(seq.int(inits$t, 10))
  y2 <- mod2$transform_variables(y2)
  
  expect_equivalent(y1$U[y1$t >= 5, , , drop = FALSE], y2$U, tol = 0.1)
  expect_equivalent(y1$lambda[y1$t >= 5, , drop = FALSE], y2$lambda, tol = 1e-5)
  
})



test_that("t_stop is working correctly", {                                                    #come back to this later 
  ## check with single parameter set
  vp <- vax_params_xvwv(vbe = 0, uptake = 1, strategy = "VoD", vea = 1,
                        t_stop = 2 / 365)
  params <- model_params(gono_params = gono_params(1)[[1]],
                         vax_params = vp)
  mod <- model(user = params, unused_user_action = "ignore")
  tt <- seq.int(0, 5) / 365
  y <- mod$run(tt, )
  y <- mod$transform_variables(y)
  
  expect_equal(diff(apply(y$cum_vaccinated, 1, sum))[tt[-length(tt)] > 2 / 365],
               c(0, 0))
  
})



test_that("aggregated time series output correctly", {                                        #come back to this later - do we need it ? 
  ## check with single parameter set
  params <- model_params(gono_params = gono_params(1)[[1]])
  mod <- model(user = params, unused_user_action = "ignore")
  tt <- seq.int(0, 5) / 365
  y <- mod$run(tt)
  y <- mod$transform_variables(y)
  expect_equal(y$tot_treated, apply(y$cum_treated, 1, sum))
  expect_equal(y$tot_attended, apply(y$cum_screened, 1, sum) + y$tot_treated)
  
})


### lets test if this runs lol
#devtools::test_active_file()

