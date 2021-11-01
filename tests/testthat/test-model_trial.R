#Call in model_trial

######### testing ########
context("trial model check")

test_that("no new infections if lambda is 0", {
  params <- gonovax::model_params(gono_params = gonovax::gono_params(1)[[1]])                #gono_params() and model_params() generate params for testing
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
  params <- model_params(gono_params = gono_params(1)[[1]]) 
  
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

### lets test if this runs lol
devtools::test_active_file()
  