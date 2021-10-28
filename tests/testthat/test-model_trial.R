#Call in model_trial

######### testing ########
context("trial model check")

test_that("no new infections if lambda is 0", {
    params <- gonovax::model_params(gonovax::gono_params(1)[[1]])                  #gono_params() and model_params() generate params for testing
    states <- c(18,19,20,21,22)                                                    #sum across activity groups and make 1 i
    for (i in states) {    
      params[[i]]<- matrix(params[[i]][1,1] + params[[i]][2,1], nrow =1, ncol = 1) 
      print(params[i])
    }
    
    params$vbe <- array(0, dim = c(1, 1, 1))                                       #adjust dimensions of vaccination strategies 
    params$vos <- array(0, dim = c(1, 1, 1))
    params$vod <- array(0, dim = c(1, 1, 1))
    
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
  params <- gonovax::model_params(gonovax::gono_params(1)[[1]])                  #gono_params() and model_params() generate params for testing
  states <- c(18,19,20,21,22)                                                    #sum across activity groups and make 1 i
  for (i in states) {    
    params[[i]]<- matrix(params[[i]][1,1] + params[[i]][2,1], nrow =1, ncol = 1) 
    print(params[i])
  }
  
  params$vbe <- array(0, dim = c(1, 1, 1))                                       #adjust dimensions of vaccination strategies 
  params$vos <- array(0, dim = c(1, 1, 1))
  params$vod <- array(0, dim = c(1, 1, 1))
  
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
  