#Call in model_trial
source("C:/Users/Dariya/OneDrive - Imperial College London/Documents/gonovax/inst/odin/model_trial.R")

#want to create an instance BUT need to pass it our model parameters
params <- gonovax::model_params(gonovax::gono_params(1)[[1]])     #gono_params is function that generates params for testing #model 
                                                                  #params sets the other model parameters
#change state variables to 1*1 matrix rather than 2*1 , sum across so all go into 1 group
states <- c(18,19,20,21,22)
for (i in states) {
params[[i]]<- matrix(params[[i]][1,1] + params[[i]][2,1], nrow =1, ncol = 1)
print(params[i])
}
params$vbe <- array(0, dim = c(1, 1, 1))
params$vos <- array(0, dim = c(1, 1, 1))
params$vod <- array(0, dim = c(1, 1, 1))

params$lambda <- 2                                                   #added lambda because we made it user defined

mod <- model_trial$new(user = params, unused_user_action = "ignore")     #hand over user defined parameters to the model
tt <- seq.int(0, 5) / 365                                            #define the time steps
y <- mod$run(t = tt)                                                #run model
y <- mod$transform_variables(y)


######### testing #########
context("trial model check")

test_that("no new infections if lambda is 0", {
    params <- gonovax::model_params(gonovax::gono_params(1)[[1]]) 
    states <- c(18,19,20,21,22)
    for (i in states) {
      params[[i]]<- matrix(params[[i]][1,1] + params[[i]][2,1], nrow =1, ncol = 1)
      print(params[i])
    }
    
    params$vbe <- array(0, dim = c(1, 1, 1))
    params$vos <- array(0, dim = c(1, 1, 1))
    params$vod <- array(0, dim = c(1, 1, 1))
    
    params$lambda <- 0  
    
    mod <- model_trial$new(user = params, unused_user_action = "ignore")     
    tt <- seq.int(0, 5) / 365                                             
    y <- mod$run(t = tt)                                                 
    y <- mod$transform_variables(y)
    
  
    expect_true(all(y$I ==0))             #incubating remain at 0
    expect_true(all(y$S ==0))             #no symptomatic infections
    expect_true(all(cum_incid == 0))
})