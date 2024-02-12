context("model (check)")

test_that("there are no infections when beta is 0", {
  
for (i in (1:10)){
  
  params <- model_params(gono_params = gono_params(1)[[1]], n_diag_rec = i)
  params$beta_t[] <- 0
  
  mod <- model$new(user = params, unused_user_action = "ignore")
  
  tt <- seq.int(0, 5) / 365
  y <- mod$run(tt)
  y <- mod$transform_variables(y)
  expect_true(all(y$I == 0))
  expect_true(all(y$S == 0))
  expect_true(all(y$cum_incid == 0))
  expect_true(all(unlist(y) >= 0))
  
  
}
  
})


test_that("there are no symptomatic infections when psi = 0", {

  
  for (i in (1:5)){  
  
  params <- model_params(gono_params = gono_params(1)[[1]], n_diag_rec = i)

  
  params$psi <- 0
  mod <- model$new(user = params, unused_user_action = "ignore")
  tt <- seq.int(0, 5) / 365 
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  expect_true(any(y$I == 0))
  expect_true(all(y$S == 0))
  expect_true(all(y$cum_diag_s == 0))
  expect_true(all(y$cum_diag_a[-1, , ] > 0))
  expect_true(all(unlist(y) >= 0))
  
  
  }
  
  
  
})

test_that("there are no asymptomatic infections when psi = 1", {

  
for (i in 1:5){  
  
  
  params <- model_params(gono_params = gono_params(1)[[1]], n_diag_rec = i)

  
  params$psi <- 1
  params$S0[, ] <- params$A0[, ]
  params$A0[, ] <- 0
  mod <- model$new(user = params, unused_user_action = "ignore")
  tt <- seq.int(0, 5)/365 
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  expect_true(any(y$I == 0))
  expect_true(all(y$A == 0))
  expect_true(all(y$cum_diag_a == 0))
  expect_true(all(y$cum_diag_s[-1, , ] > 0))
  expect_true(all(unlist(y) >= 0))
  
}
  
  
})

test_that("there are no infections when A0 = 0", {

for (i in 1:5){  
  
  params <- model_params(gono_params = gono_params(1)[[1]], n_diag_rec = i)
  params$A0[, ] <- 0
  mod <- model$new(user = params, unused_user_action = "ignore")
  tt <- seq.int(0, 5) / 365
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  
  expect_true(all(y$I == 0))
  expect_true(all(y$A == 0))
  expect_true(all(y$S == 0))
  expect_true(all(unlist(y) >= 0))
  
}
  
})


test_that("no-one is treated when mu and eta = 0", {

  for (i in 1:5){
  
  params <- model_params(gono_params = gono_params(1)[[1]], n_diag_rec = i)
  params$mu <- params$eta_h_t[] <- params$eta_l_t[] <-  0
  mod <- model$new(user = params, unused_user_action = "ignore")
  
  tt <- seq.int(0, 5) / 365
  y <- mod$run(tt)
  y <- mod$transform_variables(y)
  expect_true(all(y$T == 0))
  expect_true(all(y$cum_treated == 0))
  expect_true(all(unlist(y) >= 0))
  
  }
  
  
})


test_that("Diagnosis waning is working correctly",{
  
  #Set A0 to 0 and move all uninfecteds to diagnosed once stratum
  
  # set birth and death rates to 0 as well
  
  params <- model_params(gono_params = gono_params(1)[[1]], n_diag_rec = 2)

  params$enr = 0
  params$exr = 0
  params$A0[, ] <- 0
  params$U0[,2] = params$U0[,1]
  params$U0[,1] = 0
  

  mod <- model$new(user = params, unused_user_action = "ignore")
  tt <- seq.int(0,2)
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  
  expect_equal(y$U[1,1,1],  y$U[1,1,2]*(1-exp(-0)))
  expect_equal(y$U[1,2,1],  y$U[1,2,2]*(1-exp(-0)))
  expect_equal(y$U[2,1,1],  y$U[1,1,2]*(1-exp(-1)), tolerance = 1e-6)
  expect_equal(y$U[2,2,1],  y$U[1,2,2]*(1-exp(-1)), tolerance = 1e-6)
  expect_equal(y$U[3,1,1],  y$U[1,1,2]*(1-exp(-2)), tolerance = 1e-6)
  expect_equal(y$U[3,2,1],  y$U[1,2,2]*(1-exp(-2)), tolerance = 1e-6)
  
  
  
  for (i in 3:5){
  
  params <- model_params(gono_params = gono_params(1)[[1]], n_diag_rec = i)
  
  params$enr = 0
  params$exr = 0
  params$A0[, ] <- 0
  params$U0[,i] = params$U0[,1]
  params$U0[,1] = 0
  
  
  mod <- model$new(user = params, unused_user_action = "ignore")
  tt <- seq.int(0,2)
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  
  expect_equal(sum(y$U[1,1,1:(i-1)]),  y$U[1,1,i]*(1-exp(-0)))
  expect_equal(sum(y$U[1,2,1:(i-1)]),  y$U[1,2,i]*(1-exp(-0)))
  expect_equal(sum(y$U[2,1,1:(i-1)]),  y$U[1,1,i]*(1-exp(-1)), tolerance = 1e-6)
  expect_equal(sum(y$U[2,2,1:(i-1)]),  y$U[1,2,i]*(1-exp(-1)), tolerance = 1e-6)
  expect_equal(sum(y$U[3,1,1:(i-1)]),  y$U[1,1,i]*(1-exp(-2)), tolerance = 1e-6)
  expect_equal(sum(y$U[3,2, 1:(i-1)]),  y$U[1,2,i]*(1-exp(-2)), tolerance = 1e-6)
  

  
  for (ii in (0:(i-1))){
    
    if (ii < (i-1)) {
      expect_equal(y$U[2,1, (i -ii)], y$U[1,1,i]*dpois(ii, 1))
      expect_equal(y$U[2,2, (i -ii)], y$U[1,2,i]*dpois(ii, 1))
      
    }
    else if (ii == i-1){
      expect_equal(y$U[2,1, (i -ii)], y$U[1,1,i]*(1 - ppois(ii-1, 1)))
      expect_equal(y$U[2,2, (i -ii)], y$U[1,2,i]*(1 - ppois(ii-1, 1)))
      
    } 
    
  }
  
  }

})






test_that("the foi is calculated correctly", {
  vei <- 0.123
  vax_params <- vax_params_xvwv(uptake = 0.5, dur = 1,
                                strategy = "VoA", vei = 0.123, n_diag_rec = 1)
  
  
  for (i in 1:5){
  
  vax_params <- vax_params_xvwv(uptake = 0.5, dur = 1,
                                          strategy = "VoA", vei = 0.123, n_diag_rec = i)
  params<- model_params(gono_params = gono_params(1)[[1]],
                                  vax_params = vax_params, n_diag_rec =  i)
  
  
  expect_true(length(params$beta_t) > 0)
  
  
  mod <- model$new(user = params, unused_user_action = "ignore")
  
  tt <- seq.int(0, 5) / 365
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  
  pL <- params$p[1]
  pH <- params$p[2]
  NL <- rowSums(y$N[, 1, ])
  NH <- rowSums(y$N[, 2, ])

  C <- y$I + y$A + y$S
  CL <- c(C[, 1, ] %*% (1 - vax_params$vei))
  CH <- c(C[, 2, ] %*% (1 - vax_params$vei))
  
  #  CL <- C[, 1, ]
  #  CH <- C[,2, ]
  eps <- params$epsilon
  beta <- params$beta_t
  
  np <- pL * NL + pH * NH
  npL <- pL * NL / np
  npH <- pH * NH / np
  
  # calculate FOI
  foi_cross <- (1 - eps) * (npL * CL / NL + npH * CH / NH)
  foi_L <- pL * beta * (eps * CL / NL + foi_cross)
  foi_H <- pH * beta * (eps * CH / NH + foi_cross)
  
  expect_equal(y$lambda[, 1], foi_L)
  expect_equal(y$lambda[, 2], foi_H)
  
  
  }
  
  
})


test_that("Bex model runs with no vaccination", {
  tt <- seq.int(0, 5) / 365
  
  for (i in 1:5){
    
  n_diag_rec <- i  
  
  params0 <- model_params(gono_params = gono_params(1)[[1]], n_diag_rec = n_diag_rec)
  mod0 <- model$new(user = params0, unused_user_action = "ignore")
  y0 <- mod0$run(t = tt)
  y0 <- mod0$transform_variables(y0)
  
  params1 <- model_params(gono_params = gono_params(1)[[1]],
                          vax_params = vax_params_xvwv(vbe = 0, n_diag_rec = n_diag_rec), n_diag_rec = n_diag_rec)
  mod1 <- model$new(user = params1, unused_user_action = "ignore")
  y1 <- mod1$run(t = tt)
  y1 <- mod1$transform_variables(y1)
  
  # check that nil vaccination gives same results as before
  expect_true(all(y1$U[, , 1:n_diag_rec, drop = FALSE] == y0$U[, , 1:n_diag_rec, drop = FALSE]))
  expect_true(all(y1$A[, , 1:n_diag_rec, drop = FALSE] == y0$A[, , 1:n_diag_rec, drop = FALSE]))
  expect_true(all(y1$S[, , 1:n_diag_rec, drop = FALSE] == y0$S[, , 1:n_diag_rec, drop = FALSE]))
  expect_true(all(y1$T[, , 1:n_diag_rec, drop = FALSE] == y0$T[, , 1:n_diag_rec, drop = FALSE]))
  
  expect_true(all(y1$N[, , (n_diag_rec+1):(2*n_diag_rec)] == 0))
  expect_true(all(apply(y1$N, c(1, 2), sum) - 6e5 < 1e-6))
  
  }
  
})

test_that("Bex model runs with vbe", {
  tt <- seq.int(0, 2) / 365

  
  
  for (i in 1:5){ 

  n_diag_rec <- i
  # with perfect efficacy
  params <- model_params(gono_params = gono_params(1)[[1]],
                         vax_params = vax_params_xvwv(vbe = 1, vea = 1, n_diag_rec = n_diag_rec), n_diag_rec = n_diag_rec)
  
  #params$beta_t[] <- 0
  #params$enr <- 0
  #params$exr <- 0
  #params$A0[, ] <- 0
  
  
  mod <- model$new(user = params, unused_user_action = "ignore")
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)

  
  
  temp = array(rep(0, 3*2*3), dim = c(3,2,3))
  
  if (n_diag_rec > 1){
    temp[,,1] = rowSums(y$U[,,1:n_diag_rec],dims = 2)
    temp[,,2] = rowSums(y$U[,,(n_diag_rec+1):(2*n_diag_rec)],dims = 2)
    temp[,,3] = rowSums(y$U[,,(2*n_diag_rec+1):(3*n_diag_rec)],dims = 2)
    
  }else{
    temp[,,1] = y$U[,,1]
    temp[,,2] = y$U[,,2]
    temp[,,3] = y$U[,,3]
  }

  
  expect_equal(temp,
               array(c(505266, 505270.782413465, 505276.010170825,
                       89303, 89292.2121354773, 89281.4753785608,
                       0, 27.9444015916977, 55.8871954695537,
                       0, 4.93136498677018, 9.862446259333,
                       0, 3.82796084566103e-05, 0.000153112453356132,
                       0, 6.75495729408945e-06, 2.70176968895868e-05),
                     dim = c(3, 2, 3)))
  # check some people are being vaccinated
  expect_true(all(y$U[-1, , n_diag_rec+1] > 0))
  # check no compartments are leaking
  expect_true(all(apply(y$N, c(1, 2), sum) - 6e5 < 1e-6))
  # check all entrants are vaccinated
  expect_equal(y$cum_offered_vbe, y$cum_vbe)
  # check there are infections in unvaccinated group
  expect_false(all(y$I[, , 1:n_diag_rec] == 0))
  expect_false(all(y$A[, , 1:n_diag_rec] == 0))
  expect_false(all(y$S[, , 1:n_diag_rec] == 0))
  expect_false(all(y$T[, , 1:n_diag_rec] == 0))
  # check there are no infections in vaccinated group
  expect_true(all(y$I[, , (n_diag_rec+1):(2*n_diag_rec)] == 0))
  expect_true(all(y$A[, , (n_diag_rec+1):(2*n_diag_rec)] == 0))
  expect_true(all(y$S[, , (n_diag_rec+1):(2*n_diag_rec)] == 0))
  expect_true(all(y$T[, , (n_diag_rec+1):(2*n_diag_rec)] == 0))
  
  
  
  
}
  
  
  
  
})


test_that("Check vaccination on screening in Bex model", {
  tt <- seq.int(0, 2) / 365
  
  
for (i in 1:5){
  
  n_diag_rec <- i
  
  params <-
    model_params(gono_params = gono_params(1)[[1]],
                           vax_params = vax_params_xvwv(vbe = 0, uptake = 1,
                                                                  strategy = "VoS", vea = 1, n_diag_rec = n_diag_rec), n_diag_rec = n_diag_rec)
  
  
  mod <- model$new(user = params, unused_user_action = "ignore")
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  
  
  temp = array(rep(0, 3*2*3), dim = c(3,2,3))
  
  if (n_diag_rec > 1){
    temp[,,1] = rowSums(y$U[,,1:n_diag_rec],dims = 2)
    temp[,,2] = rowSums(y$U[,,(n_diag_rec+1):(2*n_diag_rec)],dims = 2)
    temp[,,3] = rowSums(y$U[,,(2*n_diag_rec+1):(3*n_diag_rec)],dims = 2)
    
  }else{
    temp[,,1] = y$U[,,1]
    temp[,,2] = y$U[,,2]
    temp[,,3] = y$U[,,3]
  }
  
  expect_equal(temp,
               array(c(505266, 504689.65898183, 504114.492415873,
                       89303, 89189.5072358554, 89076.2208839755,
                       0, 609.068445893679, 1217.40743007721,
                       0, 107.642503971552, 215.142050698242,
                       0, 0.000834155028318817, 0.0033338806826574,
                       0, 0.000147420094481985, 0.000589147144744636),
                     dim = c(3L, 2L, 3L)))
  
  
  # check some people are being vaccinated
  expect_true(all(y$U[-1, , (n_diag_rec+1):(2*n_diag_rec)] > 0))
  expect_true(all(y$cum_vaccinated[-1, , 1] > 0))
  expect_true(all(y$cum_vaccinated[, , (n_diag_rec+1):(2*n_diag_rec)] == 0))
  # check all those treated were vaccinated
  expect_true(all(y$cum_vaccinated[, , 1:n_diag_rec] == y$cum_screened[, , 1:n_diag_rec]))
  expect_true(all(y$cum_vaccinated[, , (2*n_diag_rec+1):(3*n_diag_rec)] == y$cum_screened[, , (2*n_diag_rec+1):(3*n_diag_rec)]))
  # check no compartments are leaking
  expect_true(all(apply(y$N, 1, sum) - 6e5 < 1e-6))
  # check there are infections in unvaccinated group
  expect_false(all(y$I[, , 1:n_diag_rec] == 0))
  expect_false(all(y$A[, , 1:n_diag_rec] == 0))
  expect_false(all(y$S[, , 1:n_diag_rec] == 0))
  expect_false(all(y$T[, , 1:n_diag_rec] == 0))
  expect_false(all(y$I[, , (2*n_diag_rec+1):(3*n_diag_rec)] == 0))
  expect_false(all(y$A[, , (2*n_diag_rec+1):(3*n_diag_rec)] == 0))
  expect_false(all(y$S[, , (2*n_diag_rec+1):(3*n_diag_rec)] == 0))
  expect_false(all(y$T[, , (2*n_diag_rec+1):(3*n_diag_rec)] == 0))
  # check there are no infections in vaccinated group
  expect_true(all(y$I[, , (1*n_diag_rec+1):(2*n_diag_rec)] == 0))
  expect_true(all(y$A[, , (1*n_diag_rec+1):(2*n_diag_rec)] == 0))
  expect_true(all(y$S[, , (1*n_diag_rec+1):(2*n_diag_rec)] == 0))
  expect_true(all(y$T[, , (1*n_diag_rec+1):(2*n_diag_rec)] == 0))
  
  
}
  
  
})

test_that("Check vaccination on diagnosis in Bex model", {
  tt <- seq.int(0, 2) / 365
  
  
  for (i in 1:5){
  
  n_diag_rec <- i
  
  # with perfect efficacy
  params <-
    model_params(gono_params = gono_params(1)[[1]],
                           vax_params = vax_params_xvwv(vbe = 0, uptake = 1,
                                                                  strategy = "VoD", vea = 1, n_diag_rec = n_diag_rec), n_diag_rec = n_diag_rec)
  mod <- model$new(user = params, unused_user_action = "ignore")
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  
  temp = array(rep(0, 3*2*3), dim = c(3,2,3))
  
  if (n_diag_rec > 1){
    temp[,,1] = rowSums(y$U[,,1:n_diag_rec],dims = 2)
    temp[,,2] = rowSums(y$U[,,(n_diag_rec+1):(2*n_diag_rec)],dims = 2)
    temp[,,3] = rowSums(y$U[,,(2*n_diag_rec+1):(3*n_diag_rec)],dims = 2)
    
  }else{
    temp[,,1] = y$U[,,1]
    temp[,,2] = y$U[,,2]
    temp[,,3] = y$U[,,3]
  }
  
  expect_equal(temp,
               array(c(505266, 505298.314652489, 505330.328904702,
                       89303, 89297.0812095395, 89291.0874630527,
                       0, 0.412133420953426, 1.56834473390679,
                       0, 0.06199394863926, 0.249176144196752,
                       0, 3.81264799241714e-07, 2.93702682328813e-06,
                       0, 5.68712159308333e-08, 4.53390060409456e-07),
                     dim = c(3L, 2L, 3L)))
  
  
  # check some people are being vaccinated
  expect_true(all(y$U[-1, , (n_diag_rec+1):(2*n_diag_rec)] > 0))
  expect_true(all(y$cum_vaccinated[-1, , 1] > 0))
  expect_true(all(y$cum_vaccinated[, , (n_diag_rec+1):(2*n_diag_rec)] == 0))
  # check all those treated were vaccinated
  expect_true(all(y$cum_vaccinated == y$cum_treated))
  # check no compartments are leaking
  expect_true(all(apply(y$N, 1, sum) - 6e5 < 1e-6))
  # check there are infections in unvaccinated group
  expect_false(all(y$I[, , 1:n_diag_rec] == 0))
  expect_false(all(y$A[, , 1:n_diag_rec] == 0))
  expect_false(all(y$S[, , 1:n_diag_rec] == 0))
  expect_false(all(y$T[, , 1:n_diag_rec] == 0))
  expect_false(all(y$I[, , (2*n_diag_rec+1):(3*n_diag_rec)] == 0))
  expect_false(all(y$A[, , (2*n_diag_rec+1):(3*n_diag_rec)] == 0))
  expect_false(all(y$S[, , (2*n_diag_rec+1):(3*n_diag_rec)] == 0))
  expect_false(all(y$T[, , (2*n_diag_rec+1):(3*n_diag_rec)] == 0))
  # check there are no infections in vaccinated group
  expect_true(all(y$I[, , (1*n_diag_rec+1):(2*n_diag_rec)] == 0))
  expect_true(all(y$A[, , (1*n_diag_rec+1):(2*n_diag_rec)] == 0))
  expect_true(all(y$S[, , (1*n_diag_rec+1):(2*n_diag_rec)] == 0))
  expect_true(all(y$T[, , (1*n_diag_rec+1):(2*n_diag_rec)] == 0))
  
  }
  
  
})

test_that("Vaccination according to diagnosis history works as expected with the Bex model",{
  
  tt <- seq.int(0, 2) / 365
  
  
  n_diag_rec <- 1
  
  
  # with perfect efficacy
  params <-
    model_params(gono_params = gono_params(1)[[1]],
                           vax_params = vax_params_xvwv(vbe = 0, uptake = 1,
                                                                  strategy = "VaH", vea = 1, n_diag_rec = n_diag_rec), n_diag_rec = n_diag_rec)
  
  #set vod to 0 for tests!
  params$vod[,,] = 0
  
  mod <- model$new(user = params, unused_user_action = "ignore")
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)

  # no diagnosis history so no vaccination!  
  expect_true(all(y$N[, , (n_diag_rec+1):(2*n_diag_rec)] == 0))
  expect_true(all(apply(y$N, c(1, 2), sum) - 6e5 < 1e-6))
  
   y
  
  
  n_diag_rec <- 3
  
  params <-
    model_params(gono_params = gono_params(1)[[1]],
                           vax_params = vax_params_xvwv(vbe = 0, uptake = 1,
                                                                  strategy = "VaH", vea = 1, n_diag_rec = n_diag_rec), n_diag_rec = n_diag_rec)
  
  #set vod to 0 for tests!
  params$vod[,,] = 0
  
  mod <- model$new(user = params, unused_user_action = "ignore")
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  

  #vaccinations
  expect_false(all(y$N[-1, , (n_diag_rec+1):(2*n_diag_rec)] == 0))

  expect_true(all(apply(y$N, c(1, 2), sum) - 6e5 < 1e-6))
  
  ## vaccines in all recent diagnosis strata
  expect_true(all(y$cum_vaccinated[-1, , (2:n_diag_rec)] > 0))
  expect_true(all(y$cum_vaccinated[-1, , ((2*n_diag_rec)+2):(3*n_diag_rec)] > 0))
  
  # vaccinations the sum of screenings among those recently diagnosed
  expect_true(all(y$cum_vaccinated[, , 2:n_diag_rec] == y$cum_screened[, , 2:n_diag_rec]))
  expect_true(all(y$cum_vaccinated[, , (2*n_diag_rec+2):(3*n_diag_rec)] == y$cum_screened[, , (2*n_diag_rec+2):(3*n_diag_rec)]))
  

  
  # For a vaccine with 0% efficacy, check diagnoses are the same between model
  # without diagnosis history and model with diagnosis history

  
  tt <- seq.int(0, 50) 

  n_diag_rec <- 1

  # with 0% efficacy
  params <-
    model_params(gono_params = gono_params(1)[[1]],
                           vax_params = vax_params_xvwv(vbe = 0, uptake = 1,
                                                                  strategy = "VaH", vea = 0, n_diag_rec = n_diag_rec), n_diag_rec = n_diag_rec)
  #set vod to 0 for tests!
  params$vod[,,] = 0
  mod <- model$new(user = params, unused_user_action = "ignore")
  y0 <- mod$run(t = tt)
  y0 <- mod$transform_variables(y0)
  
  
  n_diag_rec <- 3
  # with 0% efficacy
  params <-
    model_params(gono_params = gono_params(1)[[1]],
                           vax_params = vax_params_xvwv(vbe = 0, uptake = 1,
                                                                  strategy = "VaH", vea = 0, n_diag_rec = n_diag_rec), n_diag_rec = n_diag_rec)
  #set vod to 0 for tests!
  params$vod[,,] = 0
  mod <- model$new(user = params, unused_user_action = "ignore")
  y3 <- mod$run(t = tt)
  y3 <- mod$transform_variables(y3)
  
  #check diagnoses are the same between models
  expect_equal(rowSums(y0$cum_diag_a[2:51,,]), rowSums(y3$cum_diag_a[2:51,,]), tolerance = 1e-6)
  expect_equal(rowSums(y0$cum_diag_s[2:51,,]), rowSums(y3$cum_diag_s[2:51,,]), tolerance = 1e-6)
  
  
})




test_that("XPVWRH model runs with no vaccination", {
  tt <- seq.int(0, 5) / 365

  
  for (i in 1:5){
    
  n_diag_rec <- i
  
  params0 <- model_params(gono_params = gono_params(1)[[1]], n_diag_rec = n_diag_rec)
  mod0 <- model$new(user = params0, unused_user_action = "ignore")
  y0 <- mod0$run(t = tt)
  y0 <- mod0$transform_variables(y0)
  
  params1 <- model_params_xpvwrh(gono_params = gono_params(1)[[1]],
                                    vax_params = vax_params_xpvwrh(vbe = 0, n_erlang = 1, n_diag_rec = n_diag_rec), n_diag_rec = n_diag_rec)
  mod1 <- model$new(user = params1, unused_user_action = "ignore")
  y1 <- mod1$run(t = tt)
  y1 <- mod1$transform_variables(y1)
  
  # check that nil vaccination gives same results as before
  expect_true(all(y1$U[, , 1:n_diag_rec, drop = FALSE] == y0$U[, , 1:n_diag_rec, drop = FALSE]))
  expect_true(all(y1$A[, , 1:n_diag_rec, drop = FALSE] == y0$A[, , 1:n_diag_rec, drop = FALSE]))
  expect_true(all(y1$S[, , 1:n_diag_rec, drop = FALSE] == y0$S[, , 1:n_diag_rec, drop = FALSE]))
  expect_true(all(y1$T[, , 1:n_diag_rec, drop = FALSE] == y0$T[, , 1:n_diag_rec, drop = FALSE]))
  
  expect_true(all(y1$N[, , (n_diag_rec+1):(2*n_diag_rec)] == 0))
  expect_true(all(apply(y1$N, c(1, 2), sum) - 6e5 < 1e-6))
  
  }
  
})


test_that("XPVWRH model runs with vbe", {
  tt <- seq.int(0, 2) / 365


  
  for (i in 1:5){
  
  
  
  n_diag_rec <-i
  # with perfect efficacy
  params <- model_params_xpvwrh(gono_params = gono_params(1)[[1]],
                                   vax_params = vax_params_xpvwrh(vbe = 1, vea = 1, n_erlang = 1, n_diag_rec = n_diag_rec), n_diag_rec = i)

  
  mod <- model$new(user = params, unused_user_action = "ignore")
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  
  
  temp = array(rep(0, 3*2*3), dim = c(3,2,3))
  
  if (n_diag_rec > 1){
    temp[,,1] = rowSums(y$U[,,1:n_diag_rec],dims = 2)
    temp[,,2] = rowSums(y$U[,,(2*n_diag_rec+1):(3*n_diag_rec)],dims = 2)
    temp[,,3] = rowSums(y$U[,,(4*n_diag_rec+1):(5*n_diag_rec)],dims = 2)
    
  }else{
    temp[,,1] = y$U[,,1]
    temp[,,2] = y$U[,,3]
    temp[,,3] = y$U[,,5]
  }
  
  
  expect_equal(temp,
               array(c(505266, 505270.782413465, 505276.010170825,
                       89303, 89292.2121354773, 89281.4753785608,
                       0, 27.9444015916977, 55.8871954695537,
                       0, 4.93136498677018, 9.862446259333,
                       0, 3.82796084566103e-05, 0.000153112453356132,
                       0, 6.75495729408945e-06, 2.70176968895868e-05),
                     dim = c(3, 2, 3)))
  # check some people are being vaccinated
  expect_true(all(y$U[-1, , 2*n_diag_rec+1] > 0))
  # check no compartments are leaking
  expect_true(all(apply(y$N, c(1, 2), sum) - 6e5 < 1e-6))
  # check all entrants are vaccinated
  expect_equal(y$cum_offered_vbe, y$cum_vbe)
  # check there are infections in unvaccinated group
  expect_false(all(y$I[, , 1:n_diag_rec] == 0))
  expect_false(all(y$A[, , 1:n_diag_rec] == 0))
  expect_false(all(y$S[, , 1:n_diag_rec] == 0))
  expect_false(all(y$T[, , 1:n_diag_rec] == 0))
  # check there are no infections in vaccinated group
  expect_true(all(y$I[, , (2*n_diag_rec+1):(3*n_diag_rec)] == 0))
  expect_true(all(y$A[, , (2*n_diag_rec+1):(3*n_diag_rec)] == 0))
  expect_true(all(y$S[, , (2*n_diag_rec+1):(3*n_diag_rec)] == 0))
  expect_true(all(y$T[, , (2*n_diag_rec+1):(3*n_diag_rec)] == 0))
  
  
  }
  
  
  
  
})


test_that("Check vaccination on screening in XPVWRH model", {
  tt <- seq.int(0, 2) / 365
  
  
  for (i in 1:5){
  
  n_diag_rec <- i
  
  params <-
    model_params_xpvwrh(gono_params = gono_params(1)[[1]],
                           vax_params = vax_params_xpvwrh(vbe = 0, r1 = 1, r2 = 1,
                                                                  strategy = "VoS", vea = 1, n_erlang = 1, n_diag_rec = n_diag_rec), n_diag_rec = n_diag_rec)
  
  
  mod <- model$new(user = params, unused_user_action = "ignore")
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  
  
  temp = array(rep(0, 3*2*3), dim = c(3,2,3))
  
  if (n_diag_rec > 1){
    temp[,,1] = rowSums(y$U[,,1:n_diag_rec],dims = 2)
    temp[,,2] = rowSums(y$U[,,(2*n_diag_rec+1):(3*n_diag_rec)],dims = 2)
    temp[,,3] = rowSums(y$U[,,(4*n_diag_rec+1):(5*n_diag_rec)],dims = 2)
    
  }else{
    temp[,,1] = y$U[,,1]
    temp[,,2] = y$U[,,3]
    temp[,,3] = y$U[,,5]
  }
  
  
  expect_equal(temp,
               array(c(505266, 504689.65898183, 504114.492415873,
                       89303, 89189.5072358554, 89076.2208839755,
                       0, 609.068445893679, 1217.40743007721,
                       0, 107.642503971552, 215.142050698242,
                       0, 0.000834155028318817, 0.0033338806826574,
                       0, 0.000147420094481985, 0.000589147144744636),
                     dim = c(3L, 2L, 3L)))
  
  
  # check some people are being vaccinated
  expect_true(all(y$U[-1, , (2*n_diag_rec+1):(3*n_diag_rec)] > 0))
  expect_true(all(y$cum_vaccinated[-1, , 1] > 0))
  expect_true(all(y$cum_vaccinated[, , (2*n_diag_rec+1):(3*n_diag_rec)] == 0))
  # check all those treated were vaccinated
  expect_true(all(y$cum_vaccinated[, , 1:n_diag_rec] == y$cum_screened[, , 1:n_diag_rec]))
  expect_true(all(y$cum_vaccinated[, , (3*n_diag_rec+1):(4*n_diag_rec)] == y$cum_screened[, , (3*n_diag_rec+1):(4*n_diag_rec)]))
  # check no compartments are leaking
  expect_true(all(apply(y$N, 1, sum) - 6e5 < 1e-6))
  # check there are infections in unvaccinated group
  expect_false(all(y$I[, , 1:n_diag_rec] == 0))
  expect_false(all(y$A[, , 1:n_diag_rec] == 0))
  expect_false(all(y$S[, , 1:n_diag_rec] == 0))
  expect_false(all(y$T[, , 1:n_diag_rec] == 0))
  expect_false(all(y$I[, , (3*n_diag_rec+1):(4*n_diag_rec)] == 0))
  expect_false(all(y$A[, , (3*n_diag_rec+1):(4*n_diag_rec)] == 0))
  expect_false(all(y$S[, , (3*n_diag_rec+1):(4*n_diag_rec)] == 0))
  expect_false(all(y$T[, , (3*n_diag_rec+1):(4*n_diag_rec)] == 0))
  # check there are no infections in vaccinated group
  expect_true(all(y$I[, , (2*n_diag_rec+1):(3*n_diag_rec)] == 0))
  expect_true(all(y$A[, , (2*n_diag_rec+1):(3*n_diag_rec)] == 0))
  expect_true(all(y$S[, , (2*n_diag_rec+1):(3*n_diag_rec)] == 0))
  expect_true(all(y$T[, , (2*n_diag_rec+1):(3*n_diag_rec)] == 0))
  
  }
  
  
  
})


test_that("Check vaccination on diagnosis in XPVWRH model", {
  tt <- seq.int(0, 2) / 365
  
  for (i in 1:5){  

  n_diag_rec <- i
  
  # with perfect efficacy
  params <-
    model_params_xpvwrh(gono_params = gono_params(1)[[1]],
                           vax_params = vax_params_xpvwrh(vbe = 0, r1 = 1, r2 = 1,
                                                                  strategy = "VoD", vea = 1, n_erlang = 1, n_diag_rec = n_diag_rec), n_diag_rec = n_diag_rec)
  mod <- model$new(user = params, unused_user_action = "ignore")
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  
  temp = array(rep(0, 3*2*3), dim = c(3,2,3))

  
  
  if (n_diag_rec > 1){
    temp[,,1] = rowSums(y$U[,,1:n_diag_rec],dims = 2)
    temp[,,2] = rowSums(y$U[,,(2*n_diag_rec+1):(3*n_diag_rec)],dims = 2)
    temp[,,3] = rowSums(y$U[,,(4*n_diag_rec+1):(5*n_diag_rec)],dims = 2)
    
  }else{
    temp[,,1] = y$U[,,1]
    temp[,,2] = y$U[,,3]
    temp[,,3] = y$U[,,5]
  }
  
  expect_equal(temp,
               array(c(505266, 505298.314652489, 505330.328904702,
                       89303, 89297.0812095395, 89291.0874630527,
                       0, 0.412133420953426, 1.56834473390679,
                       0, 0.06199394863926, 0.249176144196752,
                       0, 3.81264799241714e-07, 2.93702682328813e-06,
                       0, 5.68712159308333e-08, 4.53390060409456e-07),
                     dim = c(3L, 2L, 3L)))
  
  
  # check some people are being vaccinated
  
  
  # no people would end up in twice diagnosed strata - because vaccine is 100% effective!
  expect_true(all(y$U[-1, , (2*n_diag_rec+1):(3*n_diag_rec-(n_diag_rec - 1))] > 0))
  
  
  
  expect_true(all(y$cum_vaccinated[-1, , 1] > 0))
  expect_true(all(y$cum_vaccinated[, , (2*n_diag_rec+1):(3*n_diag_rec)] == 0))
  # check all those treated were vaccinated
  expect_true(all(y$cum_vaccinated == y$cum_treated))
  # check no compartments are leaking
  expect_true(all(apply(y$N, 1, sum) - 6e5 < 1e-6))
  # check there are infections in unvaccinated group
  expect_false(all(y$I[, , 1:n_diag_rec] == 0))
  expect_false(all(y$A[, , 1:n_diag_rec] == 0))
  expect_false(all(y$S[, , 1:n_diag_rec] == 0))
  expect_false(all(y$T[, , 1:n_diag_rec] == 0))
  expect_false(all(y$I[, , (3*n_diag_rec+1):(4*n_diag_rec)] == 0))
  expect_false(all(y$A[, , (3*n_diag_rec+1):(4*n_diag_rec)] == 0))
  expect_false(all(y$S[, , (3*n_diag_rec+1):(4*n_diag_rec)] == 0))
  expect_false(all(y$T[, , (3*n_diag_rec+1):(4*n_diag_rec)] == 0))
  # check there are no infections in vaccinated group
  expect_true(all(y$I[, , (2*n_diag_rec+1):(3*n_diag_rec)] == 0))
  expect_true(all(y$A[, , (2*n_diag_rec+1):(3*n_diag_rec)] == 0))
  expect_true(all(y$S[, , (2*n_diag_rec+1):(3*n_diag_rec)] == 0))
  expect_true(all(y$T[, , (2*n_diag_rec+1):(3*n_diag_rec)] == 0))
  
  }
  
  
})


test_that("Vaccination according to history works as expected with the XPVWRH
          model",{
  
  tt <- seq.int(0, 2) / 365
  
  
  n_diag_rec <- 1
  
  
  # with imperfect efficacy
  params <-
    model_params_xpvwrh(gono_params = gono_params(1)[[1]],
                           vax_params = vax_params_xpvwrh(vbe = 0, r1 = 1, r2 = 1,
                                                                  strategy = "VaH", vea = 0.5, n_erlang = 1, n_diag_rec = n_diag_rec), n_diag_rec = n_diag_rec)
  
  #set vod to 0 for tests!
  params$vod[,,] = 0
  
  mod <- model$new(user = params, unused_user_action = "ignore")
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  
  # no diagnosis history so no vaccination!  
  expect_true(all(y$N[, , (2*n_diag_rec+1):(3*n_diag_rec)] == 0))
  expect_true(all(apply(y$N, c(1, 2), sum) - 6e5 < 1e-6))
  
  
  
  
  n_diag_rec <- 3
  

  
  ## Need to set some initial people up so that they can be recently diagnosed! or else they just end up vaccinated on diagnosis
  ## 17 July to pick up
  
  
  tt <- seq.int(0, 10)
  
  

  
  params0 <- model_params(gono_params = gono_params(1)[[1]], n_diag_rec = n_diag_rec)
  mod0 <- model$new(user = params0, unused_user_action = "ignore")
  y0 <- mod0$run(t = tt)
  y0 <- mod0$transform_variables(y0)
  
  n_par <- 1
  
  parameter_table <- read.csv(system.file("extdata/gono_params_t.csv",package='gonovax'))
  head(parameter_table)
  
  gono_params <- lapply(seq_len(n_par),
                        function(i) transform_fixed(parameter_table[i, ]))
  y0 <- run_onevax_xpvwrh(tt, gono_params, n_diag_rec = n_diag_rec)
  # get final (equilibrium) model state to use as starting point of run with vaccination
  init_params <- lapply(y0, restart_params)
  
  

  params <-
    model_params_xpvwrh(gono_params = gono_params[[1]], init_params= init_params[[1]],
                           vax_params = vax_params_xpvwrh(vbe = 0, r1 = 1, r2 = 1,
                                                                  strategy = "VaH", vea = 0, n_erlang = 1, n_diag_rec = n_diag_rec), n_diag_rec = n_diag_rec)
  
  #set vod to 0 for tests!
  
  
  #12 Dec
  ## not working because having VoRD scales down diag_rec!
  
  
  #set vod to 0 for tests!
  params$vod[,,] = 0

  tt <- seq.int(0, 2)/365
  
  
  mod <- model$new(user = params, unused_user_action = "ignore")
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  
  
  #vaccinations
  expect_false(all(y$N[-1, , (2*n_diag_rec+1):(3*n_diag_rec)] == 0))
  
  expect_true(all(apply(y$N, c(1, 2), sum) - 6e5 < 1e-6))
  
  ## vaccines in all unvaccinated recent diagnosis strata
  expect_true(all(y$cum_vaccinated[-1, , (2:n_diag_rec)] > 0))
  
  ## vaccines in all waned recent diagnosis strata
  expect_true(all(y$cum_vaccinated[-1, , ((3*n_diag_rec)+2):(4*n_diag_rec)] > 0))
  
  # vaccinations the sum of screenings among those recently diagnosed
  expect_true(all(y$cum_vaccinated_screen[, , 2:n_diag_rec] == y$cum_screened[, , 2:n_diag_rec]))
  expect_true(all(y$cum_vaccinated_screen[, , (3*n_diag_rec+2):(4*n_diag_rec)] == y$cum_screened[, , (3*n_diag_rec+2):(4*n_diag_rec)]))
  
  y
  
  

  # with 0% efficacy!
  
  tt <- seq.int(0, 50) 
  
  
  n_diag_rec <- 1
  
  
 vax_params = vax_params_xpvwrh(vbe = 0, r1 = 1, r2 = 1,
                                strategy = "VaH", vea = 0, n_diag_rec = n_diag_rec)
  
  
  params <-
    model_params_xpvwrh(gono_params = gono_params(1)[[1]],
                           vax_params = vax_params_xpvwrh(vbe = 0, r1 = 1, r2 = 1,
                                                                  strategy = "VaH", vea = 0, n_diag_rec = n_diag_rec), n_diag_rec = n_diag_rec)
  #set vod to 0 for tests!
  params$vod[,,] = 0
  mod <- model$new(user = params, unused_user_action = "ignore")
  y0 <- mod$run(t = tt)
  y0 <- mod$transform_variables(y0)
  
  
  n_diag_rec <- 3
  
  
  params <-
    model_params_xpvwrh(gono_params = gono_params(1)[[1]],
                           vax_params = vax_params_xpvwrh(vbe = 0, r1 = 1, r2 = 1,
                                                                  strategy = "VaH", vea = 0, n_diag_rec = n_diag_rec), n_diag_rec = n_diag_rec)
  #set vod to 0 for tests!
  params$vod[,,] = 0
  mod <- model$new(user = params, unused_user_action = "ignore")
  y3 <- mod$run(t = tt)
  y3 <- mod$transform_variables(y3)
  
  
  #check diagnoses are the same between models
  expect_equal(rowSums(y0$cum_diag_a[2:51,,]), rowSums(y3$cum_diag_a[2:51,,]), tolerance = 1e-6)
  expect_equal(rowSums(y0$cum_diag_s[2:51,,]), rowSums(y3$cum_diag_s[2:51,,]), tolerance = 1e-6)
  
  
})

test_that("can initialise after time 0", {
  
  for (i in 1:5){
    
  n_diag_rec <- i
  
  ## check with single parameter set
  params <- model_params(gono_params = gono_params(1)[[1]], n_diag_rec = n_diag_rec)
  mod <- model$new(user = params, unused_user_action = "ignore")
  
  tt <- seq.int(0, 5)
  y <- mod$run(tt)
  y <- mod$transform_variables(y)
  
  inits <- restart_params(y, n_vax = n_diag_rec)
  
  expect_true(all(y$U[length(tt), , ] == inits$U0[, 1:n_diag_rec]))
  expect_true(all(y$I[length(tt), , ] == inits$I0[, 1:n_diag_rec]))
  expect_true(all(y$A[length(tt), , ] == inits$A0[, 1:n_diag_rec]))
  expect_true(all(y$S[length(tt), , ] == inits$S0[, 1:n_diag_rec]))
  expect_true(all(y$T[length(tt), , ] == inits$T0[, 1:n_diag_rec]))
  expect_true(5 == inits$t)
  
  ## check that restarting works properly
  tt1 <- seq.int(0, 10)
  y1 <- mod$run(tt1)
  y1 <- mod$transform_variables(y1)
  
  params2 <- model_params(gono_params = gono_params(1)[[1]],
                          init_params = inits, n_diag_rec = n_diag_rec)
  mod2 <- model$new(user = params2, unused_user_action = "ignore")
  y2 <- mod2$run(seq.int(inits$t, 10))
  y2 <- mod2$transform_variables(y2)
  
  expect_equivalent(y1$U[y1$t >= 5, , , drop = FALSE], y2$U, tol = 0.1)
  expect_equivalent(y1$lambda[y1$t >= 5, , drop = FALSE], y2$lambda, tol = 1e-5)
  
  }
  
})


test_that("t_stop is working correctly", {
  ## check with single parameter set
  
  
  for (i in 1:5){
  
  n_diag_rec <- i
  vp <- vax_params_xvwv(vbe = 0, uptake = 1, strategy = "VoD", vea = 1,
                        t_stop = 2 / 365, n_diag_rec = n_diag_rec)
  params <- model_params(gono_params = gono_params(1)[[1]],
                         vax_params = vp, n_diag_rec = n_diag_rec)
  mod <- model$new(user = params, unused_user_action = "ignore")
  tt <- seq.int(0, 5) / 365
  y <- mod$run(tt, )
  y <- mod$transform_variables(y)
  
  expect_equal(diff(apply(y$cum_vaccinated, 1, sum))[tt[-length(tt)] > 2 / 365],
               c(0, 0))
  
  }
  
})

test_that("aggregated time series output correctly", {
  ## check with single parameter set
  
  
  for (i in 1:5){
  
  n_diag_rec <- i
  
  params <- model_params(gono_params = gono_params(1)[[1]], n_diag_rec = n_diag_rec)
  mod <- model$new(user = params, unused_user_action = "ignore")
  tt <- seq.int(0, 5) / 365
  y <- mod$run(tt)
  y <- mod$transform_variables(y)
  expect_equal(y$tot_treated, apply(y$cum_treated, 1, sum))
  expect_equal(y$tot_attended, apply(y$cum_screened, 1, sum) + y$tot_treated)
  
  }
  
})


test_that("time-varying eta works as expected", {
  
  
  for (i in 1:5){
  
  n_diag_rec <- i
  
  gono_pars <- gono_params(1)[[1]]
  params <- model_params(gono_params = gono_pars, n_diag_rec = n_diag_rec)
  params$tt <- c(0, 1, 2)
  gono_pars$eta <- 1
  params$eta_l_t <- params$eta_h_t <- gono_pars$eta * c(1, 2, 2)
  params$beta_t <- rep(gono_pars$beta[1], 3)
  mod <- model$new(user = params, unused_user_action = "ignore")
  tt <- seq(0, 2, by = 1 / 12)
  y <- mod$run(tt)
  y <- mod$transform_variables(y)
  plot(tt[-1] + 2009, diff(rowSums(y$cum_screened)))
  
  expect_equal(y$eta[, 1], approx(params$tt, params$eta_l_t, tt)$y)
  expect_equal(y$eta[, 2], approx(params$tt, params$eta_h_t, tt)$y)
  
  # check can vary wrt group
  params$eta_l_t[] <- gono_pars$eta
  mod <- model$new(user = params, unused_user_action = "ignore")
  y <- mod$run(tt)
  y <- mod$transform_variables(y)
  matplot(apply(y$cum_screened, 2, diff), type = "l")
  
  expect_equal(y$eta[, 1], approx(params$tt, params$eta_l_t, tt)$y)
  expect_equal(y$eta[, 2], approx(params$tt, params$eta_h_t, tt)$y)
  
  # check can switch off screening in a group
  params$eta_l_t[] <- 0
  mod <- model$new(user = params, unused_user_action = "ignore")
  y1 <- mod$run(tt)
  y1 <- mod$transform_variables(y1)
  expect_equal(sum(y1$cum_screened[, 1, ]), 0)
  expect_true(all(y1$cum_screened[-1, 2, ] > 0))
  }
  
  
})


test_that("for n_diag_rec > 1, the number treated = the number recorded
          as diagnosed", {
            
            
            for (i in 2:5){

            n_diag_rec = 5
          
            params <- model_params(gono_params = gono_params(1)[[1]], n_diag_rec = n_diag_rec)
            params$wd = matrix(data =  rep(0, n_diag_rec*n_diag_rec), nrow = n_diag_rec, ncol = n_diag_rec)
            
            
            params$enr = 0
            params$exr = 0
            

            
            mod <- model$new(user = params, unused_user_action = "ignore")
            
            tt <- seq.int(0, 5) / 365
            y <- mod$run(tt)
            y <- mod$transform_variables(y)
            
            

            expect_equivalent(y$cum_treated[, 1, n_diag_rec-1], y$N[, 1, n_diag_rec], tol = 1e-8)
            expect_equivalent(y$cum_treated[, 2, n_diag_rec-1], y$N[, 2, n_diag_rec], tol = 1e-8)
            
            
            for (ii in 2:n_diag_rec){
              
              if (ii < n_diag_rec){
                expect_equivalent(diff(y$cum_treated[, 2, ii-1]), diff(rowSums(y$N[, 2, ii:n_diag_rec])),  tol = 1e-12)
              } else{
                expect_equivalent(diff(y$cum_treated[, 2, ii-1]), diff(y$N[, 2, n_diag_rec]),  tol = 1e-12)
                
              }
              
            }
            

            

            
            
            }
          
            
          })






