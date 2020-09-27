context("model (check)")

test_that("there are no infections when beta is 0", {
  params <- model_params(gono_params = gono_params(1))
  params$beta_t[] <- 0
  mod <- model(user = params)

  tt <- seq.int(0, 5) / 365
  y <- mod$run(tt)
  y <- mod$transform_variables(y)
  expect_true(all(y$I == 0))
  expect_true(all(y$S == 0))
  expect_true(all(y$cum_incid == 0))
})

test_that("there are no symptomatic infections when psi = 0", {
  params <- model_params(gono_params = gono_params(1))
  params$psi <- 0
  mod <- model(user = params)
  tt <- seq.int(0, 5) / 365
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  expect_true(any(y$I == 0))
  expect_true(all(y$S == 0))
  expect_true(all(y$cum_diag_s == 0))
  expect_true(all(y$cum_diag_a[-1, , ] > 0))
})

test_that("there are no asymptomatic infections when psi = 1", {
  params <- model_params(gono_params = gono_params(1))
  params$psi <- 1
  params$S0[, ] <- params$A0[, ]
  params$A0[, ] <- 0
  mod <- model(user = params)
  tt <- seq.int(0, 5) / 365
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  expect_true(any(y$I == 0))
  expect_true(all(y$A == 0))
  expect_true(all(y$cum_diag_a == 0))
  expect_true(all(y$cum_diag_s[-1, , ] > 0))
})

test_that("there are no infections when A0 = 0", {
  params <- model_params(gono_params = gono_params(1))
  params$A0[, ] <- 0
  mod <- model(user = params)
  tt <- seq.int(0, 5) / 365
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)

  expect_true(all(y$I == 0))
  expect_true(all(y$A == 0))
  expect_true(all(y$S == 0))
})

test_that("no-one is treated when mu and eta = 0", {
  params <- model_params(gono_params = gono_params(1))
  params$mu <- params$eta_t[] <- 0
  mod <- model(user = params)

  tt <- seq.int(0, 5) / 365
  y <- mod$run(tt)
  y <- mod$transform_variables(y)
  expect_true(all(y$T == 0))
  expect_true(all(y$cum_treated == 0))
})

test_that("the foi is calculated correctly", {
  params <- model_params(gono_params = gono_params(1))
  expect_true(length(params$beta) > 0)
  mod <- model(user = params)
  tt <- seq.int(0, 5) / 365
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  # unpack parameters
  pL <- params$p[1]
  pH <- params$p[2]
  NL <- sum(y$N[1, 1, ])
  NH <- sum(y$N[1, 2, ])
  C <- y$I + y$A + y$S
  CL <- C[, 1, ]
  CH <- C[, 2, ]
  eps <- params$epsilon
  beta <- params$beta
  np <- pL * NL + pH * NH
  # calculate FOI
  foi_LL <- pL * beta * (eps + (1 - eps) * pL * NL / np) * CL / NL
  foi_LH <- pL * pH * beta * (1 - eps) / np * CH
  foi_HL <- pH * pL * beta * (1 - eps) / np * CL
  foi_HH <- pH * beta * (eps + (1 - eps) * pH * NH / np) * CH / NH
  # test
  expect_equal(y$foi[, 1, 1], foi_LL)
  expect_equal(y$foi[, 1, 2], foi_LH)
  expect_equal(y$foi[, 2, 1], foi_HL)
  expect_equal(y$foi[, 2, 2], foi_HH)
})

test_that("Bex model runs with no vaccination", {
  tt <- seq.int(0, 5) / 365
  params0 <- model_params(gono_params = gono_params(1))
  mod0 <- model(user = params0)
  y0 <- mod0$run(t = tt)
  y0 <- mod0$transform_variables(y0)

  params1 <- model_params(gono_params = gono_params(1),
                           vax_params = vax_params1(ve = 0))
  mod1 <- model(user = params1)
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

test_that("Bex model runs with vbe", {
  tt <- seq.int(0, 2) / 365
  # with perfect efficacy
  params <- model_params(gono_params = gono_params(1),
                            vax_params = vax_params1(ve = 1, eff = 1))
  mod <- model(user = params)
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  expect_equal(y$U,
               array(c(505266, 505270.782413465, 505276.010170825,
                       89303, 89292.2121354773, 89281.4753785608,
                       0, 27.9444015916977, 55.8871954695537,
                       0, 4.93136498677018, 9.862446259333,
                       0, 3.82796084566103e-05, 0.000153112453356132,
                       0, 6.75495729408945e-06, 2.70176968895868e-05),
               dim = c(3, 2, 3)))
  # check some people are being vaccinated
  expect_true(all(y$U[-1, , 2] > 0))
  # check no compartments are leaking
  expect_true(all(apply(y$N, c(1, 2), sum) - 6e5 < 1e-6))
  # check there are infections in unvaccinated group
  expect_false(all(y$I[, , 1] == 0))
  expect_false(all(y$A[, , 1] == 0))
  expect_false(all(y$S[, , 1] == 0))
  expect_false(all(y$T[, , 1] == 0))
  # check there are no infections in vaccinated group
  expect_true(all(y$I[, , 2] == 0))
  expect_true(all(y$A[, , 2] == 0))
  expect_true(all(y$S[, , 2] == 0))
  expect_true(all(y$T[, , 2] == 0))
})

test_that("Check vaccination on screening in Bex model", {
  tt <- seq.int(0, 2) / 365
  # with perfect efficacy
  params <- model_params(gono_params = gono_params(1),
                           vax_params = vax_params1(ve = 0, vs = 1, eff = 1))
  mod <- model(user = params)
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  expect_equal(y$U,
               array(c(505266, 504689.65898183, 504114.492415873,
                       89303, 89189.5072358554, 89076.2208839755,
                       0, 609.068445893679, 1217.40743007721,
                       0, 107.642503971552, 215.142050698242,
                       0, 0.000834155028318817, 0.0033338806826574,
                       0, 0.000147420094481985, 0.000589147144744636),
                     dim = c(3L, 2L, 3L)))
  # check some people are being vaccinated
  expect_true(all(y$U[-1, , 2] > 0))
  expect_true(all(y$cum_vaccinated[-1, , 1] > 0))
  expect_true(all(y$cum_vaccinated[, , 2] == 0))
  # check all those treated were vaccinated
  expect_true(all(y$cum_vaccinated[, , 1] == y$cum_screened[, , 1]))
  # check no compartments are leaking
  expect_true(all(apply(y$N, 1, sum) - 6e5 < 1e-6))
  # check there are infections in unvaccinated group
  expect_false(all(y$I[, , 1] == 0))
  expect_false(all(y$A[, , 1] == 0))
  expect_false(all(y$S[, , 1] == 0))
  expect_false(all(y$T[, , 1] == 0))
  # check there are no infections in vaccinated group
  expect_true(all(y$I[, , 2] == 0))
  expect_true(all(y$A[, , 2] == 0))
  expect_true(all(y$S[, , 2] == 0))
  expect_true(all(y$T[, , 2] == 0))
})

test_that("Check vaccination on diagnosis in Bex model", {
  tt <- seq.int(0, 2) / 365
  # with perfect efficacy
  params <- model_params(gono_params = gono_params(1),
                           vax_params = vax_params1(ve = 0, vd = 1, eff = 1))
  mod <- model(user = params)
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  expect_equal(y$U,
               array(c(505266, 505298.314652489, 505330.328904702,
                       89303, 89297.0812095395, 89291.0874630527,
                       0, 0.412133420953426, 1.56834473390679,
                       0, 0.06199394863926, 0.249176144196752,
                       0, 3.81264799241714e-07, 2.93702682328813e-06,
                       0, 5.68712159308333e-08, 4.53390060409456e-07),
                     dim = c(3L, 2L, 3L)))
  # check some people are being vaccinated
  expect_true(all(y$U[-1, , 2] > 0))
  expect_true(all(y$cum_vaccinated[-1, , 1] > 0))
  expect_true(all(y$cum_vaccinated[, , 2] == 0))
  # check all those treated were vaccinated
  expect_true(all(y$cum_vaccinated == y$cum_treated))
  # check no compartments are leaking
  expect_true(all(apply(y$N, 1, sum) - 6e5 < 1e-6))
  # check there are infections in unvaccinated group
  expect_false(all(y$I[, , 1] == 0))
  expect_false(all(y$A[, , 1] == 0))
  expect_false(all(y$S[, , 1] == 0))
  expect_false(all(y$T[, , 1] == 0))
  # check there are no infections in vaccinated group
  expect_true(all(y$I[, , 2] == 0))
  expect_true(all(y$A[, , 2] == 0))
  expect_true(all(y$S[, , 2] == 0))
  expect_true(all(y$T[, , 2] == 0))
})

test_that("can initialise after time 0", {

  ## check with single parameter set
  params <- model_params(gono_params = gono_params(1))
  mod <- model(user = params)

  tt <- seq.int(0, 5) / 365
  y <- mod$run(tt)
  y <- mod$transform_variables(y)

  inits <- restart_params(y, n_vax = 1)

  expect_true(all(y$U[length(tt), , ] == inits$U0[, 1]))
  expect_true(all(y$I[length(tt), , ] == inits$I0[, 1]))
  expect_true(all(y$A[length(tt), , ] == inits$A0[, 1]))
  expect_true(all(y$S[length(tt), , ] == inits$S0[, 1]))
  expect_true(all(y$T[length(tt), , ] == inits$T0[, 1]))
})

test_that("t_stop is working correctly", {
  ## check with single parameter set
  params <- model_params(gono_params = gono_params(1),
                         vax_params = vax_params1(ve = 0, vd = 1, eff = 1,
                                                  t_stop = 2 / 365))
  mod <- model(user = params)
  tt <- seq.int(0, 5) / 365
  y <- mod$run(tt, )
  y <- mod$transform_variables(y)

  expect_equal(diff(apply(y$cum_vaccinated, 1, sum))[tt[-length(tt)] > 2 / 365],
               c(0, 0))

})

test_that("aggregated time series output correctly", {
  ## check with single parameter set
  params <- model_params(gono_params = gono_params(1))
  mod <- model(user = params)
  tt <- seq.int(0, 5) / 365
  y <- mod$run(tt)
  y <- mod$transform_variables(y)
  expect_equal(y$tot_treated, apply(y$cum_treated, 1, sum))
  expect_equal(y$tot_attended, apply(y$cum_screened, 1, sum) + y$tot_treated)

})

test_that("time-varying eta works as expected", {
  params <- model_params(gono_params = gono_params(1))
  params$tt <- c(0, 1, 2)
  params$eta_t <- params$eta * c(1, 2, 2)
  params$beta_t <- rep(params$beta, 3)
  mod <- model(user = params)
  tt <- seq.int(0, 2, by = 1 / 12) 
  y <- mod$run(tt)
  y <- mod$transform_variables(y)
  plot(diff(rowSums(y$cum_screened)))
  
}) 