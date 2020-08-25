context("model (check)")

test_that("there are no infections when beta is 0", {
  params <- model_params(gono_params = gono_params(1))
  params$beta <- 0
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
  params$mu <- params$eta <- 0
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
  tt <- seq.int(0, 5) / 365
  # with perfect efficacy
  params <- model_params(gono_params = gono_params(1),
                            vax_params = vax_params1(ve = 1, eff = 1))
  mod <- model(user = params)
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  expect_equal(y$U,
               array(c(
               c(505266, 505270.782451745, 505276.01032394,
                 505281.584362789, 505287.421644038, 505293.45181186),
               c(89303, 89292.2121422323, 89281.4754055862,
                 89270.817547954, 89260.2713849845, 89249.867453731),
               c(0, 27.9444015916977, 55.8871954695537,
                 83.8283817260638, 111.767960453719, 139.705931745003),
               c(0, 4.93136498677018, 9.862446259333,
                 14.7932438340113, 19.7237577271268,24.6539879550006)), 
               dim = c(6, 2, 2)))
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
  tt <- seq.int(0, 5) / 365
  # with perfect efficacy
  params <- model_params(gono_params = gono_params(1),
                           vax_params = vax_params1(ve = 0, vs = 1, eff = 1))
  mod <- model(user = params)
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  expect_equal(y$U,
               array(c(505266, 504689.659815985, 504114.495749757,
                       503540.407180298, 502967.309913545, 502395.132403755,
                       89303, 89189.5073832754, 89076.2214731302,
                       88963.1699414225, 88850.3854164374, 88737.8981445319, 
                       0, 609.06844589368, 1217.4074300772,
                       1825.01835049653, 2431.90249401664, 3038.06105382112,
                       0, 107.642503971552, 215.142050698235,
                       322.498913374029, 429.713402534724, 536.78586712584),
                     dim = c(6L, 2L, 2L)))
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
  tt <- seq.int(0, 5) / 365
  # with perfect efficacy
  params <- model_params(gono_params = gono_params(1),
                           vax_params = vax_params1(ve = 0, vd = 1, eff = 1))
  mod <- model(user = params)
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  expect_equal(y$U,
               array(c(505266, 505298.314652871, 505330.328907639,
                       505362.048127296, 505393.477204426, 505424.620601876,
                       89303, 89297.0812095964, 89291.0874635061,
                       89285.032255699, 89278.9271648226, 89272.7819663311,
                       0, 0.412133420953426, 1.56834473675643,
                       3.36401902671277, 5.71134154251699, 8.53549599287203,
                       0, 0.0619939486392602, 0.249176151802269,
                       0.575817426283736, 1.06316731455961, 1.73200597450277),
                     dim = c(6L, 2L, 2L)))
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
