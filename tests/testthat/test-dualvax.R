context("dualvax (check)")

test_that("there are no infections when beta is 0", {
  params <- dualvax_params(gono_params = gono_params(1))
  params$beta[] <- 0
  mod <- dualvax(user = params)

  tt <- seq.int(0, 10) / 365
  y <- mod$run(tt)
  y <- mod$transform_variables(y)
  expect_true(all(y$I == 0))
  expect_true(all(y$S == 0))
  expect_true(all(y$cum_incid == 0))
})

test_that("there are no symptomatic infections when psi = 0", {
  params <- dualvax_params(gono_params = gono_params(1))
  params$psi[] <- 0
  mod <- dualvax(user = params)
  tt <- seq.int(0, 10) / 365
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  expect_true(any(y$I == 0))
  expect_true(all(y$S == 0))
})

test_that("there are no infections when A0 = 0", {
  params <- dualvax_params(gono_params = gono_params(1))
  params$A0[, , ] <- 0
  mod <- dualvax(user = params)
  tt <- seq.int(0, 10) / 365
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)

  expect_true(all(y$I == 0))
  expect_true(all(y$A == 0))
  expect_true(all(y$S == 0))
})

test_that("no-one is treated when mu and eta = 0", {
  params <- dualvax_params(gono_params = gono_params(1))
  params$mu[] <- params$eta[] <- 0
  mod <- dualvax(user = params)

  tt <- seq.int(0, 10) / 365
  y <- mod$run(tt)
  y <- mod$transform_variables(y)
  expect_true(all(y$T == 0))
  expect_true(all(y$cum_treated == 0))
})

test_that("the foi is calculated correctly", {
  params <- dualvax_params(gono_params = gono_params(2))
  mod <- dualvax(user = params)
  tt <- seq.int(0, 10) / 365
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  # unpack parameters
  pL <- params$p[1]
  pH <- params$p[2]
  NL <- sum(y$N[1, 1, 1, ])
  NH <- sum(y$N[1, 1, 2, ])
  C <- y$I + y$A + y$S
  CL <- t(C[, , 1, ])
  CH <- t(C[, , 2, ])
  eps <- params$epsilon
  beta <- params$beta
  np <- pL * NL + pH * NH
  # calculate FOI
  foi_LL <- pL * beta * (eps + (1 - eps) * pL * NL / np) * CL / NL
  foi_LH <- pL * pH * beta * (1 - eps) / np * CH
  foi_HL <- pH * pL * beta * (1 - eps) / np * CL
  foi_HH <- pH * beta * (eps + (1 - eps) * pH * NH / np) * CH / NH
  # test
  expect_equal(y$foi[, , 1, 1], t(foi_LL))
  expect_equal(y$foi[, , 1, 2], t(foi_LH))
  expect_equal(y$foi[, , 2, 1], t(foi_HL))
  expect_equal(y$foi[, , 2, 2], t(foi_HH))
})

test_that("Bex model runs with no vaccination", {
  tt <- seq.int(0, 10) / 365
  params0 <- dualvax_params(gono_params = gono_params(1))
  mod0 <- dualvax(user = params0)
  y0 <- mod0$run(t = tt)
  y0 <- mod0$transform_variables(y0)

  params1 <- dualvax_params(gono_params = gono_params(1),
                           vax_params = vax_params1(ve = 0))
  mod1 <- dualvax(user = params1)
  y1 <- mod1$run(t = tt)
  y1 <- mod1$transform_variables(y1)

  # check that nil vaccination gives same results as before
  expect_true(all(y1$U[, , , 1, drop = FALSE] == y0$U))
  expect_true(all(y1$I[, , , 1, drop = FALSE] == y0$I))
  expect_true(all(y1$A[, , , 1, drop = FALSE] == y0$A))
  expect_true(all(y1$S[, , , 1, drop = FALSE] == y0$S))
  expect_true(all(y1$T[, , , 1, drop = FALSE] == y0$T))

  expect_true(all(y1$N[, , , 2] == 0))
  expect_true(all(apply(y1$N, c(1, 2), sum) - 6e5 < 1e-6))
})

test_that("Bex model runs with vbe", {
  tt <- seq.int(0, 10) / 365
  # with perfect efficacy
  params <- dualvax_params(gono_params = gono_params(1),
                            vax_params = vax_params1(ve = 1, eff = 1))
  mod <- dualvax(user = params)
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  # check some people are being vaccinated
  expect_true(all(y$U[-1, , , 2] > 0))
  # check no compartments are leaking
  expect_true(all(apply(y$N, c(1, 2), sum) - 6e5 < 1e-6))
  # check there are infections in unvaccinated group
  expect_false(all(y$I[, , , 1] == 0))
  expect_false(all(y$A[, , , 1] == 0))
  expect_false(all(y$S[, , , 1] == 0))
  expect_false(all(y$T[, , , 1] == 0))
  # check there are no infections in vaccinated group
  expect_true(all(y$I[, , , 2] == 0))
  expect_true(all(y$A[, , , 2] == 0))
  expect_true(all(y$S[, , , 2] == 0))
  expect_true(all(y$T[, , , 2] == 0))
})

test_that("Check vaccination on screening in Bex model", {
  tt <- seq.int(0, 10) / 365
  # with perfect efficacy
  params <- dualvax_params(gono_params = gono_params(1),
                           vax_params = vax_params1(ve = 0, vs = 1, eff = 1))
  mod <- dualvax(user = params)
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  # check some people are being vaccinated
  expect_true(all(y$U[-1, , , 2] > 0))
  expect_true(all(y$cum_vaccinated[-1, , , 1] > 0))
  expect_true(all(y$cum_vaccinated[ , , , 2] == 0))
  # check all those treated were vaccinated
  expect_true(all(y$cum_vaccinated[, , , 1] == y$cum_screened[, , , 1]))
  # check no compartments are leaking
  expect_true(all(apply(y$N, c(1, 2), sum) - 6e5 < 1e-6))
  # check there are infections in unvaccinated group
  expect_false(all(y$I[, , , 1] == 0))
  expect_false(all(y$A[, , , 1] == 0))
  expect_false(all(y$S[, , , 1] == 0))
  expect_false(all(y$T[, , , 1] == 0))
  # check there are no infections in vaccinated group
  expect_true(all(y$I[, , , 2] == 0))
  expect_true(all(y$A[, , , 2] == 0))
  expect_true(all(y$S[, , , 2] == 0))
  expect_true(all(y$T[, , , 2] == 0))
})

test_that("Check vaccination on diagnosis in Bex model", {
  tt <- seq.int(0, 10) / 365
  # with perfect efficacy
  params <- dualvax_params(gono_params = gono_params(1),
                           vax_params = vax_params1(ve = 0, vd = 1, eff = 1))
  mod <- dualvax(user = params)
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  # check some people are being vaccinated
  expect_true(all(y$U[-1, , , 2] > 0))
  expect_true(all(y$cum_vaccinated[-1, , , 1] > 0))
  expect_true(all(y$cum_vaccinated[ , , , 2] == 0))
  # check all those treated were vaccinated
  expect_true(all(y$cum_vaccinated == y$cum_treated))
  # check no compartments are leaking
  expect_true(all(apply(y$N, c(1, 2), sum) - 6e5 < 1e-6))
  # check there are infections in unvaccinated group
  expect_false(all(y$I[, , , 1] == 0))
  expect_false(all(y$A[, , , 1] == 0))
  expect_false(all(y$S[, , , 1] == 0))
  expect_false(all(y$T[, , , 1] == 0))
  # check there are no infections in vaccinated group
  expect_true(all(y$I[, , , 2] == 0))
  expect_true(all(y$A[, , , 2] == 0))
  expect_true(all(y$S[, , , 2] == 0))
  expect_true(all(y$T[, , , 2] == 0))
})
