context("dualvax (check)")

test_that("there are no infections when beta is 0", {
  params <- dualvax_params(gono_params = gono_params(1))
  params$beta[] <- 0
  mod <- dualvax(user = params)

  tt <- seq(0, 1, 1 / 365)
  y <- mod$run(tt)
  y <- mod$transform_variables(y)
  expect_true(all(y$I == 0))
  expect_true(all(y$S == 0))
})

test_that("there are no symptomatic infections when psi = 0", {
  params <- dualvax_params(gono_params = gono_params(1))
  params$psi[] <- 0
  mod <- dualvax(user = params)
  tt <- seq(0, 1, 1 / 365)
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)
  expect_true(any(y$I == 0))
  expect_true(all(y$S == 0))
})

test_that("there are no infections when A0 = 0", {
  params <- dualvax_params(gono_params = gono_params(1))
  params$A0[, , ] <- 0
  mod <- dualvax(user = params)
  tt <- seq(0, 1, 1 / 365)
  y <- mod$run(t = tt)
  y <- mod$transform_variables(y)

  expect_true(all(y$I == 0))
  expect_true(all(y$A == 0))
  expect_true(all(y$S == 0))
})

test_that("the foi is calculated correctly", {
  params <- dualvax_params(gono_params = gono_params(2))
  mod <- dualvax(user = params)
  tt <- seq(0, 10 / 365, 1 / 365)
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
