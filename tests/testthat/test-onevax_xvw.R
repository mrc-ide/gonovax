context("onevax_xvw")

test_that("run_onevax_xvw works correctly", {
  tt <- seq(0, 5)
  gp <- gono_params(1:2)
  y1 <- run_onevax_xvw(tt, gp, vea = 0, dur = 1e3)[[1]]

  # check no-one is offered vaccine with v switched off
  expect_true(all(y1$cum_vaccinated == 0))
  expect_true(all(y1$cum_offered == 0))

  y2 <- run_onevax_xvw(tt, gp, vea = 0, dur = 1e3, vbe = 1)[[1]]
  # check 100% vbe vaccinates all new entrants
  expect_equal(diff(rowSums(y2$cum_vaccinated[, , 1])), rep(12e3, max(tt)))
  # and no-one else
  expect_equal(sum(y2$cum_vaccinated[, , 2:3]), 0)
  # only adolescents are offered vaccine and all accept
  expect_equal(y2$cum_offered, y2$cum_vbe)

  # check can restart
  y2 <- run_onevax_xvw(tt, gp, vea = 0, dur = 1e3, vbe = 1)
  init_params <- lapply(y2, restart_params)
  y3 <- run_onevax_xvw(seq(max(tt), length.out = 2, by = 1),
                       gp, init_params,
                       vea = 0, dur = 1e3, vbe = 1)
  for (i in seq_along(y2)) {
    expect_equal(y2[[i]]$U[length(tt), , ], y3[[i]]$U[1, , ])
    expect_equal(y2[[i]]$I[length(tt), , ], y3[[i]]$I[1, , ])
    expect_equal(y2[[i]]$A[length(tt), , ], y3[[i]]$A[1, , ])
    expect_equal(y2[[i]]$S[length(tt), , ], y3[[i]]$S[1, , ])
    expect_equal(y2[[i]]$T[length(tt), , ], y3[[i]]$T[1, , ])
  }

  uptake <- c(0.5, 1) ## i.e. 50% for 1st param set, 100% for second
  # check VoD is working correctly
  y3e <- run_onevax_xvw(tt, gp, vea = 1, dur = 1e3, strategy = "VoD",
                        uptake = uptake)

  for (i in seq_along(y3e)) {
    # no-one in stratum V or W is vaccinated again
    expect_equal(sum(y3e[[i]]$cum_vaccinated[, , "V"]), 0)
    expect_equal(sum(y3e[[i]]$cum_vaccinated[, , "W"]), 0)
    # uptake % of treated are vaccinated
    expect_equal(y3e[[i]]$cum_vaccinated[, , "X"] / uptake[i],
                 y3e[[i]]$cum_treated[, , "X"])
    # all those treated are offered vaccine
    expect_equal(y3e[[i]]$cum_offered[, , "X"], y3e[[i]]$cum_treated[, , "X"])
    # only those in X are offered vaccine
    expect_equal(sum(y3e[[i]]$cum_offered[, , c("V", "W")]), 0)
    # efficacy is perfect
    expect_equal(sum(y3e[[i]]$cum_treated[, , "V"]), 0)
    # no-one is lost
    expect_equal(apply(y3e[[i]]$N, 1, sum), rep(6e5, 6), tolerance = 1e-5)
  }

  # check VoA is working correctly
  y4e <- run_onevax_xvw(tt, gp, vea = 1, dur = 1e3, strategy = "VoA",
                        uptake = uptake)
  for (i in seq_along(y4e)) {
    # no-one in stratum V or W is vaccinated again
    expect_equal(sum(y4e[[i]]$cum_vaccinated[, , "V"]), 0)
    expect_equal(sum(y4e[[i]]$cum_vaccinated[, , "W"]), 0)

    # uptake % of treated are vaccinated
    expect_equal(y4e[[i]]$cum_vaccinated[, , "X"] / uptake[i],
                 y4e[[i]]$cum_screened[, , "X"] + y4e[[i]]$cum_treated[, , "X"])

    # all those in X treated  or screened are offered vaccine
    expect_equal(y4e[[i]]$cum_offered[, , "X"],
                 y4e[[i]]$cum_treated[, , "X"] + y4e[[i]]$cum_screened[, , "X"])
    # only those in X are offered vaccine
    expect_equal(sum(y4e[[i]]$cum_offered[, , c("V", "W")]), 0)

    # efficacy is perfect
    expect_equal(sum(y4e[[i]]$cum_treated[, , "V"]), 0)
    # no-one is lost
    expect_equal(apply(y4e[[i]]$N, 1, sum), rep(6e5, 6), tolerance = 1e-5)
  }

  # check vaccination targeting
  y5e <- run_onevax_xvw(tt, gp, vea = 1, dur = 1e3, strategy = "VoD(L)+VoA(H)",
                        uptake = uptake)
  for (i in seq_along(y5e)) {
    # no-one in stratum V or W is vaccinated again
    expect_equal(sum(y5e[[i]]$cum_vaccinated[, , "V"]), 0)
    expect_equal(sum(y5e[[i]]$cum_vaccinated[, , "W"]), 0)

    # only treated L are offered vaccine
    expect_equal(y5e[[i]]$cum_offered[, "L", "X"],
                 y5e[[i]]$cum_treated[, "L", "X"])
    # all attending H are offered vaccine
    expect_equal(y5e[[i]]$cum_offered[, "H", "X"],
                 y5e[[i]]$cum_treated[, "H", "X"] +
                   y5e[[i]]$cum_screened[, "H", "X"])

    # check that vaccinated = offered * uptake
    expect_equal(y5e[[i]]$cum_offered * uptake[i], y5e[[i]]$cum_vaccinated)

    # only those in X are offered vaccine
    expect_equal(sum(y4e[[i]]$cum_offered[, , c("V", "W")]), 0)

    # efficacy is perfect
    expect_equal(sum(y4e[[i]]$cum_treated[, , "V"]), 0)
    # no-one is lost
    expect_equal(apply(y4e[[i]]$N, 1, sum), rep(6e5, 6), tolerance = 1e-5)
  }

  # check length of uptake vector must be 1 or length(gp)
  expect_error(run_onevax_xvw(tt, gp, vea = 1, dur = 1e3, strategy = "VbE",
                              uptake = c(0, 0.5, 1)))
  expect_error(run_onevax_xvw(tt, gp, vea = c(0, 1, 2), dur = 1e3,
                              strategy = "VbE", uptake = 1))
  expect_error(run_onevax_xvw(tt, gp, vea = 1, dur = c(0, 1e2, 1e3),
                              strategy = "VbE", uptake = 1))
})

test_that("vaccine effects work as expected", {
  tt <- seq(0, 5)
  gp <- gono_params(1:2)

  ## check perfect protection against symptoms works
  y1 <- run_onevax_xvw(tt, gp, ves = 1, dur = 1,
                       uptake = 1, strategy = "VoA")[[1]]

  expect_true(all(y1$S[, , 2] == 0))
  expect_true(all(y1$S[-1, , 3] > 1e-6))
  expect_true(all(y1$A[-1, , ] > 1e-6))

  ## check protection against duration works
  # for perfect ved (= 1), the maximum rate of movement A->U is the rate of
  # care seeking (mu)
  # if same number of people start in A as in S and parameters
  # lambda and eta = 0 (i.e other routes into A & S and out of A respectively)
  # the rate of movement of people A -> U and S -> T should be the same

  # replace lambda, eta = 0
  elements_to_replace <- c("eta_l_t", "eta_h_t", "beta_t")
  for (i in 1:2) {
    for (element in elements_to_replace) {
      gp[[i]][[element]] <- rep(0, times = length(gp[[1]][[element]]))
    }
  }

  # set starting conditions
  pars <- c(demographic_params(), gp)
  n_vax <- 3
  U0 <- I0 <- A0 <- S0 <- T0 <- array(0, c(2, n_vax))

  #half of population in A0 and half in S0
  #everyone begins vaccinated
  N0 <- pars$N0 * outer(pars$q, c(0, 1, 0))
  A0 <- pars$N0 * outer(pars$q, c(0, 0.5, 0))
  S0 <- pars$N0 * outer(pars$q, c(0, 0.5, 0))

  init_params <- list(U0 = U0, I0 = I0, A0 = A0, S0 = S0, T0 = T0)
  init_params_list <- list(init_params, init_params)

  y2 <- run_onevax_xvw(tt, gono_params = gp, init_params = init_params_list,
                       ved = 1, dur = 1e99)

  for (i in 1:2){
    expect_equal(y2[[i]]$A[, , 2], y2[[i]]$S[, , 2]) # 2 = V stratum
  }

})

test_that("can set initial coverage", {
  tt <- seq(0, 5)
  gp <- gono_params(1:2)
  cov <- 0.33

  ## check perfect protection against symptoms works
  y <- run_onevax_xvw(tt, gp, ves = 1, dur = 1e3, coverage = 0.33)[[1]]
  expect_equivalent(apply(y$N[, , -1], c(1, 2), sum) / apply(y$N, c(1, 2), sum),
                    rep(cov, length(tt) * 2))

  expect_error(run_onevax_xvw(tt, gp, ves = 1, dur = 1e3, coverage = c(0.1, 0)),
               "'coverage' must be a scalar")
})

test_that("can set n_AU conditional statement works as expected", {
  #when mu = 0, defaults to nu
  #and when ved = 0, and mu not = 0, mu's cancel out, and rate is nu
  #so in both cases movement A->U should be the same

  tt <- seq(0, 5)
  gp_mu_og <- gono_params(1:2)
  gp_mu_zero <- gono_params(1:2)

  # replace parameters with zero
  elements_to_replace <- c("eta_l_t", "eta_h_t", "beta_t")
  for (i in 1:2) {
    for (element in elements_to_replace) {
      gp_mu_og[[i]][[element]] <- rep(0,
                                      times = length(gp_mu_og[[1]][[element]]))
    }
  }

  elements_to_replace <- c("eta_l_t", "eta_h_t", "beta_t", "mu")
  for (i in 1:2) {
    for (element in elements_to_replace) {
      gp_mu_zero[[i]][[element]] <- rep(0,
                                      times = length(gp_mu_zero[[1]][[element]]))
    }
  }

  # set starting conditions
  pars <- c(demographic_params(), gp)
  n_vax <- 3
  U0 <- I0 <- A0 <- S0 <- T0 <- array(0, c(2, n_vax))

  #all of population in A0 because we are only interested in movement out of A
  N0 <- pars$N0 * outer(pars$q, c(0, 1, 0))
  A0 <- pars$N0 * outer(pars$q, c(0, 1, 0))

  init_params <- list(U0 = U0, I0 = I0, A0 = A0, S0 = S0, T0 = T0)
  init_params_list <- list(init_params, init_params)

  y1 <- run_onevax_xvw(tt = tt, gono_params = gp_mu_zero,
                             init_params = init_params_list,
                             vea = 0, vei = 0, ved = 0, ves = 0,
                             dur = 1e99)

  y2 <- run_onevax_xvw(tt = tt, gono_params = gp_mu_og,
                             init_params = init_params_list,
                             vea = 0, vei = 0, ved = 0, ves = 0,
                             dur = 1e99)

  for (i in 1:2){
    expect_equal(y1[[i]], y2[[i]])
  }

  #adding in ved when mu = 0, makes no difference
  y3 <- run_onevax_xvw(tt = tt, gono_params = gp_mu_zero,
                             init_params = init_params_list,
                             vea = 0, vei = 0, ved = 0.5, ves = 0,
                             dur = 1e99)
  for (i in 1:2){
    expect_equal(y3[[i]], y2[[i]])
  }

  #adding in ved when mu is not 0, does make a difference
  y4 <- run_onevax_xvw(tt = tt, gono_params = gp_mu_og,
                             init_params = init_params_list,
                             vea = 0, vei = 0, ved = 0.5, ves = 0,
                             dur = 1e99)

  for (i in 1:2){
    expect_error(expect_equal(y3[[i]], y4[[i]]))
  }

})
