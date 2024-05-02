context("parameters")

test_that("cannot input bad parameters", {
  p <- read_csv(gonovax_file("extdata/gono_params.csv"))[1, ]

  expect_error(check_gono_params(replace(p, "prev_Asl", -1)),
               "'prev_Asl' must be between 0 and 1")
  expect_error(check_gono_params(replace(p, "prev_Ash", 2)),
               "'prev_Ash' must be between 0 and 1")
  expect_error(check_gono_params(replace(p, "epsilon", 2)),
               "'epsilon' must be between 0 and 1")
  expect_error(check_gono_params(replace(p, "psi", 2)),
               "'psi' must be between 0 and 1")
  expect_error(check_gono_params(replace(p, "sigma", 0)),
               "'sigma' must be greater than 0")
  expect_error(check_gono_params(replace(p, "nu", 0)),
               "'nu' must be greater than 0")
  expect_error(check_gono_params(replace(p, "mu", 0)),
               "'mu' must be greater than 0")
  expect_error(check_gono_params(replace(p, "rho", 0)),
               "'rho' must be greater than 0")
})

test_that("can select specific parameter sets", {
  p <- read_csv(gonovax_file("extdata/gono_params.csv"))
  gp <- lapply(seq_len(nrow(p)), function(i) transform0(p[i, ]))
  # check that null argument returns all params
  expect_equivalent(gono_params(), gp)
  # check that can extract a single parameter
  expect_equivalent(gono_params(1)[[1]], transform0(p[1, ]))
  # check that can extract multiple parameters
  expect_equivalent(gono_params(2:3), gp[2:3])
  # check that negative parameters will not be returned
  expect_equivalent(gono_params(c(-1, 1)), gp[1])
  # check that cannot extend beyond parameter set
  expect_equivalent(gono_params(c(1, nrow(p) + 1)), gp[1])
})

test_that("vax_map works correctly", {
  # check error when v not length 2
  expect_error(create_vax_map(n_vax = 2, v = rep(0, 3), i_u = 1, i_v = 2))
  expect_error(create_vax_map(n_vax = 2, v = 1, i_u = 1, i_v = 2))

  # check error when v not in 0/1
  expect_error(create_vax_map(n_vax = 2, v = c(2, 2), i_u = 1, i_v = 2))
  expect_error(create_vax_map(n_vax = 2, v = c(0.5, 0.5), i_u = 1, i_v = 2))
  expect_error(create_vax_map(n_vax = 2, v = c(-1, -1), i_u = 1, i_v = 2))

  # check error length(i_u) != length(i_v)
  expect_error(create_vax_map(n_vax = 3, v = c(1, 1), i_u = c(1, 3),
                              i_v = 2))

  # indices of y [group, -to, -from]
  # test onevax XVW map
  y <- create_vax_map(n_vax = 3, v = c(1, 1), i_u = 1, i_v = 2)
  expect_equal(y[1, , ], y[2, , ]) ## same for L and H
  expect_equal(y[, 1, ], -y[, 2, ]) ## everyone leaving X appears in V
  expect_equal(y[1, , 1], c(1, -1, 0)) ## people vaccinated from X -> V
  expect_true(all(y[, , 2] == 0)) ## no-one vaccinated from V
  expect_true(all(y[, , 3] == 0)) ## no-one vaccinated from W
  y <- create_vax_map(n_vax = 3, v = c(0, 1), i_u = 1, i_v = 2)
  expect_true(all(y[1, , ] == 0)) ## no-one vaccinated from L

  # test onevax_xvwv map
  y <- create_vax_map(n_vax = 3, v = c(1, 1), i_u = c(1, 3), i_v = c(2, 2))
  expect_equal(y[1, , ], y[2, , ]) ## same for L and H
  ## everyone leaving X and R appears in V
  expect_equal(y[, 1, ] + y[, 3, ], -y[, 2, ])
  expect_equal(y[1, , 1], c(1, -1, 0)) ## people vaccinated from X -> V
  expect_equal(y[1, , 3], c(0, -1, 1)) ## people vaccinated from R -> V
  expect_true(all(y[, , 2] == 0)) ## no-one vaccinated from V

  y <- create_vax_map(n_vax = 3, v = c(0, 1), i_u = c(1, 3), i_v = c(2, 2))
  expect_true(all(y[1, , ] == 0)) ## no-one vaccinated from L

  # test onevax_xvwr map
  y <- create_vax_map(n_vax = 4, v = c(1, 1), i_u = c(1, 3), i_v = c(2, 4))
  expect_equal(y[1, , ], y[2, , ]) ## same for L and H
  expect_equal(y[, 1, ], -y[, 2, ])   ## everyone leaving X appears in V
  expect_equal(y[, 3, ], -y[, 4, ])  ## everyone leaving W appears in R
  expect_equal(y[1, , 1], c(1, -1, 0, 0)) ## people vaccinated from X -> V
  expect_equal(y[1, , 3], c(0, 0, 1, -1)) ## people vaccinated from W -> R
  expect_true(all(y[, , c(2, 4)] == 0)) ## no-one vaccinated from V or R

  y <- create_vax_map(n_vax = 4, v = c(0, 1), i_u = c(1, 3), i_v = c(2, 4))
  expect_true(all(y[1, , ] == 0)) ## no-one vaccinated from L

  # test onevax_xvwrh map
  y <- create_vax_map(n_vax = 5, v = c(1, 1), i_u = c(1, 3), i_v = c(2, 4))
  expect_equal(y[1, , ], y[2, , ]) ## same for L and H
  expect_equal(y[, 1, ], -y[, 2, ])   ## everyone leaving X appears in V
  expect_equal(y[, 3, ], -y[, 4, ])  ## everyone leaving W appears in R
  expect_equal(y[1, , 1], c(1, -1, 0, 0, 0)) ## people vaccinated from X -> V
  expect_equal(y[1, , 3], c(0, 0, 1, -1, 0)) ## people vaccinated from W -> R
  expect_true(all(y[, , c(2, 4, 5)] == 0)) ## no-one vaccinated from V, R or H

  y <- create_vax_map(n_vax = 5, v = c(0, 1), i_u = c(1, 3), i_v = c(2, 4))
  expect_true(all(y[1, , ] == 0)) ## no-one vaccinated from L

})

test_that("waning map works correctly", {

  # check error negative waning
  expect_error(create_waning_map(n_vax = 2, i_v = 2, i_w = 1, z = -1))
  # check error more than one waning compartment
  expect_error(create_waning_map(n_vax = 3, i_v = 2, i_w = c(1, 3), z = 1))
  # check error length(z) > length(i_v)
  expect_error(create_waning_map(n_vax = 3, i_v = 2, i_w = 3, z = c(1, 0.5)))

  # test onevax map
  y <- create_waning_map(n_vax = 2, i_v = 2, i_w = 1, z = 1)
  expect_equal(rowSums(y), c(1, -1))
  expect_equal(colSums(y), c(0, 0))

  # test onevax waning map
  y <- create_waning_map(n_vax = 3, i_v = 2, i_w = 3, z = 1)
  expect_equal(rowSums(y), c(0, -1, 1))
  expect_equal(colSums(y), c(0, 0, 0))

  y <- create_waning_map(n_vax = 4, i_v = c(2, 4), i_w = 3, z = c(1, 2))
  expect_equal(rowSums(y), c(0, -1, 3, -2))
  expect_equal(colSums(y), c(0, 0, 0, 0))

})

test_that("transform works as expected", {
  p <- read_csv(gonovax_file("extdata/gono_params.csv"))[1, ]
  p$beta2009 <- p$beta
  p$phi_beta <- 0.01
  p$phi_eta <- 0.02
  p$eta_h <- p$eta
  p$omega <- 0.6
  gp <- transform(p)

  # check ratio is fixed between eta_l and eta_h
  expect_equal(gp$eta_l_t / gp$eta_h_t, rep(p$omega, length(gp$tt)))

  i_2020 <- which(gp$tt == gonovax_year(2020))

  ## check increasing to 2020
  expect_equal(diff(gp$beta_t[seq_len(i_2020)]),
               rep(p$beta2009 * p$phi_beta, i_2020 - 1L))
  expect_equal(diff(gp$eta_h_t[seq_len(i_2020)]),
               rep(p$eta_h * p$phi_eta, i_2020 - 1L))
  expect_equal(diff(gp$eta_l_t[seq_len(i_2020)]),
               rep(p$eta_h * p$phi_eta * p$omega, i_2020 - 1L))

  ## check stable after 2020
  expect_true(all(diff(gp$beta_t[-seq_len(i_2020)]) == 0))
  expect_true(all(diff(gp$eta_h_t[-seq_len(i_2020)]) == 0))
  expect_true(all(diff(gp$eta_l_t[-seq_len(i_2020)]) == 0))

  #Check for fix_par_t = FALSE
  gp2 <- transform(p, FALSE)

  expect_equal(diff(gp2$beta_t), rep(p$beta2009 * p$phi_beta,
                                     length(gp2$tt) - 1L))
  expect_equal(diff(gp2$eta_h_t), rep(p$eta_h * p$phi_eta, length(gp2$tt) - 1L))
  expect_equal(diff(gp2$eta_l_t), rep(p$eta_h * p$phi_eta * p$omega,
                                      length(gp2$tt) - 1L))

  ## cannot input bad params
  expect_error(transform(replace(p, "beta2009", 0)),
               "'beta2009' must be greater than 0")
  expect_error(transform(replace(p, "phi_beta", 0)),
               "'phi_beta' must be greater than 0")
  expect_error(transform(replace(p, "phi_eta", 0)),
               "'phi_eta' must be greater than 0")
  expect_error(transform(replace(p, "eta_h", 0)),
               "'eta_h' must be greater than 0")
  expect_error(transform(replace(p, "omega", 1.1)),
               "'omega' must be between 0 and 1")

})

test_that("transform_fixed works as expected", {
  p <- read_csv(gonovax_file("extdata/gono_params_t.csv"))[1, ]

  gp <- transform_fixed(p)

  # check ratio is fixed between eta_l and eta_h
  omega <- gp$eta_l_t / gp$eta_h_t
  expect_equal(omega[1], omega[2])

  ## check stable
  expect_true(all(diff(gp$beta_t) == 0))
  expect_true(all(diff(gp$eta_h_t) == 0))
  expect_true(all(diff(gp$eta_l_t) == 0))

  ## cannot input bad params
  expect_error(transform_fixed(replace(p, "beta", 0)),
               "'beta' must be greater than 0")
  expect_error(transform_fixed(replace(p, "eta_l", 0)),
               "'eta_l' must be greater than 0")
  expect_error(transform_fixed(replace(p, "eta_h", 0)),
               "'eta_h' must be greater than 0")

})

test_that("set_strategy works as expected", {

  # and different strategies vod and vos maps as expected
  expect_equal(set_strategy(),
               list(vod = c(0, 0), vos = c(0, 0), vbe = c(0, 0), vopn = c(0,0)))
  expect_equal(set_strategy(include_vbe =  TRUE),
               list(vod = c(0, 0), vos = c(0, 0), vbe = c(1, 1), vopn = c(0,0)))
  expect_equal(set_strategy("VoD"),
               list(vod = c(1, 1), vos = c(0, 0), vbe = c(0, 0), vopn = c(0,0)))
  expect_equal(set_strategy("VoD(H)"),
               list(vod = c(0, 1), vos = c(0, 0), vbe = c(0, 0), vopn = c(0,0)))
  expect_equal(set_strategy("VoA"),
               list(vod = c(1, 1), vos = c(1, 1), vbe = c(0, 0), vopn = c(0,0)))
  expect_equal(set_strategy("VoA(H)"),
               list(vod = c(0, 1), vos = c(0, 1), vbe = c(0, 0), vopn = c(0,0)))
  expect_equal(set_strategy("VoD(L)+VoA(H)"),
               list(vod = c(1, 1), vos = c(0, 1), vbe = c(0, 0), vopn = c(0,0)))
  expect_equal(set_strategy("VoS"),
               list(vod = c(0, 0), vos = c(1, 1), vbe = c(0, 0), vopn = c(0,0)))
  expect_equal(set_strategy("VoN"),
               list(vod = c(1, 1), vos = c(0, 0), vbe = c(0, 0), vopn = c(1,1)))
  expect_equal(set_strategy("VaH"),
               list(vod = c(1, 1), vos = c(1, 1), vbe = c(0, 0), vopn = c(0,0)))
  expect_equal(set_strategy("VaH+VoN"),
               list(vod = c(1, 1), vos = c(1, 1), vbe = c(0, 0), vopn = c(1,1)))
  
  expect_error(set_strategy("hello"), "strategy not recognised")

})

test_that("initial params works as expected", {
  pars <- c(demographic_params(), gono_params(1)[[1]])
  prev0 <- c(pars$prev_Asl, pars$prev_Ash)

  y <- initial_params(pars)
  expect_equal(sum(unlist(y)), pars$N0)
  expect_equivalent(y$A0 / (y$A0 + y$U0), prev0, tol = 1e-5)
  expect_equal(sum(y$I0 + y$S0 + y$T0), 0)
  expect_equivalent(y$U0, round(pars$N0 * (1 - prev0) * pars$q))

  cov <- c(0.4, 0.6, 0)
  y1 <- initial_params(pars, n_vax = 3, coverage = cov)
  expect_equal(sum(unlist(y1)), pars$N0)
  # check initial prevalence in unvaccinated
  expect_equivalent((y1$A0[, 1] / (y1$A0[, 1] + y1$U0[, 1])), prev0, tol = 1e-5)
  # check initial prevalence is zero in vaccinated
  expect_equal(sum(y1$A0[, -1]), 0)
  expect_equal(sum(y1$I0 + y1$S0 + y1$T0), 0)
  expect_equal(colSums(y1$U0 + y1$A0) / pars$N0, cov)

  expect_error(initial_params(pars, coverage = 2))
  expect_error(initial_params(pars, coverage = -1))
  expect_error(initial_params(pars, n_vax = 3, coverage = c(1, 0)))
  expect_error(initial_params(pars, n_vax = 3, coverage = c(1, 1, 0)))
})
