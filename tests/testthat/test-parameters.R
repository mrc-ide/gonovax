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
  expect_equal(gono_params(), gp)
  # check that can extract a single parameter
  expect_equal(gono_params(1)[[1]], transform0(p[1, ]))
  # check that can extract multiple parameters
  expect_equal(gono_params(2:3), gp[2:3])
  # check that negative parameters will not be returned
  expect_equal(gono_params(c(-1, 1)), gp[1])
  # check that cannot extend beyond parameter set
  expect_equal(gono_params(c(1, nrow(p) + 1)), gp[1])
})

test_that("vax_map works correctly", {
  # check error when input not matrix
  expect_error(create_vax_map(n_vax = 2, v = rep(0, 3), i_u = 1, i_v = 2))

  # check error when uptake not in 0-1
  v <- matrix(c(rep(2, 2), rep(0, 2)), nrow = 2, byrow = TRUE)
  expect_error(create_vax_map(n_vax = 2, v = v, i_u = 1, i_v = 2))

  v <- matrix(c(rep(-1, 2), rep(0, 2)), nrow = 2, byrow = TRUE)
  expect_error(create_vax_map(n_vax = 2, v = v, i_u = 1, i_v = 2))

  # check error length(i_u) != length(i_v)
  expect_error(create_vax_map(n_vax = 3, v = c(0.1, 0.2), i_u = c(1, 3),
                              i_v = 2))
  # check error length(i_u) != length(i_v)
  expect_error(create_vax_map(n_vax = 2, v = c(0.1, 0.2), i_u = 1, i_v = 3))

  # indices of y [group, -to, -from]
  # test onevax map

  #change vbe input to matrix format
  v <- 0.1
  v <- matrix(c(rep(v, 2), rep(0, 2)), nrow = 2, byrow = TRUE)

  y <- create_vax_map(n_vax = 2, v = v, i_u = 1, i_v = 2)
  expect_equal(y[1, , 1], c(0.1, -0.1))
  expect_true(all(y[, , 2] == 0))

  # test onevax_xvwv waning map
  y <- create_vax_map(n_vax = 3, v = v, i_u = c(1, 3), i_v = c(2, 2))

  expect_equal(y[1, , 1], c(0.1, -0.1, 0))
  expect_equal(y[1, , 3], c(rep(0, 3)))
  expect_true(all(y[, , 2] == 0))

  p <- set_strategy(strategy = "VoD", primary_uptake = 0.1)
  y <- create_vax_map(n_vax = 3, v = p$vod, i_u = c(1, 3), i_v = c(2, 2))

  expect_equal(y[1, , 1], c(0.1, -0.1, 0))
  expect_equal(y[1, , 3], c(0, -0.1, 0.1))
  expect_true(all(y[, , 2] == 0))

  # test onevax_xvwr waning map
  y <- create_vax_map(n_vax = 4, v = v, i_u = c(1, 3), i_v = c(2, 4))

  expect_equal(y[1, , 1], c(0.1, -0.1, 0, 0))
  expect_true(all(y[, , 2] == 0))
  expect_true(all(y[, , 3] == 0))
  expect_true(all(y[, , 4] == 0))

  y <- create_vax_map(n_vax = 4, v = p$vod, i_u = c(1, 3), i_v = c(2, 4))

  expect_equal(y[1, , 1], c(0.1, -0.1, 0, 0))
  expect_equal(y[2, , 1], c(0.1, -0.1, 0, 0))
  expect_equal(y[1, , 3], c(0, 0, 0.1, -0.1))
  expect_equal(y[2, , 3], c(0, 0, 0.1, -0.1))

  expect_true(all(y[, , 2] == 0))
  expect_true(all(y[, , 4] == 0))

  # test onevax_xvwrh waning map
  y <- create_vax_map(n_vax = 5, v = v, i_u = c(1, 3), i_v = c(2, 4))

  expect_equal(y[1, , 1], c(0.1, -0.1, 0, 0, 0))
  expect_true(all(y[, , 2] == 0))
  expect_true(all(y[, , 3] == 0))
  expect_true(all(y[, , 4] == 0))
  expect_true(all(y[, , 5] == 0))

  y <- create_vax_map(n_vax = 5, v = p$vod, i_u = c(1, 3), i_v = c(2, 4))

  expect_equal(y[1, , 1], c(0.1, -0.1, 0, 0, 0))
  expect_equal(y[2, , 1], c(0.1, -0.1, 0, 0, 0))
  expect_equal(y[1, , 3], c(0, 0, 0.1, -0.1, 0))
  expect_equal(y[2, , 3], c(0, 0, 0.1, -0.1, 0))

  expect_true(all(y[, , 2] == 0))
  expect_true(all(y[, , 4] == 0))
  expect_true(all(y[, , 5] == 0))

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

  # different primary and booster uptakes work as expected
  # and different strategies vod and vos maps as expected

  i <- 0.234
  j <- 0.468
  expect_equal(set_strategy("VbE", i, j),
               list(vod = matrix(c(rep(0, 4)), nrow = 2),
                    vos = matrix(c(rep(0, 4)), nrow = 2)))
  expect_equal(set_strategy("VoD", i, j),
               list(vod = matrix(c(i, i, j, j), nrow = 2, byrow = TRUE),
                    vos = matrix(c(rep(0, 4)), nrow = 2)))
  expect_equal(set_strategy("VoD(H)", i, j),
               list(vod = matrix(c(0, i, 0, j), nrow = 2, byrow = TRUE),
                    vos = matrix(c(rep(0, 4)), nrow = 2)))
  expect_equal(set_strategy("VoA", i, j),
               list(vod = matrix(c(i, i, j, j), nrow = 2, byrow = TRUE),
                    vos = matrix(c(i, i, j, j), nrow = 2, byrow = TRUE)))
  expect_equal(set_strategy("VoA(H)", i, j),
               list(vod = matrix(c(0, i, 0, j), nrow = 2, byrow = TRUE),
                    vos = matrix(c(0, i, 0, j), nrow = 2, byrow = TRUE)))
  expect_equal(set_strategy("VoD(L)+VoA(H)", i, j),
               list(vod = matrix(c(i, i, j, j), nrow = 2, byrow = TRUE),
                    vos = matrix(c(0, i, 0, j), nrow = 2, byrow = TRUE)))
  expect_equal(set_strategy("VoS", i, j),
               list(vod = matrix(c(rep(0, 4)), nrow = 2),
                    vos = matrix(c(i, i, j, j), nrow = 2, byrow = TRUE)))

  expect_error(set_strategy("hello", i), "strategy not recognised")
  expect_error(set_strategy("VbE", c(i, i)), "uptake must be length 1")
  
  # booster_uptake defaults to primary_uptake

  expect_equal(set_strategy("VoA", primary_uptake = 0.5),
  set_strategy("VoA", primary_uptake = 0.5, booster_uptake = 0.5))

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
  
test_that("create_vax_map and set_strategy working as expected", {

  # check errors generated when stop() if loops activated

  expect_error(set_strategy("VbE", primary_uptake = c(0, 0.5)))
  expect_error(set_strategy("VbE", primary_uptake = 1,
                            booster_uptake = c(0, 0.5)))

  # primary and booster vaccination can be different and map correctly

  primary_uptake <- 0.75
  booster_uptake <- 0.5

  v <- set_strategy("VoA", primary_uptake, booster_uptake)
  vax_map <- create_vax_map(n_vax = 5, v$vos, i_u = c(1, 3), i_v = c(2, 4))

  x_vaxmap <- vax_map[, 1, 1]
  w_vaxmap <- vax_map[, 3, 3]

  expect_equal(x_vaxmap[1], x_vaxmap[2])
  expect_equal(w_vaxmap[1], w_vaxmap[2])
  expect_true(x_vaxmap[1] != w_vaxmap[1])
  expect_true(x_vaxmap[1] == primary_uptake)
  expect_true(w_vaxmap[1] == booster_uptake)
  expect_equal(x_vaxmap, w_vaxmap * 1.5)

  # primary and booster vaccination can be different and map correctly for
  # strategies where only high activity groups receive vaccination

  v2 <- set_strategy("VoD(H)", primary_uptake, booster_uptake)
  v3 <- set_strategy("VoA(H)", primary_uptake, booster_uptake)
  v4 <- set_strategy("VoD(L)+VoA(H)", primary_uptake, booster_uptake)

  vax_map_v2_vod <- create_vax_map(n_vax = 5, v2$vod, i_u = c(1, 3),
                                   i_v = c(2, 4))
  vax_map_v3_vod <- create_vax_map(n_vax = 5, v3$vod, i_u = c(1, 3),
                                   i_v = c(2, 4))
  vax_map_v3_vos <- create_vax_map(n_vax = 5, v3$vos, i_u = c(1, 3),
                                   i_v = c(2, 4))
  vax_map_v4_vos <- create_vax_map(n_vax = 5, v4$vos, i_u = c(1, 3),
                                   i_v = c(2, 4))

  # low activity groups 0 according to strategy
  expect_equal(rowSums(vax_map_v2_vod[1, , c(1, 3)]), c(rep(0, 5)))
  expect_equal(rowSums(vax_map_v3_vod[1, , c(1, 3)]), c(rep(0, 5)))
  expect_equal(rowSums(vax_map_v3_vos[1, , c(1, 3)]), c(rep(0, 5)))
  expect_equal(rowSums(vax_map_v4_vos[1, , c(1, 3)]), c(rep(0, 5)))

  # high activity groups not 0 according to strategy, and match primary and
  # booster uptake values
  expect_equal(rowSums(vax_map_v2_vod[2, , c(1, 3)]),
               c(primary_uptake, -primary_uptake,
                 booster_uptake, -booster_uptake, 0))

  expect_equal(rowSums(vax_map_v3_vod[2, , c(1, 3)]),
               c(primary_uptake, -primary_uptake,
                 booster_uptake, -booster_uptake, 0))

  expect_equal(rowSums(vax_map_v3_vos[2, , c(1, 3)]),
               c(primary_uptake, -primary_uptake,
                 booster_uptake, -booster_uptake, 0))

  expect_equal(rowSums(vax_map_v4_vos[2, , c(1, 3)]),
               c(primary_uptake, -primary_uptake,
                 booster_uptake, -booster_uptake, 0))

  # check vbe working correctly

  vbe <- 1
  v <- matrix(c(rep(vbe, 2), rep(0, 2)), nrow = 2, byrow = TRUE)

  vax_map_vbe <- create_vax_map(n_vax = 5, v = v, i_u = c(1, 3),
                                i_v = c(2, 4))

  expect_equal(vax_map_vbe[, 1, 1], c(vbe, vbe))
  expect_equal(vax_map_vbe[, 3, 3], c(0, 0))

})
