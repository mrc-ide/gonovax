context("parameters")

test_that("can select specific parameter sets", {
  p <- read_csv(gonovax_file("extdata/gono_params.csv"))
  gp <- lapply(seq_len(nrow(p)), function(i) transform0(p[i, ]))
  # check that null argument returns all params
  expect_equal(gono_params(), gp)
  # check that can extract a single parameter
  expect_equal(gono_params(1)[[1]], transform0(p[1, ]))
  # check that can extract multiple parameters
  expect_equal(gono_params(2:3), gp[2:3])
  # check that negative paramters will not be returned
  expect_equal(gono_params(c(-1, 1)), gp[1])
  # check that cannot extend beyond parameter set
  expect_equal(gono_params(c(1, nrow(p) + 1)), gp[1])
})

test_that("vax_map works correctly", {
  # check error when input more than n_group vaccine uptake
  expect_error(create_vax_map(n_vax = 2, v = rep(0, 3), i_u = 1, i_v = 2))
  # check error when uptake not in 0-1
  expect_error(create_vax_map(n_vax = 2, v = c(0, 2), i_u = 1, i_v = 2))
  expect_error(create_vax_map(n_vax = 2, v = c(0, -1), i_u = 1, i_v = 2))
  # check error length(i_u) != length(i_v)
  expect_error(create_vax_map(n_vax = 3, v = c(0.1, 0.2), i_u = c(1, 3),
                              i_v = 2))
  # check error length(i_u) != length(i_v)
  expect_error(create_vax_map(n_vax = 2, v = c(0.1, 0.2), i_u = 1, i_v = 3))

  # indices of y [group, -to, -from]
  # test onevax map
  y <- create_vax_map(n_vax = 2, v = c(0.1, 0.2), i_u = 1, i_v = 2)
  expect_equal(y[1, , 1], c(0.1, -0.1))
  expect_equal(y[2, , 1], c(0.2, -0.2))
  expect_true(all(y[, , 2] == 0))

  # test onevax_xvwv waning map
  y <- create_vax_map(n_vax = 3, v = c(0.1, 0.2), i_u = c(1, 3), i_v = c(2, 2))
  expect_equal(y[1, , 1], c(0.1, -0.1, 0))
  expect_equal(y[2, , 1], c(0.2, -0.2, 0))
  expect_equal(y[1, , 3], c(0, -0.1, 0.1))
  expect_equal(y[2, , 3], c(0, -0.2, 0.2))
  expect_true(all(y[, , 2] == 0))

  # test onevax_xvwr waning map
  y <- create_vax_map(n_vax = 4, v = c(0.1, 0.2), i_u = c(1, 3), i_v = c(2, 4))
  expect_equal(y[1, , 1], c(0.1, -0.1, 0, 0))
  expect_equal(y[2, , 1], c(0.2, -0.2, 0, 0))
  expect_equal(y[1, , 3], c(0, 0, 0.1, -0.1))
  expect_equal(y[2, , 3], c(0, 0, 0.2, -0.2))
  expect_true(all(y[, , 2] == 0))
  expect_true(all(y[, , 4] == 0))

  # test onevax_xvwr waning map with single v input
  y <- create_vax_map(n_vax = 4, v = 0.1, i_u = c(1, 3), i_v = c(2, 4))
  expect_equal(y[1, , 1], c(0.1, -0.1, 0, 0))
  expect_equal(y[2, , 1], c(0.1, -0.1, 0, 0))
  expect_equal(y[1, , 3], c(0, 0, 0.1, -0.1))
  expect_equal(y[2, , 3], c(0, 0, 0.1, -0.1))
  expect_true(all(y[, , 2] == 0))
  expect_true(all(y[, , 4] == 0))

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
  p$gamma_l <- 0.6
  gp <- transform(p)

  # check ratio is fixed between eta_l and eta_h
  expect_equal(gp$eta_l_t / gp$eta_h_t, rep(p$gamma_l, length(gp$tt)))

  i_2020 <- which(gp$tt == gonovax_year(2020))

  ## check increasing to 2020
  expect_equal(diff(gp$beta_t[seq_len(i_2020)]),
               rep(p$beta2009 * p$phi_beta, i_2020 - 1L))
  expect_equal(diff(gp$eta_h_t[seq_len(i_2020)]),
               rep(p$eta_h * p$phi_eta, i_2020 - 1L))
  expect_equal(diff(gp$eta_l_t[seq_len(i_2020)]),
               rep(p$eta_l * p$phi_eta * p$gamma_l, i_2020 - 1L))

  ## check stable after 2020
  expect_true(all(diff(gp$beta_t[-seq_len(i_2020)]) == 0))
  expect_true(all(diff(gp$eta_h_t[-seq_len(i_2020)]) == 0))
  expect_true(all(diff(gp$eta_l_t[-seq_len(i_2020)]) == 0))

  #Check for fix_par_t = FALSE
  gp2 <- transform(p, FALSE)

  expect_equal(diff(gp2$beta_t), rep(p$beta2009 * p$phi_beta,
                                   length(gp2$tt) - 1L))
  expect_equal(diff(gp2$eta_h_t), rep(p$eta_h * p$phi_eta, length(gp2$tt) - 1L))
  expect_equal(diff(gp2$eta_l_t), rep(p$eta_h * p$phi_eta * p$gamma_l,
                                      length(gp2$tt) - 1L))

})
