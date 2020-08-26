
test_that("can select a specific parameter sets", {
  gp <- read_csv(gonovax_file("extdata/gono_params.csv"))
  # check that null argument returns all params
  expect_equal(gono_params(), gp)
  # check that can extract a single parameter
  expect_equal(gono_params(1), gp[1, ])
  # check that can extract multiple parameters
  expect_equal(gono_params(2:3), gp[2:3, ])
  # check that negative paramters will not be returned
  expect_equal(gono_params(c(-1, 1)), gp[1, ])
  # check that cannot extend beyond parameter set
  expect_equal(gono_params(c(1, nrow(gp) + 1)), gp[1, ])
})

test_that("vax_map works correctly", {
  # check error when input more than n_group vaccine uptake
  expect_error(create_vax_map(2, rep(0, 3)))
  # check error when uptake not in 0-1
  expect_error(create_vax_map(2, c(0, 2)))
  expect_error(create_vax_map(2, c(0, -1)))

  # indices of y [group, out, in]
  # test onevax map
  y <- create_vax_map(2, c(0.1, 0.2))
  expect_equal(y[1, 1, 1], 0.1)
  expect_equal(y[2, 1, 1], 0.2)
  expect_equal(sum(y[, , -1]), 0)
  expect_true(all(apply(y, c(1, 3), sum) == 0))

  # test onevax waning map
  y <- create_vax_map(5, c(0.1, 0.2))
  expect_equal(y[1, 1, 1], 0.1)
  expect_equal(y[2, 1, 1], 0.2)
  expect_equal(sum(y[, , -1]), 0)
  expect_true(all(apply(y, c(1, 3), sum) == 0))
})

test_that("waning map works correctly", {

  expect_error(create_waning_map(2, 1, -1))

  # test onevax map
  y <- create_waning_map(2, 1, 1)
  expect_equal(rowSums(y), c(1, -1))
  expect_equal(colSums(y), c(0, 0))

  # test onevax waning map
  y <- create_waning_map(5, 5, 1)
  expect_equal(rowSums(y), c(0, -1, 0, 0, 1))
  expect_equal(colSums(y), c(0, 0, 0, 0, 0))
})
