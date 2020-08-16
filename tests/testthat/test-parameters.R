
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
