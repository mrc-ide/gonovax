
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

test_that("novax_equilib returns equlibrium conditions", {

  y <- novax_equilib(1)
  tt <- y$t
  params <- dualvax_params(gono_params = gono_params(1))
  mod <- dualvax(user = params)

  y1 <- mod$run(tt)
  y1 <- mod$transform_variables(y1)

  expect_equal(y, y1)
# check that can load all parameters at once
  y2 <- novax_equilib()
  expect_equal(length(y2), 982)
  expect_equal(y2[[1]], y1)

})
