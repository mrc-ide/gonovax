
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
  params <- model_params(gono_params = gono_params(1))
  mod <- model(user = params)

  y1 <- mod$run(tt)
  y1 <- mod$transform_variables(y1)

  expect_equal(y, y1)
# check that can load all parameters at once
  y2 <- novax_equilib()
  expect_equal(length(y2), 982)
  expect_equal(y2[[1]], y1)

})

test_that("novax_baseline returns baseline conditions", {

  y <- novax_baseline(1:3, 1:10)
  expect_equal(names(y), c("incid", "cum_incid"))
  expect_equal(dim(y$incid), c(10, 3))
  expect_equal(dim(y$cum_incid), c(10, 3))
  expect_equal(colSums(y$incid), y$cum_incid[10, ])

  expect_error(novax_baseline(1, 11), "t must be an integer between 1 and 10")
})
