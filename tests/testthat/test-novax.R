test_that("run_novax works correctly", {
  tt <- c(0, 1)
  y1 <- run_novax_int(1, tt)
  y2 <- run_novax_int(2, tt)

  y1e <- run_novax_int(1, tt, equilib = TRUE)
  y2e <- run_novax_int(2, tt, equilib = TRUE)

  expect_error(run_novax_int(1:2, tt, equilib = TRUE),
               label = "if length(n) > 1, equilib must be FALSE")

  y1m <- run_novax(1, tt)
  y12m <- run_novax(1:2, tt)
  y1me <- run_novax(1, tt, equilib = TRUE)
  y12me <- run_novax(1:2, tt, equilib = TRUE)

  expect_equal(y1, y1m[[1]])
  expect_equal(y1e, y1me[[1]])
  expect_equal(y1, y12m[[1]])
  expect_equal(y2, y12m[[2]])
  expect_equal(y1e, y12me[[1]])
  expect_equal(y2e, y12me[[2]])

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

  y2 <- novax_baseline(1:100, t = 10)
  expect_equal(length(y2[[1]]), 100)

  expect_error(novax_baseline(1, 11), "t must be an integer between 1 and 10")
})
