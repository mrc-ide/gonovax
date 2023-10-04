context("mcmc (tools)")

test_that("mcmc_thin with no args is a no-op", {
  results <- example_mcmc()$mcmc
  expect_identical(mcmc_thin(results), results)
})


test_that("discarding burnin drops beginnings of chain", {
  results <- example_mcmc()$mcmc
  res <- mcmc_thin(results, 5)
  i <- 6:11
  expect_identical(res$pars, results$pars[i, ])
  expect_identical(res$probabilities, results$probabilities[i, ])
})


test_that("thinning drops all over chain", {
  results <- example_mcmc()$mcmc
  res <- mcmc_thin(results, thin = 2)
  i <- seq(1, 11, by = 2)
  expect_identical(res$pars, results$pars[i, ])
  expect_identical(res$probabilities, results$probabilities[i, ])
})


test_that("burnin and thin can be used together", {
  results <- example_mcmc()$mcmc
  i <- seq(3, 11, by = 2)
  res <- mcmc_thin(results, 2, 2)
  expect_identical(res$pars, results$pars[i, ])
  expect_identical(res$probabilities, results$probabilities[i, ])
})


test_that("can't discard the whole chain (or more)", {
  results <- example_mcmc()$mcmc
  expect_error(mcmc_thin(results, 11),
               "'burnin' must be less than 10 for your results")
  expect_error(mcmc_thin(results, 100),
               "'burnin' must be less than 10 for your results")
})


test_that("Can thin when no state/trajectories present", {
  results <- example_mcmc()$mcmc

  i <- seq(3, 11, by = 2)
  res <- mcmc_thin(results, 2, 2)
  expect_identical(res$pars, results$pars[i, ])
  expect_identical(res$probabilities, results$probabilities[i, ])
})


test_that("can combine chains", {
  results <- example_mcmc2()$results

  results1 <- results[[1]]
  results2 <- results[[2]]
  results3 <- results[[3]]

  res <- mcmc_combine(results1, results2, results3)

  n_mcmc <- nrow(results1$pars)
  n_par <- ncol(results1$pars)

  n_mcmc3 <- n_mcmc * 3

  expect_equal(dim(res$pars), c(n_mcmc3, n_par))
  expect_equal(dim(res$probabilities), c(n_mcmc3, 3))

  i <- seq_len(n_mcmc) + n_mcmc
  expect_equal(res$pars[i, ], results2$pars)
  expect_equal(res$probabilities[i, ], results2$probabilities)
})


test_that("can combine chains with list interface", {
  results <- example_mcmc2()$results
  expect_identical(
                   mcmc_combine(results[[1]], results[[2]], results[[3]]),
                   mcmc_combine(samples = results))
})


test_that("can drop burnin from combined chains", {
  results <- example_mcmc2()$results
  combined <- mcmc_combine(samples = results)
  res <- mcmc_thin(combined, burnin = 5)
  expect_equal(res$chain, rep(1:3, each = 6))
  expect_equal(res$iteration, rep(5:10, 3))

  ## Same performed either way:
  expect_identical(
                   res,
                   mcmc_combine(samples = lapply(results, mcmc_thin,
                                                 burnin = 5)))
})


test_that("can thin combined chains", {
  results <- example_mcmc2()$results
  combined <- mcmc_combine(samples = results)
  res <- mcmc_thin(combined, burnin = 2, thin = 2)
  expect_equal(res$chain, rep(1:3, each = 5))
  expect_equal(res$iteration, rep(seq(2, 10, by = 2), 3))

  ## Same performed either way:
  expect_identical(
                   res,
                   mcmc_combine(samples = lapply(results, mcmc_thin, 2, 2)))
})

test_that("mcmc_combine works as expected", {
  z1 <- example_mcmc()

  set.seed(2)
  z2 <- mcmc(z1$pars, n_steps = 10, progress = TRUE)
  z <- mcmc_combine(z1$mcmc, z2)

  # check error cases
  class(z2) <- NULL
  expect_error(mcmc_combine(z1$mcmc, z2),
               "All elements of '...' must be 'gonovax_mcmc' objects")

  z3 <- mcmc(z1$pars, n_steps = 4, progress = TRUE)
  expect_error(mcmc_combine(z1$mcmc, z3),
               "All chains must have the same length")
  expect_error(mcmc_combine(z, z3), "Chains have already been combined")
  expect_error(mcmc_combine(z3), "At least 2 samples objects must be provided")
  z4 <- z3
  colnames(z4$pars) <- paste0("par", seq_len(ncol(z4$pars)))
  expect_error(mcmc_combine(z3, z4), "All parameters must have the same names")
})



test_that("can sample from a mcmc", {
  results <- example_mcmc()$mcmc
  sub <- mcmc_sample(results, 8, burnin = 2)
  expect_equal(nrow(sub$pars), 8)
  expect_true(all(sub$iteration >= 2))
})


test_that("sampling is with replacement", {
  results <- example_mcmc()$mcmc
  sub <- mcmc_sample(results, 50, burnin = 2)
  expect_equal(nrow(sub$pars), 50)
  expect_true(all(sub$iteration >= 2))
  expect_true(any(duplicated(sub$iteration)))
})


test_that("can sample from a combined chain", {
  results <- mcmc_combine(samples = example_mcmc2()$results)
  sub <- mcmc_sample(results, 50, burnin = 2)
  expect_equal(nrow(sub$pars), 50)
  expect_true(all(1:3 %in% sub$chain))
  expect_true(all(sub$iteration >= 2))
})
