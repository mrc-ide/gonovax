
test_that("compare function works as expected", {
  pars <- gono_params(1)
  ll <- compare(pars)
  expect_equal(ll, -3833.782, tol = 1e-6)
  data <- cache$data
  y <- run_novax(1, seq(0, 12))[[1]]
  attended <- diff(y$tot_attended)
  diagnosed <- diff(y$tot_treated)

  lldiagnosed <- Map(dnorm, x = data$diagnosed, mean = diagnosed,
                     sd = sqrt(diagnosed), log = TRUE)
  year <- gonovax_year(2019)
  llattended <- dnorm(data$attended, attended, sqrt(attended), TRUE)[year]

  expect_equal(ll, sum(unlist(c(lldiagnosed, llattended)), na.rm = TRUE),
               tol = 1e-6)

  ## plot comparison
  par(mfrow = c(1, 2), bty = "n", mar = c(3, 3, 1, 1), mgp = c(1.7, 0.5, 0))
  plot(data$year, data$attended, type = "l", ylim = c(0, 3e5))
  lines(data$year, attended, col = "darkred")
  plot(data$year, data$diagnosed, type = "l", ylim = c(0, 3e4))
  lines(data$year, diagnosed, col = "darkred")

})

test_that("mcmc runs", {
  p <- unlist(gono_params(1))
  proposal <- cov(gono_params()) * 0.1
  pars <- Map(mcstate::pmcmc_parameter,
              name = names(p),
              initial = unname(p),
              min = 0,
              max = c(1, 1, 1, 1, 1, Inf, Inf, Inf, Inf, Inf),
              discrete = FALSE)
  pars <- mcstate::pmcmc_parameters$new(pars, proposal)
  set.seed(1)
  z <- mcmc(pars, n_steps = 10, progress = TRUE)
  expect_equal(z$probabilities[, "log_likelihood"],
               c(-3833.78202033004, -3833.78202033004, -3570.14662139024,
                 -3570.14662139024, -3570.14662139024, -3554.31011602949,
                 -3554.31011602949, -3554.31011602949, -3554.31011602949,
                 -3554.31011602949, -3543.37004716602), tol = 1e-1)

  ## what do results look like?
  p <- z$pars[nrow(z$pars), ]
  model_pars <- model_params(p)
  mod <- model(user = model_pars)
  y <- mod$run(seq(0, 12))
  data <- cache$data
  attended <- diff(y[, "tot_attended"])
  diagnosed <- diff(y[, "tot_treated"])

  par(mfrow = c(1, 2), bty = "n", mar = c(3, 3, 1, 1), mgp = c(1.7, 0.5, 0))
  plot(data$year, data$attended, type = "l", ylim = c(0, 3e5))
  lines(data$year, attended, col = "darkred")
  plot(data$year, data$diagnosed, type = "l", ylim = c(0, 3e4))
  lines(data$year, diagnosed, col = "darkred")

  par(mfrow = c(2, 5))
  lapply(colnames(z$pars),
         function(y) plot(z$iteration, z$pars[, y], type = "l",
                          ylab = y, xlab = "iter"))

})

test_that("mcmc_combine works as expected", {
  p <- unlist(gono_params(1))
  proposal <- cov(gono_params()) * 0.1
  params <- Map(mcstate::pmcmc_parameter,
              name = names(p),
              initial = unname(p),
              min = 0,
              max = c(1, 1, 1, 1, 1, Inf, Inf, Inf, Inf, Inf),
              discrete = FALSE)
  pars <- mcstate::pmcmc_parameters$new(params, proposal)
  set.seed(1)
  z1 <- mcmc(pars, n_steps = 5, progress = TRUE)
  z2 <- mcmc(pars, n_steps = 5, progress = TRUE)
  z <- mcmc_combine(z1, z2)

  # check error cases
  class(z2) <- NULL
  expect_error(mcmc_combine(z1, z2),
               "All elements of '...' must be 'gonovax_mcmc' objects")

  z3 <- mcmc(pars, n_steps = 4, progress = TRUE)
  expect_error(mcmc_combine(z1, z3), "All chains must have the same length")
  expect_error(mcmc_combine(z, z3), "Chains have already been combined")
  expect_error(mcmc_combine(z1), "At least 2 samples objects must be provided")

  pars <- mcstate::pmcmc_parameters$new(params[-1], proposal[-1, -1])
  z4 <-  mcmc(pars, n_steps = 5, progress = TRUE)
  expect_error(mcmc_combine(z1, z4), "All parameters must have the same names")

})

test_that("mcmc runs with multiple chains", {
  p <- unlist(gono_params(1))
  proposal <- cov(gono_params()) * 0.1
  pars <- Map(mcstate::pmcmc_parameter,
              name = names(p),
              initial = unname(p),
              min = 0,
              max = c(1, 1, 1, 1, 1, Inf, Inf, Inf, Inf, Inf),
              discrete = FALSE)
  pars <- mcstate::pmcmc_parameters$new(pars, proposal)
  set.seed(1)
  z <- mcmc(pars, n_steps = 4, progress = TRUE, n_chains = 2)
  expect_equal(z$probabilities[, "log_likelihood"],
               c(-3833.78202033004, -3256.29981357909, -3256.29981357909,
                 -3256.29981357909, -3256.29981357909, -3833.78202033004,
                 -3833.78202033004, -3363.93250303873, -3363.93250303873,
                 -3161.02224125879), tol = 1e-1)

  ## no progress
  set.seed(1)
  z1 <- mcmc(pars, n_steps = 4, progress = FALSE, n_chains = 2)
  expect_equal(z, z1)

})