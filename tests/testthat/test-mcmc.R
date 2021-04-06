context("mcmc")

test_that("compare function works as expected", {
  pars <- read_csv(gonovax_file("extdata/gono_params.csv"))[1, ]
  ll <- compare_basic(pars)
  expect_equal(ll, -16506.11, tolerance = 1e-6)
  data <- cache$data
  y <- run(seq(0, nrow(data)), gono_params(1)[[1]])
  attended <- diff(y$tot_attended)
  diagnosed <- diff(y$tot_treated)

  lldiagnosed <- Map(dnorm, x = data$diagnosed, mean = diagnosed,
                     sd = sqrt(diagnosed), log = TRUE)
  year <- gonovax_year(2019)
  llattended <- dnorm(data$attended, attended, sqrt(attended), TRUE)[year]

  expect_equal(ll, sum(unlist(c(lldiagnosed, llattended)), na.rm = TRUE),
               tolerance = 1e-6)

  ## plot comparison
  par(mfrow = c(1, 2), bty = "n", mar = c(3, 3, 1, 1), mgp = c(1.7, 0.5, 0))
  plot(data$year, data$attended, type = "l", ylim = c(0, 3e5))
  lines(data$year, attended, col = "darkred")
  plot(data$year, data$diagnosed, type = "l", ylim = c(0, 3e4))
  lines(data$year, diagnosed, col = "darkred")

})

test_that("mcmc runs", {

  z <- example_mcmc()
  expect_equal(z$mcmc$probabilities[, "log_likelihood"],
               c(-16506.1135588395, -16506.1135588395, -12302.0329021242,
                 -10386.6624356711, -10386.6624356711, -10386.6624356711,
                 -8292.725909926, -8292.725909926, -7002.84058594194,
                 -5408.29357213292, -5152.0304181118), tolerance = 1e-1)

  ## what do results look like?
  p <- z$mcmc$pars[nrow(z$mcmc$pars), ]
  mod <- model(user = model_params(transform0(p)),
               unused_user_action = "ignore")
  data <- cache$data
  y <- mod$run(seq(0, nrow(data)))

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

test_that("mcmc runs with multiple chains", {
  pars <- example_mcmc()$pars
  set.seed(1)
  z <- mcmc(pars, n_steps = 4, progress = TRUE, n_chains = 2)
  expect_equal(z$probabilities[, "log_likelihood"],
               c(-16506.1135588395, -16506.1135588395, -12302.0329021242,
                 -10386.6624356711, -10386.6624356711, -16506.1135588395,
                 -16506.1135588395, -9123.11679554351, -9123.11679554351,
                 -5985.5066767372), tolerance = 1e-1)

  ## no progress
  set.seed(1)
  z1 <- mcmc(pars, n_steps = 4, progress = FALSE, n_chains = 2)
  expect_equal(z, z1)

})
