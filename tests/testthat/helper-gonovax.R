test_cache <- new.env()

example_mcmc <- function() {
  if (is.null(test_cache$example_mcmc)) {
    p <- read_csv(gonovax_file("extdata/gono_params.csv"))
    proposal <- cov(p) * 0.1
    params <- Map(mcstate::pmcmc_parameter,
                  name = names(p),
                  initial = unname(p[1, ]),
                  min = 0,
                  max = c(1, 1, 1, 1, 1, Inf, Inf, Inf, Inf, Inf),
                  discrete = FALSE)
    z <- list()
    z$pars <- mcstate::pmcmc_parameters$new(params, proposal)
    set.seed(1)
    z$mcmc <- mcmc(z$pars, n_steps = 10, progress = TRUE)
    test_cache$example_mcmc <- z
  }
  test_cache$example_mcmc
}

example_mcmc2 <- function() {
  if (is.null(test_cache$example_mcmc2)) {
    p <- read_csv(gonovax_file("extdata/gono_params.csv"))
    proposal <- cov(p) * 0.1
    params <- Map(mcstate::pmcmc_parameter,
                  name = names(p),
                  initial = unname(p[2, ]),
                  min = 0,
                  max = c(1, 1, 1, 1, 1, Inf, Inf, Inf, Inf, Inf),
                  discrete = FALSE)
    z <- list()
    z$pars <- mcstate::pmcmc_parameters$new(params, proposal)
    set.seed(1)
    z$results <- list(
      mcmc(z$pars, 10),
      mcmc(z$pars, 10),
      mcmc(z$pars, 10))
    test_cache$example_mcmc2 <- z
  }
  test_cache$example_mcmc2
}
