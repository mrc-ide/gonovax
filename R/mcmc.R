##' Run a mcmc sampler
##'
##' This is a basic Metropolis-Hastings MCMC sampler.  The
##' `model` is run with a set of parameters to evaluate the
##' likelihood. A new set of parameters is proposed, and these
##' likelihoods are compared, jumping with probability equal to their
##' ratio. This is repeated for `n_steps` proposals.
##'
##' This function is adapted from the `pmcmc` in the mcstate package
##' https://github.com/mrc-ide/mcstate
##'
##' @title Run a mcmc sampler
##'
##' @param pars A [`mcmc_parameters`] object containing
##'   information about parameters ( parameter ranges, priors, proposal kernel,
##'   observation functions).
##'
##' @param n_steps Number of MCMC steps to run
##'
##' @param compare likelihood function to compare data to epidemic trajectory
##' should return a single value representing the log-likelihood. Default is
##' compare_basic()
##'
##' @param progress Logical, indicating if a progress bar should be
##'   displayed, using [`progress::progress_bar`].
##'
##' @param n_chains Optional integer, indicating the number of chains
##'   to run. If more than one then we run a series of chains and
##'   merge them with [mcmc_combine()]. Chains are run in series,
##'   with the same model.
##'
##' @return A `gonovax_mcmc` object containing `pars`
##'   (sampled parameters) and `probabilities` (log prior, log
##'   likelihood and log posterior values for these
##'   probabilities).
##'
##' @export
##' @importFrom stats runif dnorm
mcmc <- function(pars, n_steps, compare = NULL, progress = FALSE,
                 n_chains = 1) {

  assert_is(pars, "pmcmc_parameters")
  assert_scalar_positive_integer(n_steps)
  assert_scalar_positive_integer(n_chains)

  compare <- compare %||% compare_basic

  if (n_chains == 1) {
    mcmc_single_chain(pars, n_steps, compare, progress)
  } else {
    samples <- vector("list", n_chains)
    for (i in seq_along(samples)) {
      if (progress) {
        message(sprintf("Running chain %d / %d", i, n_chains))
      }
      samples[[i]] <- mcmc_single_chain(pars, n_steps, compare, progress)
    }
    mcmc_combine(samples = samples)
  }
}

mcmc_single_chain <- function(pars, n_steps, compare, progress) {

  history_pars <- history_collector(n_steps)
  history_probabilities <- history_collector(n_steps)

  curr_pars <- pars$initial()
  curr_lprior <- pars$prior(curr_pars)
  curr_llik <- compare(curr_pars)
  curr_lpost <- curr_lprior + curr_llik

  history_pars$add(curr_pars)
  history_probabilities$add(c(curr_lprior, curr_llik, curr_lpost))

  tick <- mcmc_progress(n_steps, progress)

  for (i in seq_len(n_steps)) {
    tick()

    prop_pars <- pars$propose(curr_pars)
    prop_lprior <- pars$prior(prop_pars)
    prop_llik <- compare(prop_pars)
    prop_lpost <- prop_lprior + prop_llik

    if (runif(1) < exp(prop_lpost - curr_lpost)) {
      curr_pars <- prop_pars
      curr_lprior <- prop_lprior
      curr_llik <- prop_llik
      curr_lpost <- prop_lpost
    }

    history_pars$add(curr_pars)
    history_probabilities$add(c(curr_lprior, curr_llik, curr_lpost))

  }

  pars_matrix <- do.call(rbind, history_pars$get())
  probabilities <- do.call(rbind, history_probabilities$get())
  colnames(probabilities) <- c("log_prior", "log_likelihood", "log_posterior")

  gonovax_mcmc(pars_matrix, probabilities)
}


## Generic history collector, collects anything at all into a list
##
## This would be more nicely done as a simple R6 class but it's a bit
## slow in testing; this version speeds up the total mcmc runtime by a
## factor of ~3x (0.4s/1000 iterations to 0.13s/1000) mostly by
## reducing the number of garbage collections considerably.
history_collector <- function(n) {
  data <- vector("list", n + 1L)
  i <- 0L
  add <- function(value) {
    i <<- i + 1L
    data[[i]] <<- value
  }

  get <- function() {
    data
  }

  list(add = add, get = get)
}


##' @title Calculate the log likelihood of the data given the parameters
##'
##' @param pars A named vector of parameters
##'
##' @return a single log likelihood
compare_basic <- function(pars) {

  ## load data
  data <- gonovax_data()

  # initially only fit to most recent attendance figures
  data[data$year < 2019, "attended"] <- NA

  ## run odin model
  y <- run(tt = c(0, data$gonovax_year),
           gono_params =  transform0(pars),
           transform = FALSE)

  ## output total treated and total attended per year from odin model
  diagnosed <- diff(y[, "tot_treated"])
  attended <- diff(y[, "tot_attended"])

  ## compare to data
  lltreated <- dnorm(data$diagnosed, diagnosed, sqrt(diagnosed), TRUE)
  llattended <- dnorm(data$attended, attended, sqrt(attended), TRUE)

  ## output log likelihood
  sum(lltreated, llattended, na.rm = TRUE)
}



##' @title Calculate the log likelihood of the data given the parameters
##' diagnoses and attendances lhoods are negative binomial
##' p_symp lhood is betabinomial
##' @param pars A named vector of parameters
##' @param transform the transform function to use in the comparison
##'
##' @return a single log likelihood
##' @importFrom stats dnbinom
##' @export

compare <- function(pars, transform) {

  ## load data
  data <- gonovax_data()
  w <- !is.na(data$n_reported)

  ## run odin model
  y <- run(tt = c(0, data$gonovax_year), gono_params =  transform(pars),
           transform = FALSE)

  ## output total treated and total attended per year from odin model

  diagnosed <- diff(y[, "tot_treated"])
  attended <- diff(y[, "tot_attended"])

  # indices of columns
  state_names <- colnames(y)
  i_diag_a <- grep("^cum_diag_a", state_names)
  i_diag_s <- grep("^cum_diag_s", state_names)

  diag_a <- pmax(diff(rowSums(y[, i_diag_a])), 1e-6)
  diag_s <- pmax(diff(rowSums(y[, i_diag_s])), 1e-6)
  p_symp <- diag_s / (diag_a + diag_s)

  # ensure p_symp is in [0, 1]
  stopifnot(max(p_symp) <= 1)
  stopifnot(min(p_symp) >= 0)

  ## compare to data
  ## in dnbinom size = shape param of gamma dist, mu = mean

  size <- 1 / pars[["k_gumcad"]]
  lltreated <- dnbinom(data$diagnosed, size = size, mu = diagnosed, log = TRUE)
  llattended <- dnbinom(data$attended, size = size, mu = attended, log = TRUE)

  ## in debetabinom size and prob are binomial parameters, rho is overdispersion
  ## with support [0, 1], 0 being a binomial
  llsymptomatic <- dbetabinom(x = data$n_symptomatic[w],
                              size = data$n_reported[w], prob = p_symp[w],
                              rho = pars[["k_grasp"]], log = TRUE)

  ## output log likelihood
  sum(lltreated, llattended, llsymptomatic, na.rm = TRUE)
}


# this is imported from mcstate

mcmc_progress <- function(n, show, force = FALSE) {
  if (show) {
    fmt <- "Step :current / :total [:bar] ETA :eta | :elapsedfull so far"
    t0 <- Sys.time()
    callback <- function(p) {
      message(sprintf("Finished %d steps in %s",
                      n, format(Sys.time() - t0, digits = 1)))
    }
    p <- progress::progress_bar$new(fmt, n, callback = callback, force = force)
    p$tick(0)
    p$tick
  } else {
    function() NULL
  }
}

gonovax_mcmc <- function(pars, probabilities, chain = NULL,
                         iteration = NULL) {

  ret <- list(chain = chain,
              iteration = iteration %||% seq.int(0, length.out = nrow(pars)),
              pars = pars,
              probabilities = probabilities)
  class(ret) <- "gonovax_mcmc"
  ret
}
