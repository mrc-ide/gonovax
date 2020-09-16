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
mcmc <- function(pars, n_steps, progress = FALSE, n_chains = 1) {

  assert_is(pars, "pmcmc_parameters")
  assert_scalar_positive_integer(n_steps)
  assert_scalar_positive_integer(n_chains)

  if (n_chains == 1) {
    mcmc_single_chain(pars, n_steps, progress)
  } else {
    samples <- vector("list", n_chains)
    for (i in seq_along(samples)) {
      if (progress) {
        message(sprintf("Running chain %d / %d", i, n_chains))
      }
      samples[[i]] <- mcmc_single_chain(pars, n_steps, progress)
    }
  mcmc_combine(samples = samples)
  }
}

mcmc_single_chain <- function(pars, n_steps, progress) {

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

#' comparison function used to calculate log-likelihood of data given model
#' simulation - this may be replaced with a user input if we want to try out
#' different fits / models

compare <- function(pars) {

  ## load data
  if (is.null(cache$data)) {
    cache$data <- read_csv(gonovax_file("extdata/observed_gumcad.csv"))
  }

  data <- cache$data
  # initially only fit to most recent attendance figures
  data[data$year < 2019, "attended"] <- NA

  ## run odin model
  pars <- model_params(gono_params = pars)
  mod <- model(user = pars, unused_user_action = FALSE)
  y <- mod$run(seq(0, nrow(data)))

  ## output total treated and total attended per year from odin model
  diagnosed <- diff(y[, "tot_treated"])
  attended <- diff(y[, "tot_attended"])

  ## compare to data
  lltreated <- dnorm(data$diagnosed, diagnosed, sqrt(diagnosed), TRUE)
  llattended <- dnorm(data$attended, attended, sqrt(attended), TRUE)

  ## output log likelihood
  sum(lltreated, llattended, na.rm = TRUE)
}

mcmc_progress <- mcstate:::pmcmc_progress

gonovax_mcmc <- function(pars, probabilities, chain = NULL,
                         iteration = NULL) {

    ret <- list(chain = chain,
              iteration = iteration %||% seq.int(0, length.out = nrow(pars)),
              pars = pars,
              probabilities = probabilities)
  class(ret) <- "gonovax_mcmc"
  ret
}

##' Combine multiple [mcmc()] samples into one object
##'
##' @title Combine mcmc samples
##'
##' @param ... Arguments representing [mcmc()] sample, i.e.,
##'   `gonovax_mcmc` objects. Alternatively, pass a list as the
##'   argument `samples`. Names are ignored.
##'
##' @param samples A list of `gonovax_mcmc` objects. This is often
##'   more convenient for programming against than `...`
##'
##' @export
mcmc_combine <- function(..., samples = list(...)) {

  if (!all(vlapply(samples, inherits, "gonovax_mcmc"))) {
    stop("All elements of '...' must be 'gonovax_mcmc' objects")
  }
  if (any(!vlapply(samples, function(x) is.null(x$chain)))) {
    stop("Chains have already been combined")
  }
  if (length(samples) < 2) {
    stop("At least 2 samples objects must be provided")
  }
  if (length(unique(lapply(samples, function(x) colnames(x$pars)))) != 1L) {
    stop("All parameters must have the same names")
  }
  if (length(unique(vnapply(samples, function(x) nrow(x$pars)))) != 1L) {
    stop("All chains must have the same length")
  }
  iteration <- lapply(samples, "[[", "iteration")
  if (length(unique(iteration)) != 1L) {
    stop("All chains must have the same iterations")
  }

  iteration <- unlist(iteration, FALSE, FALSE)

  chain <- rep(seq_along(samples), each = nrow(samples[[1]]$pars))

  pars <- do.call(rbind, lapply(samples, "[[", "pars"))
  probabilities <- do.call(rbind, lapply(samples, "[[", "probabilities"))

  gonovax_mcmc(pars, probabilities, chain, iteration)
}
