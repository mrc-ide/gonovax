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
##' @param model A gonovax model object
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
mcmc <- function(pars, model, n_steps, progress = FALSE,
                  n_chains = 1) {
  assert_is(pars, "mcmc_parameters")
  assert_is(model, "gonovax_model")
  assert_scalar_positive_integer(n_steps)
  assert_scalar_positive_integer(n_chains)
  
  if (n_chains == 1) {
    mcmc_single_chain(pars, model, n_steps, progress)
  } else {
    samples <- vector("list", n_chains)
    for (i in seq_along(samples)) {
      if (progress) {
        message(sprintf("Running chain %d / %d", i, n_chains))
      }
      samples[[i]] <- mcmc_single_chain(pars, model, n_steps, progress)
    }
    mcmc_combine(samples = samples)
  }
}

mcmc_single_chain <- function(pars, model, n_steps,
                                
                               progress) {
  
  history_pars <- history_collector(n_steps)
  history_probabilities <- history_collector(n_steps)
  history_state <- history_collector(n_steps)
  history_trajectories <- history_collector(n_steps)
  
  curr_pars <- pars$initial()
  curr_lprior <- pars$prior(curr_pars)
  curr_llik <- model$run(pars$model(curr_pars))
  curr_lpost <- curr_lprior + curr_llik
  
  history_pars$add(curr_pars)
  history_probabilities$add(c(curr_lprior, curr_llik, curr_lpost))

  tick <- mcmc_progress(n_steps, progress)
  
  for (i in seq_len(n_steps)) {
    tick()
    
    prop_pars <- pars$propose(curr_pars)
    prop_lprior <- pars$prior(prop_pars)
    prop_llik <- model$run(pars$model(prop_pars))
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
  
  pars_matrix <- set_colnames(list_to_matrix(history_pars$get()),
                              names(curr_pars))
  probabilities <- set_colnames(list_to_matrix(history_probabilities$get()),
                                c("log_prior", "log_likelihood", "log_posterior"))
  
  predict <- state <- trajectories <- NULL
  
  
  mcstate_mcmc(pars_matrix, probabilities, state, trajectories, predict)
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
