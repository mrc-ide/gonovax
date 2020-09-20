##' all functions in mcmc-tools.R are adapted from pmcmc-tools.R functions in
##' github.com/mrc-ide/mcstate

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
  iteration <- unlist(iteration, FALSE, FALSE)

  chain <- rep(seq_along(samples), each = nrow(samples[[1]]$pars))

  pars <- do.call(rbind, lapply(samples, "[[", "pars"))
  probabilities <- do.call(rbind, lapply(samples, "[[", "probabilities"))

  gonovax_mcmc(pars, probabilities, chain, iteration)
}



##' Thin results of running [mcmc()].`mcmc_thin` takes every `thin`'th sample,
##'  while `mcmc_sample` randomly selects a total of `n_sample` samples.
##'
##' @title Thin a mcmc chain
##' @param object Results of running [mcmc()]
##'
##' @param burnin Optional integer number of iterations to discard as
##'   "burn-in". If given then samples `1:burnin` will be
##'   excluded from your results. Remember that the first sample
##'   represents the starting point of the chain. It is an error if
##'   this is not a positive integer or is greater than or equal to
##'   the number of samples (i.e., there must be at least one sample
##'   remaining after discarding burnin).
##'
##' @param thin Optional integer thinning factor. If given, then every
##'   `thin`'th sample is retained (e.g., if `thin` is 10
##'   then we keep samples 1, 11, 21, ...).
##'
##' @export
mcmc_thin <- function(object, burnin = NULL, thin = NULL) {
  assert_is(object, "gonovax_mcmc")
  i <- rep_len(TRUE, length(object$iteration))

  if (!is.null(burnin)) {
    assert_scalar_positive_integer(burnin)
    burnin_max <- max(object$iteration)
    if (burnin >= burnin_max) {
      stop(sprintf("'burnin' must be less than %d for your results",
                   burnin_max))
    }
    i <- i & object$iteration >= burnin
  }

  if (!is.null(thin)) {
    assert_scalar_positive_integer(thin)
    offset <- min(object$iteration[i])
    i <- i & ((object$iteration - offset) %% thin == 0)
  }

  mcmc_filter(object, i)
}


##' @param n_sample The number of samples to draw from `object` *with
##'   replacement*. This means that `n_sample` can be larger than the
##'   total number of samples taken (though it probably should not)
##'
##' @export
##' @rdname mcmc_thin
mcmc_sample <- function(object, n_sample, burnin = NULL) {
  object <- mcmc_thin(object, burnin)
  i <- sample.int(nrow(object$pars), n_sample, replace = TRUE)
  mcmc_filter(object, i)
}


mcmc_filter <- function(object, i) {
  if (!is.null(object$chain)) {
    object$chain <- object$chain[i]
  }
  object$iteration <- object$iteration[i]
  object$pars <- object$pars[i, , drop = FALSE]
  object$probabilities <- object$probabilities[i, , drop = FALSE]

  object
}
