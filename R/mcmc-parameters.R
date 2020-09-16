##' Describe a single parameter for use within the mcmc. Note that
##' the name is not set here, but will end up being naturally defined
##' when used with [`mcmc_parameters`], which collects
##' these together for use with [mcmc()]. This function is adapted from the
##' `pmcmc_parameters` in the mcstate package https://github.com/mrc-ide/mcstate
##'
##' @title Describe single mcmc parameter
##'
##' @param name Name for the parameter (a string)
##'
##' @param initial Initial value for the parameter
##'
##' @param min Optional minimum value for the parameter (otherwise
##'   `-Inf`). If given, then `initial` must be at least this
##'   value.
##'
##' @param max Optional max value for the parameter (otherwise
##'   `Inf`). If given, then `initial` must be at most this
##'   value.
##'
##' @param discrete Logical, indicating if this parameter is
##'   discrete. If `TRUE` then the parameter will be rounded
##'   after a new parameter is proposed.
##'
##' @param prior A prior function (if not given an improper flat prior
##'   is used - be careful!). It must be a function that takes a
##'   single argument, being the value of this parameter. If given,
##'   then `prior(initial)` must evaluate to a finite value.
##'
##' @export
##' @examples
##' mcmc_parameter("a", 0.1)
mcmc_parameter <- mcstate::pmcmc_parameter


##' @title mcmc_parameters
##'
##' @description Construct parameters for use with
##'   [mcmc()]. This creates a utility object that is used
##'   internally to work with parameters. Most users only need to
##'   construct this object, but see the examples for how it can be
##'   used.
##'
##' @export
mcmc_parameters <- mcstate::pmcmc_parameters

## create function to reflect proposal boundaries at pars_min and pars_max
## this ensures the proposal is symmetrical and we can simplify the MH step
reflect_proposal <- mcstate:::reflect_proposal

reflect_proposal_both <- mcstate:::reflect_proposal_both

reflect_proposal_one <- mcstate:::reflect_proposal_one
