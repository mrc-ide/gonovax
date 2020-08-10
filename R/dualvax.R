##' @name dualvax
##' @title Model of gonorrhoea with dual vaccines
##'This is an odin model.
##' @export dualvax
##' 
NULL

##' @name Parameters for the "[dualvax]" model.
##' 
##' @title Parameters for the dulavax model
##' 
##' @param gono_params A dataframe of natural history parameters
##' @param vax_params A vector of vaccination params
##' @param init_params A list of starting conditions
##' @return A list of inputs to the model many of which are fixed and
##'   represent data. These correspond largely to `user()` calls
##'   within the odin code, though some are also used in processing
##'   just before the model is run.
##' @export
dualvax_parameters <- function(gono_params = NULL,
                               init_params = NULL,
                               vax_params = NULL) {
  ret <- gonovax_parameters()
  inits <- init_params %||% dualvax_initial(ret)
  c(ret, inits)
}

##' Create initial conditions for the dualvax model
##'
##' @title Initial conditions for the dualvax model
##'
##' @param pars A parameter list created by [dual_parameters()]; from
##'   this list we will use the `N0` `q`, `prev_Asl` and `prev_Ash`
##'   elements.
##'
##' @return A list of initial conditions
##' @export
##' 
dualvax_initial <- function(pars) {
  n_par <- length(pars$beta)
  N0 <- matrix(pars$N0 * c(pars$ql, 1 - pars$ql), byrow = TRUE, 
               nrow = n_par, ncol = 2 )
  A0 <- round(cbind(N0[, 1] * pars$prev_Asl, N0[, 2] * pars$prev_Ash))

  I0 <- S0 <- T0 <- matrix(0, nrow = nrow(A0), ncol = ncol(A0))

  list(U0 = N0 - A0,
       I0 = I0,
       A0 = A0,
       S0 = S0,
       T0 = T0)
}



