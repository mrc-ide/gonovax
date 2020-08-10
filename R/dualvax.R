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
##' @param gono_params A vector of natural history parameters
##' @param vax_params A vector of vaccination params
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
##' @param pars A parameter list created by [basic_parameters()]; from
##'   this list we will use the `N0` `q`, `prev_Asl` and `prev_Ash`
##'   elements.
##'
##' @return A list of initial conditions
##' @export
##' 
dualvax_initial <- function(pars) {

  N0 <- round(pars$N0 * c(pars$ql, 1- pars$ql))
  A0 <- round(N0 * c(pars$prev_Asl, pars$prev_Ash))
  I0 <- S0 <- T0 <- c(0, 0)
  
  list(U0 = N0 - A0,
       I0 = I0,
       A0 = A0,
       S0 = S0,
       T0 = T0)
}



