##' @name dualvax
##' @title Model of gonorrhoea with dual vaccines
##'This is an odin model.
##' @export dualvax
NULL

##' @name dualvax_params
##' @title Parameters for the dualvax model
##' @param gono_params A dataframe of natural history parameters
##' @param vax_params A vector of vaccination params
##' @param init_params A list of starting conditions
##' @return A list of inputs to the model many of which are fixed and
##'   represent data. These correspond largely to `user()` calls
##'   within the odin code, though some are also used in processing
##'   just before the model is run.
##' @export
dualvax_params <- function(gono_params = NULL,
                               init_params = NULL,
                               vax_params = NULL) {
  ret <- c(demographic_params(), gono_params)
  init_params <- init_params %||% dualvax_initial(ret)
  vax_params <- vax_params %||% vax_params0()
  c(ret, init_params, vax_params)
}

##' Create initial conditions for the dualvax model
##'
##' @title Initial conditions for the dualvax model
##'
##' @param pars A parameter list created by [dualvax_params()]; from
##'   this list we will use the `N0`, `q`, `prev_Asl` and `prev_Ash`
##'   elements.
##' @return A list of initial conditions
##' @export
dualvax_initial <- function(pars) {
  # separate into 1:low and 2:high activity groups
  N0 <- round(pars$N0 * pars$q)
  U0 <- I0 <- A0 <- S0 <- T0 <- array(0, c(length(pars$beta), 2, 1))

  # set initial asymptomatic prevalence in each group
  A0[, 1, 1] <- round(N0[1] * pars$prev_Asl)
  A0[, 2, 1] <- round(N0[2] * pars$prev_Ash)

  # set initial uninfecteds
  U0[, 1, 1] <- N0[1] - A0[, 1, 1]
  U0[, 2, 1] <- N0[2] - A0[, 2, 1]

  list(U0 = U0, I0 = I0, A0 = A0, S0 = S0, T0 = T0)
}

vax_params0 <- function(x) {
  list(ve = 1,
       eff = 0,
       w = as.matrix(0))
}

vax_params1 <- function(x) {
  z <- 1 / 4
  list(ve = 1,
       eff = 0,
       w = rbind(c(0, z),
                c(-z, 0))
  )
}
