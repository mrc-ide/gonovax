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
  vax_params <- vax_params %||% vax_params0()
  ret <- c(demographic_params(), gono_params, vax_params)
  init_params <- init_params %||% dualvax_initial(ret)
  c(ret, init_params)
}

##' Create initial conditions for the dualvax model
##' @name dualvax_initial
##' @title Initial conditions for the dualvax model
##' @param pars A parameter list created by [dualvax_params()]; from
##'   this list we will use the `N0`, `q`, `prev_Asl` and `prev_Ash`
##'   elements.
##' @return A list of initial conditions
##' @export
dualvax_initial <- function(pars) {
  # separate into 1:low and 2:high activity groups
  N0 <- round(pars$N0 * pars$q)
  U0 <- I0 <- A0 <- S0 <- T0 <- array(0, c(length(pars$beta), 2, pars$n_vax))
  # set initial asymptomatic prevalence in each group
  A0[, 1, 1] <- round(N0[1] * pars$prev_Asl)
  A0[, 2, 1] <- round(N0[2] * pars$prev_Ash)

  # set initial uninfecteds
  U0[, 1, 1] <- N0[1] - A0[, 1, 1]
  U0[, 2, 1] <- N0[2] - A0[, 2, 1]

  list(U0 = U0, I0 = I0, A0 = A0, S0 = S0, T0 = T0)
}

##' Create parameters to start the dualvax model from equilibrium
##' @name dualvax_initial_t
##' @title Create parameters to start the dualvax model from equilibrium
##' @param n_vax an integer indicating the number of vaccine compartments
##' @param comps a named list of arrays, as output by the model
##' (including U, I, A, S, T)
##' @return A list of initial conditions
##' @export
dualvax_initial_t <- function(n_vax, comps) {

  dim_comps <- dim(comps[["U"]])

  i_t <- dim_comps[1]
  n_par <- dim_comps[2]

  n_vax_input <- dim_comps[4]
  i_vax <- seq_len(min(n_vax,  n_vax_input))

  U0 <- I0 <- A0 <- S0 <- T0 <- array(0, c(n_par, 2, n_vax))
  # set compartments in each group
  U0[, , i_vax] <- comps$U[i_t, , , i_vax]
  I0[, , i_vax] <- comps$I[i_t, , , i_vax]
  A0[, , i_vax] <- comps$A[i_t, , , i_vax]
  S0[, , i_vax] <- comps$S[i_t, , , i_vax]
  T0[, , i_vax] <- comps$T[i_t, , , i_vax]

  list(U0 = U0, I0 = I0, A0 = A0, S0 = S0, T0 = T0)
}

vax_params0 <- function() {
  list(n_vax = 1,
       ve = 1,
       vs = 0,
       vd = 0,
       eff = 0,
       v = as.matrix(0),
       w = as.matrix(0))
}

vax_params1 <- function(z = 1 / 4, ve = 1, vs = 0, vd = 0, eff = 0.31) {
  list(n_vax = 2,
       ve    = c(1 - ve, ve),
       vs    = vs,
       vd    = vd,
       eff   = c(0, eff),
       v     = rbind(c(1, 0),
                    c(-1, 0)),
       w     = rbind(c(0,  z),
                     c(0, -z))
  )
}
