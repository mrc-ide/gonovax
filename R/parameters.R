demographic_params <- function() {
  list(N0 = 6e5,
       q = c(0.85, 0.15), # proportion in group L
       p = c(0.6, 15.6), # partner change rate in group L/H
       enr = 12000,
       exr = 1 / 50
  )
}
##' @name gono_params
##' @title Posterior parameters of gonorrhoea natural history
##' @param n an integer vector (or value) containing the indices of the required
##' parameter sets (1:982). If `n = NULL` the full parameter set is returned
##' @return A data frame of parameters
##' @export

gono_params <- function(n = NULL) {
  if (is.null(cache$gono_params)) {
    cache$gono_params <-
      read_csv(gonovax_file("extdata/gono_params.csv"))
  }

  pars <- cache$gono_params
  n_pars <- nrow(pars)
  # if n not supplied, return all parameters
  i  <- n %||% seq_len(n_pars)
  # limit to parameter sets available
  i <- i[(i > 0) & (i <= n_pars)]

  pars[i, ]
}

##' Create initial conditions for the model
##' @name initial_params
##' @title Initial conditions for the model
##' @param pars A parameter list containing `N0`, `q`, `prev_Asl` and `prev_Ash`
##'   elements.
##' @param n_vax an integer indicating the number of vaccine compartments
##' @return A list of initial conditions
##' @export
initial_params <- function(pars, n_vax = 1) {
  # separate into 1:low and 2:high activity groups
  N0 <- round(pars$N0 * pars$q)
  n_par <- length(pars$prev_Asl)
  U0 <- I0 <- A0 <- S0 <- T0 <- array(0, c(n_par, 2, n_vax))
  # set initial asymptomatic prevalence in each group
  A0[, 1, 1] <- round(N0[1] * pars$prev_Asl)
  A0[, 2, 1] <- round(N0[2] * pars$prev_Ash)
  
  # set initial uninfecteds
  U0[, 1, 1] <- N0[1] - A0[, 1, 1]
  U0[, 2, 1] <- N0[2] - A0[, 2, 1]
  
  list(U0 = U0, I0 = I0, A0 = A0, S0 = S0, T0 = T0)
}

##' Create initial conditions based on previous model run
##' @name restart_params
##' @title Create initial conditions to start the model from the end of a run
##' @param y a transformed model run output
##' @param n_vax an integer indicating the number of vaccine compartments,
##' consistent with the input
##' @return A list of initial conditions to restart a model with n_vax
##' vaccination levels
##' @export
restart_params <- function(y, n_vax = NULL) {
  
  dim_y <- dim(y[["U"]])
  
  i_t <- dim_y[1]
  n_par <- dim_y[2]
  n_vax <- n_vax %||% dim_y[4]
  
  n_vax_input <- dim_y[4]
  i_vax <- seq_len(min(n_vax,  n_vax_input))
  
  U0 <- I0 <- A0 <- S0 <- T0 <- array(0, c(n_par, 2, n_vax))
  # set compartments in each group
  U0[, , i_vax] <- y$U[i_t, , , i_vax]
  I0[, , i_vax] <- y$I[i_t, , , i_vax]
  A0[, , i_vax] <- y$A[i_t, , , i_vax]
  S0[, , i_vax] <- y$S[i_t, , , i_vax]
  T0[, , i_vax] <- y$T[i_t, , , i_vax]
  
  list(U0 = U0, I0 = I0, A0 = A0, S0 = S0, T0 = T0)
}

##' @name model_params
##' @title Parameters for the dualvax model
##' @param gono_params A dataframe of natural history parameters
##' @param demographic_params A dataframe of demographic parameters
##' @param vax_params A vector of vaccination params
##' @param init_params A list of starting conditions
##' @return A list of inputs to the model many of which are fixed and
##'   represent data. These correspond largely to `user()` calls
##'   within the odin code, though some are also used in processing
##'   just before the model is run.
##' @export
model_params <- function(gono_params = NULL,
                         demographic_params = NULL,
                         init_params = NULL,
                         vax_params = NULL) {
  gono_params <- gono_params %||% gono_params(1)
  demographic_params <- demographic_params %||% demographic_params()
  ret <- c(demographic_params, gono_params)
  vax_params <- vax_params %||% vax_params0()
  init_params <- init_params %||% initial_params(ret, n_vax = vax_params$n_vax)
  c(ret, init_params, vax_params)
}