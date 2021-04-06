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
##' @return A list of parameters
##' @export

gono_params <- function(n = NULL) {
  if (is.null(cache$gono_params)) {
    gp <- read_csv(gonovax_file("extdata/gono_params.csv"))
    cache$gono_params <-
      lapply(seq_len(nrow(gp)), function(i) transform0(gp[i, ]))
  }

  pars <- cache$gono_params
  n_pars <- length(pars)
  # if n not supplied, return all parameters
  i  <- n %||% seq_len(n_pars)
  # limit to parameter sets available
  i <- i[(i > 0) & (i <= n_pars)]

  pars[i]
}

transform0 <- function(pars) {
    # reformat time-varying pars
    pars <- as.list(pars)
    pars$tt <- c(0, 1e3)
    pars$beta_t <-  rep(pars$beta, 2)
    pars$eta_l_t <- rep(pars$eta, 2)
    pars$eta_h_t <- pars$eta_l_t
    pars$beta <- pars$eta <- NULL

    pars
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
  U0 <- I0 <- A0 <- S0 <- T0 <- array(0, c(2, n_vax))

  # set initial asymptomatic prevalence in each group
  A0[, 1] <- round(N0 * c(pars$prev_Asl, pars$prev_Ash))

  # set initial uninfecteds
  U0[, 1] <- N0 - A0[, 1]

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
  n_vax <- n_vax %||% dim_y[3]

  n_vax_input <- dim_y[3]
  i_vax <- seq_len(min(n_vax,  n_vax_input))

  U0 <- I0 <- A0 <- S0 <- T0 <- array(0, c(2, n_vax))
  # set compartments in each group
  U0[, i_vax] <- y$U[i_t, , i_vax]
  I0[, i_vax] <- y$I[i_t, , i_vax]
  A0[, i_vax] <- y$A[i_t, , i_vax]
  S0[, i_vax] <- y$S[i_t, , i_vax]
  T0[, i_vax] <- y$T[i_t, , i_vax]

  list(U0 = U0, I0 = I0, A0 = A0, S0 = S0, T0 = T0, t = y$t[i_t])
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

##' @name create_vax_map
##' @title Create mapping for movement between strata due to vaccination
##' @param n_vax Integer in (0, 5) denoting total number of strata
##' @param v numeric indicating % of population vaccinated on screening a scalar
##' will apply to both activity groups, a vector of length 2 will apply to each
##' activity group separately
##' @param i_u indices of strata eligible for vaccination
##' @param i_v indices of strata being vaccinated
##' @return an array of the mapping

create_vax_map <- function(n_vax, v, i_u, i_v) {

  # ensure vaccine input is of correct length
  n_group <- 2
  stopifnot(length(v) %in% c(1, n_group))
  stopifnot(all((v >= 0) & (v <= 1)))
  stopifnot(length(i_v) == length(i_u))
  stopifnot(max(i_u, i_v) <= n_vax)

  # set up vaccination matrix
  vax_map <- array(0, dim = c(n_group, n_vax, n_vax))

  for (i in seq_along(i_u)) {
    vax_map[, i_u[i], i_u[i]] <-  v
    vax_map[, i_v[i], i_u[i]] <- -v
  }

  vax_map
  }

##' @name create_waning_map
##' @title Create mapping for movement between strata due to vaccine waning
##' @param n_vax Integer in (0, 5) denoting total number of strata
##' @param i_v indices of strata being vaccinated
##' @param i_w Integer in (0, 5) denoting which stratum receives waned vaccinees
##' @param z Scalar denoting rate of waning
##' @return an array of the mapping

create_waning_map <- function(n_vax, i_v, i_w, z) {

  stopifnot(z > 0)
  stopifnot(length(z) %in% c(1, length(i_v)))
  stopifnot(length(i_w) == 1)
  # set up waning map
  w <- array(0, dim = c(n_vax, n_vax))

  w[i_w, i_v] <-  z

  for (i in i_v) {
    w[i, i] <- -w[i_w, i]
  }
  w
}

##' @name set_strategy
##' @title Create mapping from vaccination strategy / uptake to model params
##' @param strategy one of "ve", "va", "vd" or "vt"
##' @param uptake single numeric between 0-1
##' @return a list with elements vd and vs
##' @export
set_strategy <- function(strategy, uptake) {
  if (strategy == "ve") {
    vs <- vd <- 0
  } else if (strategy == "vd") {
    vd <- uptake
    vs <- 0
  } else if (strategy == "va") {
    vd <- uptake
    vs <- uptake
  } else if (strategy == "vt") {
    vd <- uptake
    vs <- c(0, uptake)
  }
  list(vd = vd, vs = vs)
}
