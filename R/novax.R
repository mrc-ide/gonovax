
##' @name novax_equilib
##' @title Equilibrium state of compartmental model with no vaccination
##' @param n an integer vector (or value) containing the indices of
##' corresponding parameter set (1:982). If `n = NULL` the equilibrium positions
##' for the full parameter set are returned
##' @return A list of transformed model outputs
##' @export
novax_equilib <- function(n = NULL) {

  if (is.null(cache$novax_equilib)) {
    cache$novax_equilib <-
      readRDS(gonovax_file("extdata/novax_equilib.rds"))
  }

  y <- cache$novax_equilib
  n_y <- length(y)
  # if n not supplied, return all parameters
  i  <- n %||% seq_len(n_y)
  # limit to parameter sets available
  i <- i[(i > 0) & (i <= n_y)]
  if (length(i) == 1) return(y[[i]])
  y[i]
}

##' @name novax_baseline
##' @title Equilibrium state of compartmental model with no vaccination
##' @param n an integer vector (or value) containing the indices of
##' corresponding parameter set (1:982). If `n = NULL` the equilibrium positions
##' for the full parameter set are returned
##' @param t an integer between 1 and 10 of the year to return
##' @return A list containing incid and cum_incid at time t for parameter sets n
##' @export
novax_baseline <- function(n = NULL, t) {

  if (!all(t %in% seq_len(10))) stop("t must be an integer between 1 and 10")

  if (is.null(cache$novax_baseline)) {
    cache$novax_baseline <-
      readRDS(gonovax_file("extdata/novax_baseline.rds"))
  }

  y <- cache$novax_baseline
  n_y <- ncol(y[[1]])
  # if n not supplied, return all parameters
  i  <- n %||% seq_len(n_y)
  # limit to parameter sets available
  i <- i[(i > 0) & (i <= n_y)]

  y$incid <- y$incid[t, i]
  y$cum_incid <- y$cum_incid[t, i]
  y
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


run_novax_int <- function(n = NULL, tt, equilib = FALSE) {
  init_params <- NULL
  if (equilib) {
    if (length(n) > 1) stop("if length(n) > 1, equilib must be FALSE")
    init_params <- restart_params(novax_equilib(n))
  }

  pars <- model_params(gono_params = gono_params(n),
                       init_params = init_params)
  mod <- model(user = pars, unused_user_action = FALSE)
  y <- mod$run(tt)
  mod$transform_variables(y)
}

##' @name run_novax
##' @title run model without vaccination for input parameter sets, either from
##' initialisation or from equilibrium
##' @param n an integer vector (or value) containing the indices of
##' corresponding parameter set (1:982). If `n = NULL` the equilibrium positions
##' for the full parameter set are returned
##' @param tt a numeric vector of times at which the model state is output
##' @param equilib a logical indicating whether to run from equilibrium, default
##' is `FALSE`, i.e. run from initial conditions
##' @return A list of transformed model outputs
##' @export
run_novax <- function(n = NULL, tt, equilib = FALSE) {
  n <- n %||% seq_len(nrow(gono_params()))
  lapply(n, run_novax_int, tt = tt, equilib = equilib)
}
