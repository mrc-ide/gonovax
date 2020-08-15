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
