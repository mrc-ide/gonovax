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
