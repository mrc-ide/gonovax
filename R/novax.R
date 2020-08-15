
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


vax_params0 <- function() {
  list(n_vax = 1,
       ve = 1,
       vs = 0,
       vd = 0,
       eff = 0,
       v = as.matrix(0),
       w = as.matrix(0))
}
