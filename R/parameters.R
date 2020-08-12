demographic_params <- function() {
  list(N0 = 6e5,
       q = c(0.85, 0.15), # proportion in group L
       p = c(0.6, 15.6), # partner change rate in group L/H
       enr = 12000,
       exr = 1 / 50
  )
}

gono_params <- function(n = NULL) {
  if (is.null(cache$gono_params)) {
    cache$gono_params <-
      read_csv(gonovax_file("extdata/gono_params.csv"))
  }
  pars <- cache$gono_params
  # select first n parameter sets, if n supplied
  i  <- n %||% nrow(pars)
  pars[seq_len(i), ]
}
