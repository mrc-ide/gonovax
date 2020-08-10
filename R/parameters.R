gonovax_parameters <- function() {
  ql <- 0.85
  ret <- list(N0 = 6e5,
       ql = ql, # proportion in group L
       p = c(0.6, 15.6), # partner change rate in group L/H
       enr = 12000 * c(ql, 1- ql),
       exr = 1/50
  )
  gono_params <- gono_params()[1, ]
  c(ret, gono_params)
}



