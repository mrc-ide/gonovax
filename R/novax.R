
##' @name vax_params0
##' @title create vaccination parameters for use in novax model (null)
##' @return A list parameters in the model input format
vax_params0 <- function() {
  v <- array(0, dim = c(2, 1, 1))
  list(n_vax = 1,
       vbe = v,
       vos = v,
       vod = v,
       vea = 0,
       ved = 0,
       ves = 0,
       vei = 0,
       w = as.matrix(0),
       vax_t = c(0, 99),
       vax_y = c(0, 0))
}
