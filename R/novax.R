##' @name vax_params0
##' @title create vaccination parameters for use in novax model (null)
##' @return A list parameters in the model input format
##' @param dh integer giving the number of each X, V(erlang), and W stratum,
##' allowing tracking of diagnosis history. e.g for a dh = 2 and erlang = 1,
##' there will be X.I, X.II, V1.I, V1.II, W.I, W.II strata. Where '.I'
##' corresponds to never-diagnosed individuals and '.II' is for individuals
##' diagnosed at least once.
vax_params0 <- function(dh = 1) {
  n_group <- 2
  n_vax <- dh
  v <- array(0, dim = c(n_group, n_vax, n_vax))
  list(n_vax = dh,
       willing = 1,
       u_vbe = 0,
       u = v,
       vbe = v,
       vos = v,
       vod = v,
       vea = 0,
       ved = 0,
       ves = 0,
       vei = 0,
       w = as.matrix(0),
       D = 0,
       vax_t = c(0, 99),
       vax_y = c(0, 0))
}
