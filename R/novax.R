##' @name vax_params0
##' @title create vaccination parameters for use in novax model (null)
##' @return A list parameters in the model input format
##' @param n_diag_rec integer giving the number of each X, V(erlang), and W
##' @param years_history number of years that diagnosis history is recorded for
##' stratum, allowing tracking of diagnosis history. e.g for a n_diag_rec = 2
##' and erlang = 1, there will be X.I, X.II, V1.I, V1.II, W.I, W.II strata.
##' Where '.I' corresponds to never-diagnosed individuals and '.II' is for
##' individuals diagnosed at least once.
vax_params0 <- function(n_diag_rec = 1, years_history = 1) {
  n_group <- 2
  n_vax <- n_diag_rec

  willing <- rep(0, n_diag_rec)
  willing[1] <- 1

  v <- array(0, dim = c(n_group, n_vax, n_vax))
  list(n_vax = n_diag_rec,
       willing = willing,
       u_vbe = 0,
       u_d = v,
       u_s = v,
       u_pn = v,
       vbe = v,
       vos = v,
       vod = v,
       vopn = v,
       vea = rep(0, n_diag_rec),
       ved = rep(0, n_diag_rec),
       ves = rep(0, n_diag_rec),
       vei = rep(0, n_diag_rec),
       w = array(0, dim = c(n_diag_rec, n_diag_rec)),
       wd = create_diagnosis_waning_map(n_vax, 1 / years_history, n_diag_rec),
       D = 0,
       vax_t = c(0, 99),
       vax_y = c(0, 0),
       hesgroupmatrix = matrix(1 , n_vax, n_vax) ,
       stratum_doses = rep(0, n_vax)
  )
}


##' @name vax_params0_repeat
##' @title create vaccination parameters for use in novax model (null)
##' @return A list parameters in the model input format
##' @param n_diag_rec integer giving the number of each X, V(erlang), and W
##' @param years_history number of years that diagnosis history is recorded for
##' stratum, allowing tracking of diagnosis history. e.g for a n_diag_rec = 2
##' and erlang = 1, there will be X.I, X.II, V1.I, V1.II, W.I, W.II strata.
##' Where '.I' corresponds to never-diagnosed individuals and '.II' is for
##' individuals diagnosed at least once.
vax_params0_repeat <- function(n_diag_rec = 1, years_history = 1, hesgroups = 1) {
  n_group <- 2
  n_vax <- n_diag_rec * hesgroups
  
  willing <- rep(0, n_diag_rec*hesgroups)
  willing[2*(1:n_diag_rec)-1] <- 1
  
  hesgroupmatrix = matrix(0 , n_vax, n_vax)
  
  for (i in 1:nrow(hesgroupmatrix)){
    for (j in 1:nrow(hesgroupmatrix)){
      if ( i %% hesgroups == j %% hesgroups){
        hesgroupmatrix[i,j] = 1}
    }
  }
  
  v <- array(0, dim = c(n_group, n_vax, n_vax))
  list(n_vax = n_diag_rec,
       willing = willing,
       u_vbe = 0,
       u_d = v,
       u_s = v,
       u_pn = v,
       vbe = v,
       vos = v,
       vod = v,
       vopn = v,
       vea = rep(0, n_diag_rec * hesgroups),
       ved = rep(0, n_diag_rec * hesgroups),
       ves = rep(0, n_diag_rec * hesgroups),
       vei = rep(0, n_diag_rec * hesgroups),
       w = array(0, dim = c(n_diag_rec * hesgroups, n_diag_rec * hesgroups)),
       wd = create_diagnosis_waning_map(n_vax, 1 / years_history, n_diag_rec * hesgroups),
       D = 0,
       vax_t = c(0, 99),
       vax_y = c(0, 0),
       hesgroupmatrix = hesgroupmatrix)
}