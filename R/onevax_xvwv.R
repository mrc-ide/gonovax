##' @name vax_params_xvwv
##' @title create vaccination parameters for use in onevax_xvwv model
##' @inheritParams vax_params_xvw
##' @return A list parameters in the model input format
vax_params_xvwv <- function(vea = 0, vei = 0, ved = 0, ves = 0,
                                      dur = 1e3, uptake = 0, strategy = NULL,
                                      vbe = 0, t_stop = 99, n_diag_rec = 1) {
  
  assert_scalar_unit_interval(vea)
  assert_scalar_unit_interval(vei)
  assert_scalar_unit_interval(ved)
  assert_scalar_unit_interval(ves)
  assert_scalar_positive(dur)
  assert_scalar_unit_interval(uptake)
  assert_scalar_unit_interval(vbe)
  assert_scalar_positive(t_stop)
  # waned vaccinees move to own stratum, but are eligible for re-vaccination
  # 1:x -> 2:v <-> 3:w
  
  
  # generate indices for all strata and
  idx <- stratum_index_xvw_trial(1, n_diag_rec)
  
  n_vax <- idx$n_vax
  
  i_v <- idx$V
  
  i_w <- idx$V + n_diag_rec
  
  i <- seq_len(idx$n_vax)
  
  
  # diagnosed from
  i_diagnosedfrom <- i[i %% n_diag_rec != 0]
  #i_eligible <-  c(1, 3)
  
  # diagnosed to
  
  i_diagnosedto <- i[i %% n_diag_rec != 1]
  
  n_group <- 2

  # create diagnosis history mapping
  diag_rec <- create_vax_map_branching(idx$n_vax, c(1, 1), i_diagnosedfrom, i_diagnosedto,
                                       set_vbe = FALSE, idx)
  

  
  ### no longer set up properly - needs fixed in line with xpvwrh changes! 17 July
  
  i_eligible_temp <- c(1:n_diag_rec, (2*(n_diag_rec) + 1): (3*(n_diag_rec)))
  i_v_temp <- c((1*(n_diag_rec) + 1): (2*(n_diag_rec)), ((1*(n_diag_rec) + 1): (2*(n_diag_rec))))
  
  
  
  
  ## Could be implemented better
  if (length(strategy) > 0){
    if ( strategy == "VaH"){
      i_eligible_temp2 <- i_eligible_temp[-c(1, (n_diag_rec+1))]
      i_v_temp2 <- i_v_temp[-c(1,(n_diag_rec+1))]
    }
    else{
      i_eligible_temp2 = i_eligible_temp
      i_v_temp2 = i_v_temp
    }
  }
  else{
    i_eligible_temp2 = i_eligible_temp
    i_v_temp2 = i_v_temp
  }
  
  
  
  # compartments to which vaccine efficacy applies
  ve <- c(rep(0, n_diag_rec), rep(1, 1* n_diag_rec), rep(0, n_diag_rec))
  ved <- min(ved, 1 - 1e-10) # ensure duration is not divided by 0
  
  # If uptake of VbE > 0 consider that all adolescents are offered vaccine
  p <- set_strategy(strategy, vbe > 0)
  
  
  
  # set up uptake matrix rows = groups, columns = vaccine strata
  u <- create_uptake_map(n_group = n_group, n_vax = n_vax,
                         primary_uptake = rep(uptake, n_diag_rec),
                         booster_uptake = rep(uptake, n_diag_rec),
                         i_eligible = i_eligible_temp, i_v = i_v_temp)
  
  
  willing = rep(0, n_vax)
  willing[1] = 1
 
  
  list(n_vax   = n_vax,
       willing = willing,
       u_s     = u,
       u_d     = u,
       u_vbe   = vbe,
       vbe     = create_vax_map(n_vax, p$vbe, i_eligible_temp, i_v_temp),
       vod     = create_vax_map(n_vax, p$vod, i_eligible_temp, i_v_temp),
       vos     = create_vax_map(n_vax, p$vos, i_eligible_temp2, i_v_temp2),
       vea     = vea * ve,
       vei     = vei * ve,
       ved     = ved * ve,
       ves     = ves * ve,
       w       = create_waning_map(n_vax, i_v, i_w, 1 / dur, n_diag_rec),
       wd      = create_Diagnosiswaning_map(n_vax, 1 , n_diag_rec),
       vax_t   = c(0, t_stop),
       vax_y   = c(1, 0),
       diag_rec = diag_rec,
       notification_param = 0
  )
}

##' @name run_onevax_xvwv
##' @title run model with single vaccine for input parameter sets, either from
##' initialisation or from equilibrium, those with waned vaccines are eligible
##' for revaccination, and return to the V compartment
##' @param gono_params list of gono params
##' @param vea scalar or numeric vector with same length as `gono_params` giving
##'  efficacy of the vaccine against acquisition (between 0-1)
##' @param vei scalar or numeric vector with same length as `gono_params` giving
##'  efficacy of the vaccine against infectiousness (between 0-1)
##' @param ved scalar or numeric vector with same length as `gono_params` giving
##'  efficacy of the vaccine against duration (between 0-1)
##' @param ves scalar or numeric vector with same length as `gono_params` giving
##'  efficacy of the vaccine against symptoms (between 0-1)
##' @param dur  scalar or numeric vector with same length as `gono_params`
##'  giving duration of the vaccine (in years)
##' @param uptake  scalar or numeric vector with same length as `gono_params`
##'  giving pc of population vaccinated as part of strategy
##'  @param n_diag_rec integer giving the number of each X, V(erlang), and W
##' stratum, allowing tracking of diagnosis history. e.g for a n_diag_rec = 2
##' and erlang = 1, there will be X.I, X.II, V1.I, V1.II, W.I, W.II strata.
##' Where '.I' corresponds to never-diagnosed individuals and '.II' is for
##' individuals diagnosed at least once.
##' @inheritParams run
##' @inheritParams vax_params_xvwv
##' @export
run_onevax_xvwv <- function(tt, gono_params, init_params = NULL, dur = 1e3,
                                      vea = 0, vei = 0, ved = 0, ves = 0, vbe = 0, n_diag_rec = 1,
                                      uptake = 0, strategy = NULL,
                                      t_stop = 99) {
  
  
  
  stopifnot(all(lengths(list(uptake, vea, vei, ved, ves, dur)) %in%
                  c(1, length(gono_params))))
  
  
  
  vax_params <- Map(vax_params_xvwv, uptake = uptake, dur = dur,
                    vea = vea, vei = vei, ved = ved, ves = ves, n_diag_rec = n_diag_rec,
                    MoreArgs = list(strategy = strategy, t_stop = t_stop,
                                    vbe = vbe))

  
  if (is.null(init_params)) {
    
    ret <- Map(run, gono_params = gono_params, vax_params = vax_params,
               MoreArgs = list(tt = tt))
    
    
  } else {
    ret <- Map(run, gono_params = gono_params, init_params = init_params,
               vax_params = vax_params,
               MoreArgs = list(tt = tt))
    
    
  }
  
  # name outputs
  ret <- lapply(ret, name_outputs, gen_labels(1, n_diag_rec))
  
  ret
}
