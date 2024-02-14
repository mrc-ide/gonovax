##' @name vax_params_xvwr
##' @title create vaccination parameters for use in onevax_xvwr model
##' @inheritParams vax_params_xvwv
##' @param vea_revax scalar indicating efficacy of revaccination against
##'  acquisition (between 0-1)
##' @param vei_revax scalar indicating efficacy of revaccination against
##'  infectiousness (between 0-1)
##' @param ved_revax scalar indicating efficacy of revaccination against
##'  duration of infection (between 0-1)
##' @param ves_revax scalar indicating efficacy of revaccination against
##'  symptoms (between 0-1)
##' @param dur_revax duration of protection for revaccination,
##'  default to same as primary
##' @param primary_uptake scalar or numeric vector with same length as
##'  'gono_params' giving proportion of population undertaking primary
##'  vaccination as part of strategy
##' @param booster_uptake scalar or numeric vector with same length as
##'  'gono_params' giving proportion of population undertaking booster
##'  vaccination after primary vaccination protection has waned
##' @return A list parameters in the model input format
vax_params_xvwr <- function(vea = 0, vei = 0, ved = 0, ves = 0,
                            vea_revax = vea, vei_revax = vei,
                            ved_revax = ved, ves_revax = ves,
                            dur = 1e3, dur_revax = dur, primary_uptake = 0,
                            booster_uptake = primary_uptake,
                            strategy = NULL,
                            vbe = 0, t_stop = 99, n_diag_rec = 1) {

  assert_scalar_unit_interval(vea)
  assert_scalar_unit_interval(vei)
  assert_scalar_unit_interval(ved)
  assert_scalar_unit_interval(ves)
  assert_scalar_unit_interval(vea_revax)
  assert_scalar_unit_interval(vei_revax)
  assert_scalar_unit_interval(ved_revax)
  assert_scalar_unit_interval(ves_revax)
  assert_scalar_positive(dur)
  assert_scalar_positive(dur_revax)
  assert_scalar_unit_interval(primary_uptake)
  assert_scalar_unit_interval(booster_uptake)
  assert_scalar_unit_interval(vbe)
  assert_scalar_positive(t_stop)
  # waned vaccinees move to own stratum, but are eligible for re-vaccination
  # re-vaccination is into a fourth stratum (r)
  # 1:x -> 2:v -> 3:w <-> 4:r
  
  
  
  
  # generate indices for all strata and
  idx <- stratum_index_xvwr(n_diag_rec)
  
  n_vax <- idx$n_vax
  
  i_v <- c(idx$V, idx$R)
  
  i_w <- idx$V + n_diag_rec
  
  i <- seq_len(idx$n_vax)
  
  # diagnosed from
  i_diagnosedfrom <- i[i %% n_diag_rec != 0]
  #i_eligible <-  c(1, 3)
  
  # diagnosed to
  
  i_diagnosedto <- i[i %% n_diag_rec != 1]
  
  
  
  #i_eligible <- c(1, 3)
  #i_v <- c(2, 4)
  
  
  i_eligible_temp <- c(1:n_diag_rec, ( (2*n_diag_rec + 1) :3*n_diag_rec))
  i_v_temp <- c((1*(n_diag_rec) + 1): (2*(n_diag_rec)), (3*(n_diag_rec) + 1): (4*(n_diag_rec)))
  
  
  ## Could be implemented better
  if (length(strategy) > 0){
    if ( strategy == "VaH" | strategy == "VaHonly"){
      
      #Remove values corresponding to no diagnosis history
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
  
  
  n_group <- 2
  
  # create diagnosis history mapping
  diag_rec <- create_vax_map_branching(idx$n_vax, c(1, 1), i_diagnosedfrom, i_diagnosedto,
                                       set_vbe = FALSE, idx)
  

  # ensure duration is not divided by 0
  ved <- min(ved, 1 - 1e-10)
  ved_revax <- min(ved_revax, 1 - 1e-10)

  # If uptake of VbE > 0 consider that all adolescents are offered vaccine
  p <- set_strategy(strategy, vbe > 0)

  
  print("hello")
  
  print(i_eligible_temp)
  print(i_v_temp)
  
  # set up uptake matrix rows = groups, columns = vaccine strata
  u <- create_uptake_map(n_group = n_group, n_vax = n_vax,
                         primary_uptake = rep(primary_uptake, n_diag_rec),
                         booster_uptake = rep(booster_uptake, n_diag_rec),
                         i_eligible = i_eligible_temp, i_v = i_v_temp)
  
  
  print("hello 2")
  
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
    vea     = c(0, vea, 0, vea_revax),
    vei     = c(0, vei, 0, vei_revax),
    ved     = c(0, ved, 0, ved_revax),
    ves     = c(0, ves, 0, ves_revax),
    w       = create_waning_map(n_vax, i_v, i_w, 1 / dur, n_diag_rec),
    wd      = create_Diagnosiswaning_map(n_vax, 1 , n_diag_rec),
    vax_t   = c(0, t_stop),
    vax_y   = c(1, 0),
    diag_rec = diag_rec
  )
}

##' @name run_onevax_xvwr
##' @title run model with single vaccine for input parameter sets, either from
##' initialisation or from equilibrium, those with waned vaccines are eligible
##' for revaccination (R), and return to the R stratum
##' @param vea_revax scalar or numeric vector with same length as `gono_params`
##'  giving efficacy of revaccination against acquisition, default to same as
##'  primary
##' @param vei_revax scalar or numeric vector with same length as `gono_params`
##'  giving efficacy of revaccination against infectiousness, default to same as
##'  primary
##' @param ved_revax scalar or numeric vector with same length as `gono_params`
##'  giving efficacy of revaccination against duration of infection, default to
##'  same as primary
##' @param ves_revax scalar or numeric vector with same length as `gono_params`
##'  giving efficacy of revaccination against symptoms, default to same as
##'  primary
##' @param dur_revax scalar or numeric vector with same length as `gono_params`
##'  giving duration of protection for revaccination, default to same as primary
##' @param primary_uptake scalar or numeric vector with same length as
##'  'gono_params' giving proportion of population undertaking primary
##'  vaccination as part of strategy
##' @param booster_uptake scalar or numeric vector with same length as
##'  'gono_params' giving proportion of population undertaking booster
##'  vaccination after primary vaccination protection has waned.
##'   Defaults to supplied value of `primary_uptake`.
##' @inheritParams run_onevax_xvwv
##' @return A list of transformed model outputs
##' @export
run_onevax_xvwr <- function(tt, gono_params, init_params = NULL,
                            dur = 1e3, vea = 0, vei = 0, ved = 0, ves = 0,
                            dur_revax = dur,
                            vea_revax = vea, vei_revax = vei,
                            ved_revax = ved, ves_revax = ves,
                            vbe = 0, primary_uptake = 0, n_diag_rec = 1,
                            booster_uptake = primary_uptake, strategy = NULL,
                            t_stop = 99) {

  stopifnot(all(lengths(list(primary_uptake, booster_uptake, vea, vei, ved,
                             ves, dur, vea_revax, vei_revax, ved_revax,
                             ves_revax, dur_revax)) %in%
                  c(1, length(gono_params))))

  vax_params <- Map(vax_params_xvwr, primary_uptake = primary_uptake,
                    booster_uptake = booster_uptake, dur = dur,
                    vea = vea, vei = vei, ved = ved, ves = ves,
                    dur_revax = dur_revax,
                    vea_revax = vea_revax, vei_revax = vei_revax,
                    ved_revax = ved_revax, ves_revax = ves_revax, n_diag_rec = n_diag_rec,
                    MoreArgs = list(strategy = strategy,
                                    t_stop = t_stop, vbe = vbe))

  if (is.null(init_params)) {
    ret <- Map(run, gono_params = gono_params, vax_params = vax_params,
               MoreArgs = list(tt = tt))
  } else {
    ret <- Map(run, gono_params = gono_params, init_params = init_params,
               vax_params = vax_params,
               MoreArgs = list(tt = tt))
  }

  # name outputs
  ret <- lapply(ret, name_outputs, c("X", "V", "W", "R"))
  ret
}
