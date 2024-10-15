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
  # unvaccinated 1:n_diag_rec (x)
  # vaccinated from n_diag_rec+1 to 2*n_diag_rec (v)
  # waned from  2*n_diag_rec+1 to 3*n_diag_rec (w)

  # generate indices for all strata and
  idx <- stratum_index_xvwv(1, n_diag_rec = n_diag_rec, strategy = strategy)

  n_vax <- idx$n_vax

  i_v <- idx$V
  i_w <- idx$W

  n_group <- 2

  # create diagnosis history mapping
  diag_rec <- create_vax_map_branching(idx$n_vax, c(1, 1), idx$diagnosedfrom,
                                       idx$diagnosedto, set_vbe = FALSE, idx)

  # compartments to which vaccine efficacy applies
  ve <- ifelse(seq_len(n_vax) %in% idx$V, 1, 0)
  ved <- min(ved, 1 - 1e-10) # ensure duration is not divided by 0

  # If uptake of VbE > 0 consider that all adolescents are offered vaccine
  p <- set_strategy(strategy, vbe > 0)

  # set up uptake matrix rows = groups, columns = vaccine strata
  u_s <- create_uptake_map_xvwv(n_group = n_group, n_vax = n_vax,
                                primary_uptake = uptake,
                                booster_uptake = uptake,
                                idx, n_diag_rec = n_diag_rec,
                                screening_or_diagnosis = "screening")

  u_d <- create_uptake_map_xvwv(n_group = n_group, n_vax = n_vax,
                                primary_uptake = uptake,
                                booster_uptake = uptake,
                                n_diag_rec = n_diag_rec,
                                idx, screening_or_diagnosis = "diagnosis")

  u_pn <- u_s

  if (sum(p$vod) > 0) {
    #vaccination on diagnosis occuring, so need to scale down diag_rec
    diag_rec[, idx$X, ] <- (1 - uptake) * diag_rec[, idx$X, ]
    diag_rec[, idx$W, ] <- (1 - uptake) * diag_rec[, idx$W, ]
  }

  willing <- rep(0, n_vax)
  willing[1] <- 1

  list(n_vax   = n_vax,
    willing = willing,
    u_s = u_s,
    u_d = u_d,
    u_pn = u_pn,
    u_vbe = vbe,
    vbe     = create_vax_map(n_vax, p$vbe, idx$vaccinatedfrom_vbe,
                             idx$vaccinatedto_vbe),
    vod     = create_vax_map(n_vax, p$vod, idx$vaccinatedfrom_vod,
                             idx$vaccinatedto_vod),
    vos     = create_vax_map(n_vax, p$vos, idx$vaccinatedfrom_vos,
                             idx$vaccinatedto_vos),
    vopn    = create_vax_map(n_vax, p$vopn, idx$vaccinatedfrom_vopn,
                             idx$vaccinatedto_vopn),
    vea = vea * ve,
    vei = vei * ve,
    ved = ved * ve,
    ves = ves * ve,
    w = create_waning_map(n_vax, i_v, i_w, 1 / dur, n_diag_rec),
    wd = create_diagnosis_waning_map(n_vax, 1, n_diag_rec),
    vax_t = c(0, t_stop),
    vax_y = c(1, 0),
    diag_rec = diag_rec,
    hesgroupmatrix = idx$hesgroupmatrix
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
##'  @param n_diag_rec integer for the number of diagnosis history substrata
##' @inheritParams run
##' @inheritParams vax_params_xvwv
##' @export
run_onevax_xvwv <- function(tt, gono_params, init_params = NULL, dur = 1e3,
                            vea = 0, vei = 0, ved = 0, ves = 0, vbe = 0,
                            n_diag_rec = 1, uptake = 0, strategy = NULL,
                            t_stop = 99) {

  stopifnot(all(lengths(list(uptake, vea, vei, ved, ves, dur)) %in%
                  c(1, length(gono_params))))

  vax_params <- Map(vax_params_xvwv, uptake = uptake, dur = dur,
                    vea = vea, vei = vei, ved = ved, ves = ves,
                    n_diag_rec = n_diag_rec,
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
