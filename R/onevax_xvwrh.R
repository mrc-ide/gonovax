##' Create initial conditions for the model
##' @name initial_params_xvwrh
##' @title Initial conditions for the model
##' @inheritParams initial_params
##' @param hes proportion of population vaccine hesitant
##' @param n_diag_rec integer for the number of diagnosis history substrata
##' @return A list of initial conditions
##' @export
initial_params_xvwrh <- function(pars, coverage = 0, hes = 0, n_diag_rec = 1) {
  assert_scalar_unit_interval(coverage)

  n_vax <- 5 * n_diag_rec
  willing <- 1 - hes
  x_init <- willing * (1 - coverage)
  v_init <- willing * coverage

  cov <- c(x_init, rep(0, n_diag_rec - 1), v_init, rep(0, n_diag_rec - 1),
           rep(0, n_diag_rec), rep(0, n_diag_rec), hes, rep(0, n_diag_rec - 1))


  stopifnot(length(cov) == n_vax)
  stopifnot(sum(cov) == 1)

  U0 <- I0 <- A0 <- S0 <- T0 <- array(0, c(2, n_vax))

  # separate into 1:low and 2:high activity groups and by coverage
  N0 <- pars$N0 * outer(pars$q, cov)
  # set initial asymptomatic prevalence in each group (X AND H)
  A0[, 1] <- round(N0[, 1] * c(pars$prev_Asl, pars$prev_Ash))

  indextemp <- 4 * n_diag_rec + 1
  A0[, indextemp] <- round(N0[, n_vax] * c(pars$prev_Asl, pars$prev_Ash))

  # set initial uninfecteds
  U0 <- round(N0) - A0

  list(U0 = U0, I0 = I0, A0 = A0, S0 = S0, T0 = T0)

}

##' @name vax_params_xvwrh
##' @title create vaccination parameters for use in onevax_xvwrh model
##' @inheritParams vax_params_xvwr
##' @param hes proportion of population vaccine hesitant
##' @return A list parameters in the model input format
vax_params_xvwrh <- function(vea = 0, vei = 0, ved = 0, ves = 0,
                             vea_revax = vea, vei_revax = vei,
                             ved_revax = ved, ves_revax = ves,
                             dur = 1e3, dur_revax = dur, primary_uptake = 0,
                             booster_uptake = primary_uptake, strategy = NULL,
                             vbe = 0, t_stop = 99, hes = 0,
                             n_diag_rec = 1) {

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
  # unvaccinated 1:n_diag_rec (x)
  # vaccinated from n_diag_rec+1 to 2*n_diag_rec (v)
  # waned from  2*n_diag_rec+1 to 3*n_diag_rec (w)
  # revaccinated from 3*n_diag_rec+1 to 4*n_diag_Rec (R)
  # hesitant from 5*n_diag_rec+1 to 6*n_diag_rec (H)
  # a proportion of all 'n' exist only in the hesitant compartment
  # there is no movement between the willing (x,v,w,r) and hesitant (h)

  # generate indices for all strata and
  idx <- stratum_index_xvwrh(1, n_diag_rec = n_diag_rec, strategy = strategy)

  n_vax <- idx$n_vax

  i_v <- c(idx$V, idx$R)
  i_w <- idx$W
  i <- seq_len(idx$n_vax)

  #number of compartments
  n_group <- 2

  # create diagnosis history mapping
  diag_rec <- create_vax_map_branching(idx$n_vax, c(1, 1), idx$diagnosedfrom,
                                       idx$diagnosedto, set_vbe = FALSE, idx)

  # ensure duration is not divided by 0
  ved <- min(ved, 1 - 1e-10)
  ved_revax <- min(ved_revax, 1 - 1e-10)

  # If uptake of VbE > 0 consider that all adolescents are offered vaccine
  p <- set_strategy(strategy, vbe > 0)

  # set up uptake matrix rows = groups, columns = vaccine strata
  u_s <- create_uptake_map_xvwr(n_group = n_group, n_vax = n_vax,
                                primary_uptake = primary_uptake,
                                booster_uptake = booster_uptake,
                                n_diag_rec = n_diag_rec,
                                idx, screening_or_diagnosis = "screening")

  u_d <- create_uptake_map_xvwr(n_group = n_group, n_vax = n_vax,
                                primary_uptake = primary_uptake,
                                booster_uptake = booster_uptake,
                                n_diag_rec = n_diag_rec,
                                idx, screening_or_diagnosis = "diagnosis")

  if (sum(p$vod) > 0) {
    #vaccination on diagnosis occuring, so need to scale down diag_rec
    diag_rec[, idx$X, ] <- (1 - primary_uptake) * diag_rec[, idx$X, ]
    diag_rec[, idx$W, ] <- (1 - booster_uptake) * diag_rec[, idx$W, ]
  }


  willing <- c((1 - hes), rep(0, n_vax - 1 - n_diag_rec), hes,
               rep(0, n_diag_rec - 1))

  list(n_vax   = n_vax,
    willing = willing,
    u_s     = u_s,
    u_d     = u_d,
    u_vbe   = vbe,
    vbe     = create_vax_map(n_vax, p$vbe, idx$vaccinatedfrom_vbe,
                             idx$vaccinatedto_vbe),
    vod     = create_vax_map(n_vax, p$vod, idx$vaccinatedfrom_vod,
                             idx$vaccinatedto_vod),
    vos     = create_vax_map(n_vax, p$vos, idx$vaccinatedfrom_vos,
                             idx$vaccinatedto_vos),
    vea     = c(0, vea, 0, vea_revax, 0),
    vei     = c(0, vei, 0, vei_revax, 0),
    ved     = c(0, ved, 0, ved_revax, 0),
    ves     = c(0, ves, 0, ves_revax, 0),
    w       = create_waning_map(n_vax, i_v, i_w, 1 / dur, n_diag_rec),
    wd      = create_diagnosis_waning_map(n_vax, 1, n_diag_rec),
    vax_t   = c(0, t_stop),
    vax_y   = c(1, 0),
    diag_rec = diag_rec
  )
}

##' @name run_onevax_xvwrh
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
##' @param hes Proportion of individuals in the population who are vaccine
##'  hesitant
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
run_onevax_xvwrh <- function(tt, gono_params, init_params = NULL,
                             dur = 1e3, vea = 0, vei = 0, ved = 0, ves = 0,
                             dur_revax = dur,
                             vea_revax = vea, vei_revax = vei,
                             ved_revax = ved, ves_revax = ves,
                             vbe = 0, primary_uptake = 0,
                             booster_uptake = primary_uptake, strategy = NULL,
                             t_stop = 99, hes = 0,
                             n_diag_rec = 1) {

  stopifnot(all(lengths(list(booster_uptake, primary_uptake, vea, vei,
                             ved, ves, dur, vea_revax, vei_revax, ved_revax,
                             ves_revax, dur_revax)) %in%
                  c(1, length(gono_params))))

  vax_params <- Map(vax_params_xvwrh, primary_uptake = primary_uptake,
                    booster_uptake = booster_uptake, dur = dur,
                    vea = vea, vei = vei, ved = ved, ves = ves,
                    dur_revax = dur_revax,
                    vea_revax = vea_revax, vei_revax = vei_revax,
                    ved_revax = ved_revax, ves_revax = ves_revax, hes = hes,
                    n_diag_rec = n_diag_rec,
                    MoreArgs = list(strategy = strategy,
                                    t_stop = t_stop, vbe = vbe))

  if (is.null(init_params)) {
    pars <- lapply(gono_params, model_params)
    init_params <- lapply(pars, initial_params_xvwrh, hes = hes)
  }

  ret <- Map(run, gono_params = gono_params, vax_params = vax_params,
             init_params = init_params,
             MoreArgs = list(tt = tt))

  # name outputs
  ret <- lapply(ret, name_outputs, c("X", "V", "W", "R", "H"))
  ret
}
