### gonovax population model whereby individuals in X become vaccinated
# under the given vaccination strategy, with proportion r1 * (1 - r2) accepting
# 1 dose --> P or proportion r1 * r2 accepting both doses --> V
# Individuals with 1 dose who failed to register for dose 2 straight away but
# then go to the clinic for screening or treatment can then get their dose 2,
# at an uptake of r2_p
# Individuals in P who do not accept dose 2 a the clinic or who do not
# need to attend the clinic wane back to X, and can then go to P or V again
# after (assumed to be immunologically naive)
# Individuals in V wane to W where they can get a booster vaccination to R
# This is a modificaiton of the onevax_xvwrh model where indiviuals pass
# linearly through the strata, now there is a branching after X

##' Create initial conditions for the model
##' @name initial_params_xpvwrh
##' @title Initial conditions for the model
##' @inheritParams initial_params
##' @param hes proportion of population vaccine hesitant
##' @param t number of years, only use when using function outside of
##' run_onevax_xpvwrh() to generate initial conditions for tests
##' @param coverage_p partial (one-dose) vaccine coverage of the population
##' already present (as a proportion)
##' @param coverage_v two-dose vaccine coverage of the population already
##' present (as a proportion)
##' @param n_erlang integer giving the number of transitions that need to be
##'  made through vaccine-protected strata until that protection has waned
##' @param n_diag_rec  number of diagnosis history strata
##' @return A list of initial conditions
##' @export
initial_params_xpvwrh <- function(pars, coverage_p = 0, coverage_v = 0, hes = 0,
                                  t = FALSE, n_erlang = 1, n_diag_rec = 1) {

  idx <- stratum_index_xpvwrh(n_erlang, n_diag_rec)


  if (coverage_p + coverage_v + hes > 1) {
    stop("sum of coverages and/or hesitancy must not exceed 1")
  } else if (coverage_p + coverage_v == 1) {
    stop("You cannot have 100% coverage for initial conditions")
  } else {

    assert_scalar_unit_interval(coverage_p)
    assert_scalar_unit_interval(coverage_v)
    n_vax <- idx$n_vax
    willing <- 1 - hes
    x_init <- willing * (1 - coverage_p - coverage_v)
    p_init <- willing * coverage_p
    v_init <- willing * coverage_v

    ## X[1], P[n_erlang], V[n_erlang], W[1], R[n_erlang], H[1]
    cov <- c(x_init, rep(0, n_diag_rec - 1), p_init,
             rep(0, n_diag_rec * n_erlang - 1), v_init,
             rep(0, n_diag_rec * n_erlang - 1), rep(0, n_diag_rec),
             rep(0, n_diag_rec * n_erlang), hes, rep(0, n_diag_rec - 1))

    stopifnot(length(cov) == n_vax)
    stopifnot(sum(cov) == 1)

    U0 <- I0 <- A0 <- S0 <- T0 <- array(0, c(2, n_vax))
    # separate into 1:low and 2:high activity groups and by coverage
    N0 <- pars$N0 * outer(pars$q, cov)
    # set initial asymptomatic prevalence in each group (X AND H)
    A0[, 1] <- round(N0[, 1] * c(pars$prev_Asl, pars$prev_Ash))

    indextemp <- 2 * n_diag_rec + 3 * n_diag_rec * n_erlang + 1
    A0[, indextemp] <-
      round(N0[, indextemp] * c(pars$prev_Asl, pars$prev_Ash))

    # set initial uninfecteds
    U0 <- round(N0) - A0

    if (t > 0) {
      list(U0 = U0, I0 = I0, A0 = A0, S0 = S0, T0 = T0, t = t)
    } else {
      list(U0 = U0, I0 = I0, A0 = A0, S0 = S0, T0 = T0)
    }

  }

}



##' @name vax_params_xpvwrh
##' @title create vaccination parameters for use in onevax_xpvwrh model
##' @inheritParams vax_params_xvwr
##' @param hes proportion of population vaccine hesitant
##' @param r1 proportion of population offered vaccine only accepting the first
##' dose, becoming partially vaccinated
##' @param r2 proportion of the population who accepted the first dose of the
##' vaccine who go on to accept the second dose, becoming fully vaccinated
##' @param dur_v duration of time spent in V stratum after completing a round of
##' primary vaccination (fully vaccinated, accepting first and second dose)
##' @param dur_p duration of time spent in the P stratum, partially vaccinated
##' (accepting only the first dose)
##' @param vea_p scalar indicating efficacy of partial vaccination against
##'  acquisition (between 0-1)
##' @param vei_p scalar indicating efficacy of partial vaccination against
##'  infectiousness (between 0-1)
##' @param ved_p scalar indicating efficacy of partial vaccination against
##'  duration (between 0-1)
##' @param ves_p scalar indicating efficacy of partial vaccination against
##'  symptoms (between 0-1)
##' @param r2_p proportion of partially vaccinated individuals who receive
##' a second dose when returning to the clinic due to screening or illness
##' @param n_erlang integer giving the number of transitions that need to be
##' made through vaccine-protected strata until that protection has waned
##' @param n_diag_rec  number of diagnosis history strata
##' @param years_history number of years that diagnosis history is recorded for
##' @return A list parameters in the model input format
vax_params_xpvwrh <- function(vea = 0, vei = 0, ved = 0, ves = 0,
                              vea_revax = vea, vei_revax = vei, ved_revax = ved,
                              ves_revax = ves, vea_p = vea, vei_p = vei,
                              ved_p = ved, ves_p = ves, dur_v = 1e3,
                              dur_p = dur_v, dur_revax = dur_v, r1 = 0, r2 = 0,
                              r2_p = 0, booster_uptake = r1 * r2,
                              strategy = NULL, vbe = 0, t_stop = 99, hes = 0,
                              n_erlang = 1, n_diag_rec = 1, years_history = 1) {

  assert_scalar_unit_interval(vea_p)
  assert_scalar_unit_interval(vei_p)
  assert_scalar_unit_interval(ved_p)
  assert_scalar_unit_interval(ves_p)
  assert_scalar_unit_interval(vea)
  assert_scalar_unit_interval(vei)
  assert_scalar_unit_interval(ved)
  assert_scalar_unit_interval(ves)
  assert_scalar_unit_interval(vea_revax)
  assert_scalar_unit_interval(vei_revax)
  assert_scalar_unit_interval(ved_revax)
  assert_scalar_unit_interval(ves_revax)
  assert_scalar_positive(dur_v)
  assert_scalar_positive(dur_p)
  assert_scalar_positive(dur_revax)
  assert_scalar_unit_interval(r1)
  assert_scalar_unit_interval(r2)
  assert_scalar_unit_interval(r2_p)
  assert_scalar_unit_interval(booster_uptake)
  assert_scalar_unit_interval(vbe)
  assert_scalar_positive(t_stop)

  # waned partially-vaccinated individuals (P) move back to the non-vaccinated
  # stratum (X) and are considered immunologically naive. They are eligible
  # for another round of partial vaccination or full vaccination
  # waned fully-vaccinated individuals (V) move to their own waned stratum (W)
  # and are eligible for re-vaccination with a booster dose, moving into a
  # separate stratum, (R)
  # a proportion of all 'n' exist only in the hesitant compartment (H)
  # There is no movement between the willing (X, P, V, W, R) and hesitant (H)


  # generate indices of strata and calculate total # strata
  idx <- stratum_index_xpvwrh(n_erlang, n_diag_rec, strategy)

  # 1:X -> 3:V -> 4:W <-> 5:R
  # and
  # 1:X <-> 2:P

  #number of strata + sexual activity groups
  n_vax <- idx$n_vax
  n_group <- 2

  # work out strata eligibility for vaccination, vaccine protection,
  # waning, and where people wane to, based on n_vax calculated from
  # n_erlang

  i <- seq_len(idx$n_vax)

  n_group <- 2

  # create diagnosis history mapping
  diag_rec <- create_vax_map_branching(idx$n_vax, c(1, 1), idx$diagnosedfrom,
                                       idx$diagnosedto, set_vbe = FALSE, idx)


  # strata individuals wane from
  # i.e All strata under protection P, V and R

  i_v <- c(idx$P, idx$V, idx$R)

  # strata individuals wane to

  i_w <- c(idx$P[-(1:n_diag_rec)], idx$X, idx$V[-(1:n_diag_rec)], idx$W,
           idx$R[-(1:n_diag_rec)], idx$W)


  # ensure duration is not divided by 0
  ved <- min(ved, 1 - 1e-10)
  ved_revax <- min(ved_revax, 1 - 1e-10)

  # If uptake of VbE > 0 consider that all adolescents are offered vaccine
  p <- set_strategy(strategy, vbe > 0)

  vod <- create_vax_map_branching(n_vax, p$vod, idx$vaccinatedfrom_vod,
                                  idx$vaccinatedto_vod, set_vbe = FALSE, idx)

  vos <- create_vax_map_branching(n_vax, p$vos, idx$vaccinatedfrom_vos,
                                  idx$vaccinatedto_vos, set_vbe = FALSE, idx)

  vopn <- create_vax_map_branching(n_vax, p$vopn, idx$vaccinatedfrom_vopn,
                                   idx$vaccinatedto_vopn, set_vbe = FALSE, idx)


  vbe_map <- create_vax_map_branching(n_vax, p$vbe, idx$vaccinatedfrom_vbe,
                                      idx$vaccinatedto_vbe, set_vbe = TRUE, idx)

  if (sum(abs(vod)) > 0) {
    #vaccination on diagnosis occuring, so need to scale down diag_rec

    diag_rec[, idx$X, ] <- (1 - r1) * diag_rec[, idx$X, ]
    diag_rec[, idx$P, ] <- (1 - r2_p) * diag_rec[, idx$P, ]
    diag_rec[, idx$W, ] <- (1 - booster_uptake) * diag_rec[, idx$W, ]
  }

  u_s <- create_uptake_map_xpvwrh(vos, r1, r2, r2_p, booster_uptake, idx,
                                  n_diag_rec = n_diag_rec,
                                  screening_or_diagnosis = "screening")
  u_d <- create_uptake_map_xpvwrh(vod, r1, r2, r2_p, booster_uptake, idx,
                                  n_diag_rec = n_diag_rec,
                                  screening_or_diagnosis = "diagnosis")

  u_pn <- create_uptake_map_xpvwrh(vopn, r1, r2, r2_p, booster_uptake, idx,
                                   n_diag_rec = n_diag_rec,
                                   screening_or_diagnosis = "screening")

  w <- create_waning_map_branching(n_vax, i_v, i_w,
                                   n_erlang / c(dur_p, dur_v, dur_revax),
                                   n_erlang, n_diag_rec)


  willing <- c((1 - hes), rep(0, n_vax - 1 - n_diag_rec), hes,
               rep(0, n_diag_rec - 1))

  list(n_vax   = n_vax,
    willing = willing,
    u_d = u_d,
    u_s = u_s,
    u_pn = u_pn,
    u_vbe = vbe,
    vbe = vbe_map,
    vod = vod,
    vos = vos,
    vopn = vopn,
    vea = set_protection(i_v, idx, n_vax, vea_p, vea, vea_revax),
    vei = set_protection(i_v, idx, n_vax, vei_p, vei, vei_revax),
    ved = set_protection(i_v, idx, n_vax, ved_p, ved, ved_revax),
    ves = set_protection(i_v, idx, n_vax, ves_p, ves, ves_revax),
    w = create_waning_map_branching(n_vax, i_v, i_w,
                                    n_erlang / c(dur_p, dur_v, dur_revax),
                                    n_erlang, n_diag_rec),
    wd =  create_diagnosis_waning_map(n_vax, 1 / years_history, n_diag_rec),
    vax_t = c(0, t_stop),
    vax_y = c(1, 0),
    diag_rec = diag_rec,
    hesgroupmatrix = idx$hesgroupmatrix
  )
}



# ##' @name create_uptake_map_xpvwrh
# ##' @title Creates uptake mapping for the branching XPVWRH model where
# ##' individuals can move from unvaccinated (X) to vaccinated (V) or partially
# ##' vaccinated (P) as well as revaccinated from waned (W) to (R) and, and
# ##' partially vaccinated (P) to fully vaccianted (V). The former
# ##' reflects the specific indices which are chosen for assigning uptakes.
# ##' @param array a vaccine map array of dimensions n_group by n_vax by n_vax
# ##' generated through create_vax_map_branching()
# ##' @param r1 proportion of population offered vaccine only accepting the first
# ##' dose, becoming partially vaccinated
# ##' @param r2 proportion of the population who accepted the first dose of the
# ##' vaccine who go on to accept the second dose, becoming fully vaccinated
# ##' @param booster_uptake proportion of the formerly fully vaccinated, waned
# ##' population who accept a booster vaccination dose
# ##' @param r2_p proportion of partially vaccinated individuals who receive
# ##' a second dose when returning to the clinic due to screening or illness
# ##' @param idx list containing indices of all X, P, V, W, R & H strata and n_vax
# ##' through vaccine-protected strata until that protection has waned
# ##' @param n_diag_rec  number of diagnosis history strata
# ##' @return an array of the uptakes of same dimensions
# 
# create_uptake_map_xpvwrh <- function(array, r1, r2, r2_p, booster_uptake,
#                                      idx, n_diag_rec = 1, n_erlang = 1,
#                                      screening_or_diagnosis) {
# 
#   for (i in 1:n_diag_rec) {
# 
#     # note, these indices are specific to the branching pattern of xpvwrh
#     ## individuals in X accept vaccination of the 1st dose at an uptake of r1
# 
#     if (screening_or_diagnosis == "screening") {
#       temp <- i
#     } else if (screening_or_diagnosis == "diagnosis") {
# 
#       if (i < n_diag_rec) {
#         temp <- i + 1
#       } else {
#         temp <- i
#       }
#     } else {
#       print("uptake map type not specified.")
#     }
# 
#     array[, i, i] <- array[, i, i] * r1
# 
#     ## individuals entering V (fully vaccinated) also then accept
#     ## vaccination with the 2nd dose at an uptake of r2 so the proportion fully
#     ## vaccinated is given by r1 * r2
#     ## idx$V[1] gives index of the top of the V erlang stack
#     array[, idx$V[temp], i] <- array[, idx$V[temp], i] * (r1 * r2)
# 
#     ## individuals entering P (partially vaccinated) do not then accept
#     ## vaccination with the 2nd dose so proportion partially vaccinated is
#     ## given by r1 * (1 - r2), where 1 - r2 is the proportion not accepting the
#     ## 2nd dose given they have recieved the 1st dose
#     ## idx$P[1] gives index of the top of the P erlang stack
#     array[, idx$P[temp], i] <- array[, idx$P[temp], i] * (r1 * (1 - r2))
#   }
# 
# 
#   ## individuals with only the 1st dose can later accept vaccination with the
#   ## 2nd dose at an uptake of r2_p
#   ## idx$P gives indices for all P erlang strata, r2_p applies to all equally
#   array[, , idx$P] <- array[, , idx$P] * r2_p
# 
# 
#   ## individuals who were fully vaccinated and whose immunity has waned (W)
#   ## can accept vaccination with a single booster dose at an uptake of
#   ## booster_uptake
#   ## idx$W gives the the index for (W)
# 
#   array[, , idx$W] <- array[, , idx$W] * booster_uptake
# 
# 
#   # values must be positive - otherwise negative values in this array will
#   # cancel those in the vos and vod arrays = incorrect vaccination
#   abs(array)
# }


##' @name create_waning_map_branching
##' @title Create mapping for movement between strata due to vaccine waning
##' where waning from the partially vaccinated stratum (P) moves individuals
##' back to a naive unvaccinated state (X), and waning from fully vaccinated
##' stratum (V) moves individuals into a separate waned stratum (W)
##' Note, this structure is specific to xpvwrh
##' @param n_vax Integer denoting total number of strata
##' @param i_v indices of strata receiving protection through vaccination
##' @param i_w Scalar in (0, 6) denoting which stratum receives waned vaccinees
##' @param z Scalar denoting rate of waning
##' @param n_erlang integer giving the number of transitions that need to be
##' made through vaccine-protected strata until that protection has waned
##' @param n_diag_rec integer giving number of diagnosis history strata
##' @return an array of the mapping

create_waning_map_branching <- function(n_vax, i_v, i_w, z, n_erlang = 1,
                                        n_diag_rec = 1) {

  stopifnot(z > 0)

  #creates vector containing rates
  z_erlang <- rep(z, each = n_erlang * n_diag_rec)

  # set up waning map
  w <- array(0, dim = c(n_vax, n_vax))

  for (i in seq_along(i_v)) {
    w[i_v[i], i_v[i]] <- -z_erlang[i]
    w[i_w[i], i_v[i]] <- z_erlang[i]
  }

  w

}


##' @name create_vax_map_branching
##' @title Create mapping for movement between strata due to vaccination where
##' vaccination uptake splits off into two types (partial and full) from the
##' naive population (X). Different to create_vax_map as this function
##' specifically maps vbe to V (3) than P(2)
##' @param n_vax Integer denoting total number of strata
##' @param v 0-1 vector of length two indicating whether activity group
##'  should be offered vaccination.
##' @param i_e indices of strata eligible for vaccination
##' @param i_p indices of strata vaccinated and protected
##' @param set_vbe Boolean which indicates that vaccination is occurring at some
##' level of uptake upon entering the model
##' @param idx list containing indices of all X, P, V, W, R & H strata and n_vax
##' through vaccine-protected strata until that protection has waned
##' @return an array of the mapping

create_vax_map_branching <- function(n_vax, v, i_e, i_p, set_vbe = FALSE, idx) {

  # ensure vaccine input is of correct length
  n_group <- 2
  n_vax <- idx$n_vax

  stopifnot(length(v) == n_group)
  stopifnot(all(v %in% c(0, 1)))


  if (length(i_e) > 0) {
    stopifnot(max(i_e, i_p) <= n_vax)
  }

  # set up vaccination matrix
  vax_map <- array(0, dim = c(n_group, n_vax, n_vax))

  if (set_vbe == TRUE) {

    for (i in 1:length(idx$X)){
      vax_map[, idx$X[i], idx$X[i]] <-  v
      vax_map[, idx$V[i], idx$X[i]] <- -v
    }

  } else {

    #repeat over stratum 1 column 1 for ease
    for (i in seq_along(i_e)) {

      vax_map[, i_e[i], i_e[i]] <-  v
      vax_map[, i_p[i], i_e[i]] <- -v

    }
  }

  vax_map
}

##' @name run_onevax_xpvwrh
##' @title run model with a two-dose vaccine for input parameter sets, either
##' from initialisation or from equilibrium, those with waned vaccines are
##' eligible for revaccination (R), and return to the R stratum, those with
##' waned partial vaccines return to the unvaccinated stratum (X) and considered
##' immunologically naive. A user defined proportion of the population is
##' vaccine hesitant and is never vaccinated. Full acciantion (V) with 2 doses
##' gives maximum protection whereas partial vaccination with 1 dose (P) gives
##' less. Individuals can get 2 doses either by committing in X or upon visiting
##' a clinic in P.
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
##' @param dur_v duration of time spent in V stratum after completing a round of
##' primary vaccination (fully vaccinated, accepting first and second dose)
##' @param dur_p duration of time spent in the P stratum, partially vaccinated
##' (accepting only the first dose)
##' @param hes Proportion of individuals in the population who are vaccine
##'  hesitant
##' @param vea_p scalar indicating efficacy of partial vaccination against
##'  acquisition (between 0-1)
##' @param vei_p scalar indicating efficacy of partial vaccination against
##'  infectiousness (between 0-1)
##' @param ved_p scalar indicating efficacy of partial vaccination against
##'  duration (between 0-1)
##' @param ves_p scalar indicating efficacy of partial vaccination against
##'  symptoms (between 0-1)
##' @param r1 scalar or numeric vector with same length as
##'  'gono_params' giving proportion of population offered vaccine only
##'   accepting the first dose, becoming partially vaccinated
##' @param r2 scalar or numeric vector with same length as
##'  'gono_params' giving proportion of the population who accepted the first
##'   dose of the vaccine who go on to accept the second dose, becoming fully
##'   vaccinated
##' @param r2_p scalar or numeric vector with same length as 'gono_params'
##' giving proportion of partially vaccinated individuals who later receive
##' a second dose when returning to the clinic due to screening or illness
##' @param booster_uptake scalar or numeric vector with same length as
##'  'gono_params' giving proportion of population undertaking booster
##'  vaccination after primary vaccination protection has waned.
##'   Defaults to supplied value of r1 * r2
##' @param n_erlang integer giving the number of erlang vaccination transitions
##'  through vaccine-protected strata until that protection has waned
##' @param n_diag_rec integer for the number of diagnosis history substrata
##' @param years_history number of years that diagnosis history is recorded for
##' @inheritParams run_onevax_xvwv
##' @return A list of transformed model outputs
##' @export
run_onevax_xpvwrh <- function(tt, gono_params, init_params = NULL,
                              dur_v = 1e3, dur_p = dur_v, vea = 0, vei = 0,
                              ved = 0, ves = 0, dur_revax = dur_v,
                              vea_revax = vea, vei_revax = vei, ved_revax = ved,
                              ves_revax = ves, vea_p = vea, vei_p = vei,
                              ved_p = ved, ves_p = ves, vbe = 0, r1 = 0, r2 = 0,
                              r2_p = 0, booster_uptake = (r1 * r2),
                              strategy = NULL, t_stop = 99, hes = 0,
                              n_erlang = 1, n_diag_rec = 1, years_history = 1) {

  stopifnot(all(lengths(list(booster_uptake, r1, r2, r2_p, vea, vei, ved, ves,
                             vea_revax, vei_revax, ved_revax, vea_p, vei_p,
                             ves_p, ved_p, dur_v, dur_p, ves_revax,
                             dur_revax)) %in% c(1, length(gono_params))))

  vax_params <- Map(vax_params_xpvwrh, r1 = r1, r2 = r2, r2_p = r2_p,
                    booster_uptake = booster_uptake, dur_v = dur_v,
                    vea = vea, vei = vei, ved = ved, ves = ves,
                    vea_p = vea_p, vei_p = vei_p, ved_p = ved_p, ves_p = ves_p,
                    dur_revax = dur_revax, dur_p = dur_p,
                    vea_revax = vea_revax, vei_revax = vei_revax,
                    ved_revax = ved_revax, ves_revax = ves_revax, hes = hes,
                    n_erlang = n_erlang, n_diag_rec = n_diag_rec,
                    years_history = years_history,
                    MoreArgs = list(strategy = strategy,
                                    t_stop = t_stop, vbe = vbe))

  if (is.null(init_params)) {
    pars <- lapply(gono_params, model_params)

    init_params <- lapply(pars, initial_params_xpvwrh, hes = hes,
                          n_erlang = n_erlang, n_diag_rec = n_diag_rec)
  } else {

    #check if init_params supplied, n_vax corresponds to the n_erlang
    # supplied to the run function
    stopifnot(length(init_params[[1]][[1]]) / 2 ==
                3 * n_diag_rec + (3 * n_diag_rec * n_erlang))
  }

  ret <- Map(run_xpvwrh, gono_params = gono_params, vax_params = vax_params,
             init_params = init_params, n_erlang = n_erlang,
             n_diag_rec = n_diag_rec, MoreArgs = list(tt = tt))

  # name outputs
  ret <- lapply(ret, name_outputs, gen_erlang_labels(n_erlang, n_diag_rec))
  ret
}


##' @name gen_erlang_labels
##' @title generates the appropriate strata labels for the number of strata
##' in the model, which depends on the value given to n_erlang
##' @param n_erlang integer giving the number of transitions that need to be
##'  made through vaccine-protected strata until that protection has waned
##' @param n_diag_rec integer for the number of diagnosis history substrata
##' @return a character vector of length n_vax containing strata labels
##' @export
gen_erlang_labels <- function(n_erlang = 1, n_diag_rec = 1) {

  idx <- seq_len(n_erlang)
  idy <- seq_len(n_diag_rec)
  output <- c()

  stratavec <- c("X", "P", "V", "W", "R", "H")
  diag_hist <- paste0(".", as.roman(seq_len(n_diag_rec)))

  output <- c(paste0("X", diag_hist),
              paste0("P", rep(seq_len(n_erlang), each = n_diag_rec), diag_hist),
              paste0("V", rep(seq_len(n_erlang), each = n_diag_rec), diag_hist),
              paste0("W", diag_hist),
              paste0("R", rep(seq_len(n_erlang), each = n_diag_rec), diag_hist),
              paste0("H", diag_hist))

  return(output)
}

##' @name set_protection
##' @title generates vector which tells the model which strata are under
##' vaccine protection and what the value of protection for that strata is
##' @param i_v indices of strata receiving protection through vaccination
##' @param idx list containing indices of all X, P, V, W, R & H strata and n_vax
##' through vaccine-protected strata until that protection has waned
##' through vaccine-protected strata until that protection has waned
##' @param n_vax integer denoting total number of strata
##' @param ve_p scalar 0-1 with degree of partial primary protection of the P(N)
##' strata, can take vea_p, vei_p, ves_p, ved_p
##' @param ve scalar 0-1 with degree of full primary protection of the V(N)
##' strata, can take vea, vei, ves, ved
##' @param ve_revax scalar 0-1 with degree of re-vaccinated protection of the
##' R(N) strata, can take vea_revax, vei_revax, ves_revax, ved_revax
##' @return vector of length n_vax with zeros corresponding to the indices of
##' strata with no protection, and the supplied degree of partial, full, and
##' boosted protection corresponding to the indices of strata with partial,
##' full and boosted vaccination status
##' @export
set_protection <- function(i_v, idx, n_vax, ve_p, ve, ve_revax) {

  # get indexes of strata under protection by type of protection
  p <- idx$P
  v <- idx$V
  r <- idx$R

  # generate empty vector as long as n_vax
  ve_vec <- rep(0, idx$n_vax)

  # assign corresponding level of protection to the correct position
  ve_vec[p] <- ve_p
  ve_vec[v] <- ve
  ve_vec[r] <- ve_revax

  ve_vec
}


##' @name restart_hes
##' @title uses XPVWRH model run in the absence of vaccination or hesitancy.
##' Saves down the number of individuals in each compartment, and moves
##' a given proportion (hes) of them from the X to the H strata to generate
##' new initial conditions in the presence of hesitancy.
##' @inheritParams restart_params
##' @param hes proportion of population vaccine hesitant
##' @param branching boolean to denote if xpvwrh branching model in use
##' @param n_erlang integer giving the number of transitions that need to be
##'  made
##' @param n_diag_rec integer for the number of diagnosis history substrata
##' @return A list of initial conditions to restart a model with n_vax
##' vaccination levels, and a populated hestitant stratum in the given
##' proportion 'hes'
##' @export

restart_hes <- function(y, n_vax = 6, hes = 0, n_erlang = 1, n_diag_rec = 1,
                        branching = FALSE) {

  dim_y <- dim(y[["U"]])

  i_t <- dim_y[1] # number of timepoints
  n_vax_input <- dim_y[3] # number is strata
  i_vax <- seq_len(n_vax_input)

  for (d  in 1:n_diag_rec) {
    if (round(rowSums(y$N[, , n_vax_input - n_diag_rec + d])[dim_y[1]],
              5) > 0) {
      stop("Provided model run already contains hesitancy > 0")
    }

    if (round(rowSums(y$N[, , n_diag_rec + d])[dim_y[1]], 5) > 0) {
      stop("Provided model run has vaccination, baseline run should have all V
          = 0")
    }

    # branched xpvwrh models have 2 primary vaccination compartments to check
    if (branching == TRUE) {
      if (round(rowSums(y$N[, , n_diag_rec +
                              (n_diag_rec * n_erlang) + d])[dim_y[1]], 5) > 0) {
        stop("Provided model run has vaccination, baseline run should have all V
          = 0")
      }
    }
  }

  #create blank array, 2activity groups by number of strata
  U0 <- I0 <- A0 <- S0 <- T0 <- array(0, c(2, n_vax_input))

  # set compartments new initial conditions based on final position of y
  U0[, i_vax] <- y$U[i_t, , i_vax]
  I0[, i_vax] <- y$I[i_t, , i_vax]
  A0[, i_vax] <- y$A[i_t, , i_vax]
  S0[, i_vax] <- y$S[i_t, , i_vax]
  T0[, i_vax] <- y$T[i_t, , i_vax]

  # move correct number from equilibrium X to H

  for (d in 1:n_diag_rec) {
    U0[, n_vax_input - n_diag_rec + d] <- h <- U0[, d] * hes
    U0[, d] <- U0[, d] - h

    I0[, n_vax_input - n_diag_rec + d] <- h <- I0[, d] * hes
    I0[, d] <- I0[, d] - h

    A0[, n_vax_input - n_diag_rec + d] <- h <- A0[, d] * hes
    A0[, d] <- A0[, d] - h

    S0[, n_vax_input - n_diag_rec + d] <- h <- S0[, d] * hes
    S0[, d] <- S0[,  d] - h

    T0[, n_vax_input - n_diag_rec + d] <- h <- T0[, d] * hes
    T0[,  d] <- T0[,  d] - h
  }

  list(U0 = U0, I0 = I0, A0 = A0, S0 = S0, T0 = T0, t = y$t[i_t])
}
