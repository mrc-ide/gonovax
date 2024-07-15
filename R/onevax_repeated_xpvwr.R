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
##' @name initial_params_repeated_xpvwr
##' @title Initial conditions for the model
##' @inheritParams initial_params
##' @param hesgroups number of different vaccine hesitancy groups 
##' (with different vaccine sentiments)
##' @param t number of years, only use when using function outside of
##' run_onevax_xpvwrh() to generate initial conditions for tests
##' @param coverage_x now ARRAY of 'coverage' of non-vaccination
##' @param coverage_p now ARRAY of partial (one-dose) vaccine coverage of 
##' the population already present (as a proportion)
##' @param coverage_v now ARRAY of two-dose vaccine coverage of 
##' the population already present (as a proportion)
##' @param n_erlang integer giving the number of transitions that need to be
##'  made through vaccine-protected strata until that protection has waned
##' @param n_diag_rec  number of diagnosis history strata
##' @return A list of initial conditions
##' @export
initial_params_repeated_xpvwr <- function(pars, hesgroups = 1,
                                          coverage_x = rep(1, hesgroups)/hesgroups, 
                                          coverage_p = rep(0, hesgroups),
                                          coverage_v = rep(0, hesgroups), 
                                          t = FALSE,
                                          n_erlang = 1, n_diag_rec = 1) {

  idx <- stratum_index_repeated_xpvwr(n_erlang, n_diag_rec)
  
  stopifnot(length(coverage_x) == hesgroups)
  stopifnot(length(coverage_p) == hesgroups)
  stopifnot(length(coverage_v) == hesgroups)
  
  #print(sum(coverage_p + coverage_v))
  
  if (sum(coverage_p + coverage_v) >= 1) {
    stop("sum of coverages must not exceed or equal 1")
  } else {

    n_vax <- idx$n_vax
    x_init <- (coverage_x - coverage_p - coverage_v)
    p_init <- coverage_p
    v_init <- coverage_v

    ## X[1], P[n_erlang], V[n_erlang], W[1], R[n_erlang], H[1]
    
    
    cov <- c( c(rbind(x_init, matrix(0, nrow = n_diag_rec - 1, ncol = length(x_init)))),
              c(rbind(p_init, matrix(0, nrow = n_diag_rec - 1, ncol = length(p_init)))),
              c(rbind(v_init, matrix(0, nrow = n_diag_rec - 1, ncol = length(v_init)))),
              rep(0, n_diag_rec * hesgroups),
              rep(0, n_diag_rec * n_erlang * hesgroups)
              )

    stopifnot(length(cov) == n_vax * hesgroups)
    stopifnot(sum(cov) == 1)

    U0 <- I0 <- A0 <- S0 <- T0 <- array(0, c(2, n_vax * hesgroups))
    # separate into 1:low and 2:high activity groups and by coverage
    N0 <- pars$N0 * outer(pars$q, cov)
    
    # set initial asymptomatic prevalence in each X group 
    
    A0[, seq_len(hesgroups)] <- round(N0[, seq_len(hesgroups) ] * c(pars$prev_Asl, pars$prev_Ash))
    
   # print("round N0")
  #  print(round(N0))
    
  #  print("A0")
  #  print(A0)
    
    # set initial uninfecteds
    U0 <- round(N0) - A0

    if (t > 0) {
      list(U0 = U0, I0 = I0, A0 = A0, S0 = S0, T0 = T0, t = t)
    } else {
      list(U0 = U0, I0 = I0, A0 = A0, S0 = S0, T0 = T0)
    }

  }

}



##' @name vax_params_repeated_xpvwr
##' @title create vaccination parameters for use in onevax_xpvwrh model
##' @inheritParams vax_params_xvwr
##' @param hesgroups proportion of population vaccine hesitant
##' @param r1 array of proportion of population offered vaccine only
##' accepting the first dose, becoming partially vaccinated
##' @param r2 array of proportion of the population who accepted the first dose 
##' of the vaccine who go on to accept the second dose,
##' becoming fully vaccinated
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
##' @param r2_p vector of proportion of partially vaccinated individuals who 
##' receive a second dose when returning to the clinic due to screening 
##' or illness
##' @param n_erlang integer giving the number of transitions that need to be
##' made through vaccine-protected strata until that protection has waned
##' @param n_diag_rec  number of diagnosis history strata
##' @param years_history number of years that diagnosis history is recorded for
##' @return A list parameters in the model input format
vax_params_repeated_xpvwr <- function(vea = 0, vei = 0, ved = 0, ves = 0,
                              vea_revax = vea, vei_revax = vei, ved_revax = ved,
                              ves_revax = ves, vea_p = vea, vei_p = vei,
                              ved_p = ved, ves_p = ves, dur_v = 1e3,
                              dur_p = dur_v, dur_revax = dur_v, r1 = 0, r2 = 0,
                              r2_p = 0, booster_uptake = mapply(function(x, y) x * y, r1, r2, SIMPLIFY = FALSE),
                              strategy = NULL, vbe = 0, t_stop = 99,
                              hesgroups = 1,
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
 # assert_scalar_unit_interval(r1)
#  assert_scalar_unit_interval(r2)
 # assert_scalar_unit_interval(r2_p)
  #assert_scalar_unit_interval(booster_uptake)
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
  idx <- stratum_index_repeated_xpvwr(n_erlang, n_diag_rec, hesgroups = hesgroups, strategy)
  
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

  
  # for an n_erlang of 3, and n_diag_rec of 2, the list of indexes returned
  # will be in the following order, where roman numerals refer to n_diag_rec,
  # and arabic numerals refer to erlang, and lower case letters refer to 
  # hesitancy group
  # X.I.a, X.II.a, X.I.b, X.II.b, 
  # P1.I.a, P1.II.a, P1.I.b, P1.II.b,
  # P2.I.a, P2.II.a, P2.I.b, P2.II.b, 
  # P3.I.a, P3.II.a, P3.I.b, P3.II.b,
  # V1.I.a, V1.II.a, V1.I.b, V1.II.b, 
  # V2.I.a, V2.II.a, V2.I.b, V2.II.b,
  # V3.I.a, V3.II.a, V3.I.b, V3.II.b,
  # W.I.a, W.II.a, W.I.b, W.II.b,
  # R1.I.a, R1.II.a, R1.I.b, R1.II.b,
  # R2.I.a, R2.II.a, R1.I.b, R1.II.b,
  # R3.I.a, R3.II.a, R3.I.b, R3.II.b
  
  

  # strata individuals wane from
  # i.e All strata under protection P, V and R

  i_v <- c(idx$P, idx$V, idx$R)

  # strata individuals wane to

  i_w <- c(idx$P[-(seq_len(n_diag_rec * hesgroups))], idx$X, idx$V[-(seq_len(n_diag_rec * hesgroups))], idx$W,
           idx$R[-(seq_len(n_diag_rec * hesgroups))], idx$W)


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

    
    
    ## trystan note - this will work when n_diag_rec = 1 and n_erlang = 1, but
    ## will probably break if either > 1
    
    diag_rec[, idx$X, ] <- (1 - r1) * diag_rec[, idx$X, ]
    diag_rec[, idx$P, ] <- (1 - r2_p) * diag_rec[, idx$P, ]
    diag_rec[, idx$W, ] <- (1 - booster_uptake) * diag_rec[, idx$W, ]
  }
  
  print("r2p")
  print(r2_p)
  print("booster uptake")
  print(booster_uptake)

  u_s <- create_uptake_map_repeated_xpvwr(vos, r1, r2, r2_p, booster_uptake, idx,
                                  n_diag_rec = n_diag_rec, hesgroups = hesgroups,
                                  screening_or_diagnosis = "screening")
  u_d <- create_uptake_map_repeated_xpvwr(vod, r1, r2, r2_p, booster_uptake, idx,
                                  n_diag_rec = n_diag_rec, hesgroups = hesgroups,
                                  screening_or_diagnosis = "diagnosis")

  u_pn <- create_uptake_map_repeated_xpvwr(vopn, r1, r2, r2_p, booster_uptake, idx,
                                   n_diag_rec = n_diag_rec, hesgroups = hesgroups,
                                   screening_or_diagnosis = "screening")

  w <- create_waning_map_branching_repeat(n_vax, i_v, i_w,
                                   rep(n_erlang / c(dur_p, dur_v, dur_revax), each = hesgroups),
                                   n_erlang, n_diag_rec)

  
  willing = c(rep(c(1, rep(0, n_diag_rec - 1)), hesgroups), 
              rep(0, n_vax - n_diag_rec * hesgroups))
  

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
    w = w,
    wd =  create_diagnosis_waning_map(n_vax, 1 / years_history, n_diag_rec),
    vax_t = c(0, t_stop),
    vax_y = c(1, 0),
    diag_rec = diag_rec
  )
}



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
                                        n_diag_rec = 1, hesgroups = 1) {

  stopifnot(z > 0)

  #creates vector containing rates
  z_erlang <- rep(z, each = n_erlang * n_diag_rec * hesgroups)

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

    vax_map[, idx$X[1], idx$X[1]] <-  v
    vax_map[, idx$V[1], idx$X[1]] <- -v

  } else {

    #repeat over stratum 1 column 1 for ease
    for (i in seq_along(i_e)) {

      vax_map[, i_e[i], i_e[i]] <-  v
      vax_map[, i_p[i], i_e[i]] <- -v

    }
  }

  vax_map
}

##' @name run_onevax_repeated_xpvwr
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
##' @param r1 list of length 1 or of length 'gono_params' with each entry 
##'   itself a list of proportions of length hesgroups only accepting the first
##'   dose, becoming partially vaccinated
##' @param r2 list of length 1 or of length 'gono_params' with each entry 
##'   itself a list of proportions of length hesgroups of the population who 
##'   accepted the first dose of the vaccine who go on to accept the second 
##'   dose, becoming fully vaccinated
##' @param r2_p list of length 1 or of length 'gono_params' with each entry 
##'   itself a list of proportions of partially vaccinated individuals who 
##'   later receive a second dose when returning to the clinic due to screening 
##'   or illness
##' @param booster_uptake list of length 1 or of length 'gono_params' with each 
##'   entry itself a list of proportions of length hesgroups undertaking booster
##'   vaccination after primary vaccination protection has waned.
##'   Defaults to supplied value of r1 * r2
##' @param n_erlang integer giving the number of erlang vaccination transitions
##'  through vaccine-protected strata until that protection has waned
##' @param n_diag_rec integer for the number of diagnosis history substrata
##' @param years_history number of years that diagnosis history is recorded for
##' @inheritParams run_onevax_xvwv
##' @return A list of transformed model outputs
##' @export
run_onevax_repeated_xpvwr <- function(tt, gono_params, init_params = NULL,
                              dur_v = 1e3, dur_p = dur_v, vea = 0, vei = 0,
                              ved = 0, ves = 0, dur_revax = dur_v,
                              vea_revax = vea, vei_revax = vei, ved_revax = ved,
                              ves_revax = ves, vea_p = vea, vei_p = vei, hesgroups = 1,
                              ved_p = ved, ves_p = ves, vbe = 0, r1 = list(rep(0, hesgroups)), r2 = list(rep(0, hesgroups)),
                              r2_p = list(rep(0, hesgroups)), booster_uptake = mapply(function(x, y) x * y, r1, r2, SIMPLIFY = FALSE),
                              strategy = NULL, t_stop = 99, hes = 0,
                              n_erlang = 1, n_diag_rec = 1, years_history = 1) {
  
  
  print("hello")

  stopifnot(all(lengths(list(vea, vei, ved, ves,
                             vea_revax, vei_revax, ved_revax, vea_p, vei_p,
                             ves_p, ved_p, dur_v, dur_p, ves_revax,
                             dur_revax)) %in% c(1, length(gono_params))))
  
  
  #stopifnot(all(lengths(list(booster_uptake, r1, r2, r2_p)) == hesgroups))
  
  stopifnot(all(lengths(list(booster_uptake, r1, r2, r2_p)) %in% c(1, length(gono_params))))
  
  print(r1)
  
  print("hello Allo")
  print(lengths(list(booster_uptake, r1, r2, r2_p)))
  
 
  
  stopifnot(all(sapply(list(r1, r2, booster_uptake, r2_p), function(x) all(lengths(x) == hesgroups))))
  
  print("hello B")
  print(sapply(list(r1, r2, booster_uptake, r2_p), function(x) lengths(x)))
  

  
  vax_params <- Map(vax_params_repeated_xpvwr, r1 = r1, r2 = r2, r2_p = r2_p,
                    booster_uptake = booster_uptake, dur_v = dur_v,
                    vea = vea, vei = vei, ved = ved, ves = ves,
                    vea_p = vea_p, vei_p = vei_p, ved_p = ved_p, ves_p = ves_p,
                    dur_revax = dur_revax, dur_p = dur_p,
                    vea_revax = vea_revax, vei_revax = vei_revax,
                    ved_revax = ved_revax, ves_revax = ves_revax,
                    n_erlang = n_erlang, n_diag_rec = n_diag_rec,
                    years_history = years_history, hesgroups = hesgroups,
                    MoreArgs = list(strategy = strategy,
                                    t_stop = t_stop, vbe = vbe))
  
  print("hello2")

  if (is.null(init_params)) {
    
    
    print("hello3")
    pars <- lapply(gono_params, model_params)

    init_params <- lapply(pars, initial_params_repeated_xpvwr, hesgroups = hesgroups,
                          n_erlang = n_erlang, n_diag_rec = n_diag_rec)
  } else {

    print("hello4")
    #check if init_params supplied, n_vax corresponds to the n_erlang
    # supplied to the run function
    stopifnot(length(init_params[[1]][[1]]) / 2 ==
                hesgroups * (2 * n_diag_rec + (3 * n_diag_rec * n_erlang)))
  }

  print("hello5")
  
  ret <- Map(run_repeated_xpvwr, gono_params = gono_params, vax_params = vax_params,
             init_params = init_params, n_erlang = n_erlang, hesgroups = hesgroups,
             n_diag_rec = n_diag_rec, MoreArgs = list(tt = tt))

  # name outputs
  ret <- lapply(ret, name_outputs, gen_erlang_labels_repeat(n_erlang, n_diag_rec, hesgroups))
  ret
}


##' @name gen_erlang_labels_repeat
##' @title generates the appropriate strata labels for the number of strata
##' in the model, which depends on the value given to n_erlang
##' @param n_erlang integer giving the number of transitions that need to be
##'  made through vaccine-protected strata until that protection has waned
##' @param n_diag_rec integer for the number of diagnosis history substrata
##' @return a character vector of length n_vax containing strata labels
##' @export
gen_erlang_labels_repeat <- function(n_erlang = 1, n_diag_rec = 1, hesgroups = 1) {

  idx <- seq_len(n_erlang)
  idy <- seq_len(n_diag_rec)
  output <- c()

  stratavec <- c("X", "P", "V", "W", "R")
  diag_hist <- paste0(".", as.roman(seq_len(n_diag_rec)))
  letters_seq <- paste0(".", letters[seq_len(hesgroups)])
  

  output <- c(output, paste0("X", diag_hist, ".", rep(letters_seq, each = n_diag_rec)))
  for (i in seq_len(n_erlang)) {
    output <- c(output, paste0("P", i, diag_hist, ".", rep(letters_seq, each = n_diag_rec)))
  }
  for (i in seq_len(n_erlang)) {
    output <- c(output, paste0("V", i, diag_hist, ".", rep(letters_seq, each = n_diag_rec)))
  }
  output <- c(output, paste0("W", diag_hist, ".", rep(letters_seq, each = n_diag_rec)))
  for (i in seq_len(n_erlang)) {
    output <- c(output, paste0("R", i, diag_hist, ".", rep(letters_seq, each = n_diag_rec)))
  }

  return(output)
}


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

create_waning_map_branching_repeat <- function(n_vax, i_v, i_w, z, n_erlang = 1,
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

