### gonovax population model whereby individuals in X become vaccinated
# under the given vacc strategy with proportion r1 accepting 1 dose --> P
# or proportion r1r2 accepting both doses --> V
# individuals in P wane back to X, and can then go to P or V again after
# individuals in V wane to W where they can get a booster vacc to R
# rate of waning for P is 18 months
# this is a modificaiton of the onevax_xvwrh model where indiviuals pass
# linearly through the strata, now there is a branching after X

# also note to self may need to edit 'restart_hes' function in the onevax_xvwrh
# code so it works for this as well! rather than creating whole new function

##' Create initial conditions for the model
##' @name initial_params_xpvwrh
##' @title Initial conditions for the model
##' @inheritParams initial_params
##' @param hes proportion of population vaccine hesitant
##' @return A list of initial conditions
##' @export
initial_params_xpvwrh <- function(pars, coverage_p = 0, coverage_v = 0,
                                  hes = 0) {
  assert_scalar_unit_interval(coverage_p)
  assert_scalar_unit_interval(coverage_v)
  n_vax <- 6
  willing <- 1 - hes
  x_init <- willing * (1 - coverage_p - coverage_v)
  p_init <- willing * coverage_p
  v_init <- willing * coverage_v
  cov <- c(x_init, p_init, v_init, 0, 0, hes)

  stopifnot(length(cov) == n_vax)
  stopifnot(sum(cov) == 1)

    U0 <- I0 <- A0 <- S0 <- T0 <- array(0, c(2, n_vax))
  # separate into 1:low and 2:high activity groups and by coverage
  N0 <- pars$N0 * outer(pars$q, cov)
  # set initial asymptomatic prevalence in each group (X AND H)
  A0[, 1] <- round(N0[, 1] * c(pars$prev_Asl, pars$prev_Ash))
  A0[, n_vax] <- round(N0[, n_vax] * c(pars$prev_Asl, pars$prev_Ash))

  # set initial uninfecteds
  U0 <- round(N0) - A0

  list(U0 = U0, I0 = I0, A0 = A0, S0 = S0, T0 = T0)

}


##' @name vax_params_xpvwrh
##' @title create vaccination parameters for use in onevax_xpvwrh model
##' @inheritParams vax_params_xvwr
##' @param hes proportion of population vaccine hesitant
##' @param r1 proportion of population offered vaccine only accepting the first
##' dose
##' @param r1r2 proportion of population offered vaccine who accept both the
##' first and second dose
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
##' @return A list parameters in the model input format
vax_params_xpvwrh <- function(vea = 0, vei = 0, ved = 0, ves = 0,
                             vea_revax = vea, vei_revax = vei,
                             ved_revax = ved, ves_revax = ves,
                             vea_p = vea, vei_p = vei, ved_p = ved, ves_p = ves,
                             dur_v = 1e3, dur_p = dur_v, dur_revax = dur_v,
                             r1 = 0, r1r2 = 0,
                             booster_uptake = r1r2, strategy = NULL,
                             vbe = 0, t_stop = 99, hes = 0) {

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
  assert_scalar_unit_interval(r1r2)
  assert_scalar_unit_interval(booster_uptake)
  assert_scalar_unit_interval(vbe)
  assert_scalar_positive(t_stop)

  # waned partially-vaccinated individuals (P) move back to the non-vaccinated
  # stratum (X) and are considered immunologically naive. They are eligible
  # for another round of partial vaccination or full vaccination
  # waned fully-vaccinated individuals (V) move to their own waned stratum (W)
  # and are eligible for re-vaccination with a booster, moving into a separate
  # stratum (R)
  # a proportion of all 'n' exist only in the hesitant compartment (H)
  # there is no movement between the willing (X, P, V, W, R) and hesitant (H)

  # 1:X -> 3:V -> 4:W <-> 5:R
  # and
  # 1:X <-> 2:P

  i_eligible <- c(1, 4)                # X(1) and W(4) eligible for vaccination
  i_w <- c(1, 4, 4)                   # Waned vaccinees move to X(1) and W(4)
  i_v <- c(2, 3, 5)                    # P(2), V(3), and R(5) are protected

  #number of compartments
  n_vax <- 6
  n_group <- 2

  # ensure duration is not divided by 0
  ved <- min(ved, 1 - 1e-10)
  ved_revax <- min(ved_revax, 1 - 1e-10)

  # If uptake of VbE > 0 consider that all adolescents are offered vaccine
  p <- set_strategy(strategy, vbe > 0)

  # set up uptake matrix rows = groups, columns = vaccine strata
  # tells us where vaccinated individuals are going to be pulled from

  u <- matrix(0, n_group, n_vax)
  u[, i_eligible[1]] <- r1r2 + r1
  u[, i_eligible[2]] <- booster_uptake


  list(n_vax   = n_vax,
       willing = c((1 - hes), 0, 0,  0, 0, hes),
       u       = u,
       u_vbe   = vbe,
       vbe     = create_vax_map_branching(n_vax, p$vbe, i_eligible, i_v,
                                          set_vbe = TRUE),
       vod     = create_vax_map_branching(n_vax, p$vod, i_eligible, i_v,
                                          r1 = r1, r1r2 = r1r2),
       vos     = create_vax_map_branching(n_vax, p$vos, i_eligible, i_v,
                                          r1 = r1, r1r2 = r1r2),
       vea     = c(0, vea_p, vea, 0, vea_revax, 0),
       vei     = c(0, vei_p, vei, 0, vei_revax, 0),
       ved     = c(0, ved_p, ved, 0, ved_revax, 0),
       ves     = c(0, ves_p, ves, 0, ves_revax, 0),
       w       = create_waning_map_branching(n_vax, i_v,
                                   i_w, 1 / c(dur_p, dur_v, dur_revax)),
       vax_t   = c(0, t_stop),
       vax_y   = c(1, 0)
  )
}

##' @name create_waning_map_branching
##' @title Create mapping for movement between strata due to vaccine waning
##' where waning from the partially vaccinated stratum (P) moves individuals
##' back to a naive unvaccinated state (X), and waning from fully vaccinated
##' stratum (V) moves individuals into a separate waned stratum (W)
##' @param n_vax Integer in (0, 6) denoting total number of strata
##' @param i_v indices of strata receving protection through vaccination
##' @param i_w Scalar in (0, 6) denoting which stratum receives waned vaccinees
##' @param z Scalar denoting rate of waning
##' @return an array of the mapping

create_waning_map_branching <- function(n_vax, i_v, i_w, z) {

  stopifnot(z > 0)
  stopifnot(length(z) %in% c(1, length(i_v)))
  stopifnot(length(i_w) == 3)

  # set up waning map
  w <- array(0, dim = c(n_vax, n_vax))

  loop <- length(i_v)

  for (i in 1:loop) {
    w[i_w[i], i_v[i]] <- z[i]
    w[i_v[i], i_v[i]] <- -z[i]
  }

  w

}


##' @name create_vax_map_branching
##' @title Create mapping for movement between strata due to vaccination where
##' vaccination uptake splits off into two types (partial and full) from the
##' naive population (X)
##' @param n_vax Integer denoting total number of strata
##' @param v 0-1 vector of length two indicating whether activity group
##'  should be offered vaccination.
##' @param i_u indices of strata eligible for vaccination
##' @param i_v indices of strata vaccinated and protected
##' @return an array of the mapping

create_vax_map_branching <- function(n_vax, v, i_u, i_v, r1 = NULL,
                                     r1r2 = NULL, set_vbe = FALSE) {

  # ensure vaccine input is of correct length
  n_group <- 2
  stopifnot(length(v) == n_group)
  stopifnot(all(v %in% c(0, 1)))
  stopifnot(max(i_u, i_v) <= n_vax)

  # set up vaccination matrix
  vax_map <- array(0, dim = c(n_group, n_vax, n_vax))

  if(set_vbe == TRUE){
  
  vax_map[, 1, 1] <-  v
  vax_map[, 3, 1] <- -v 
    
  }

  # tweak eligibility , repeat over stratum 1 column 1 for ease
  i_u <- c(1, i_u)

  for (i in seq_along(i_u)) {
    vax_map[, i_u[i], i_u[i]] <-  v
    vax_map[, i_v[i], i_u[i]] <- -v

  }

  #obtain proportions
 
  tot <- r1r2 + r1
  
  ifelse(r1r2 > 0, prop_full <- r1r2/tot, prop_full <- 0)
  ifelse(r1 > 0, prop_part <- r1/tot, prop_part <- 0)

  vax_map[, 2, 1] <- vax_map[, 2, 1] * prop_part
  vax_map[, 3, 1] <- vax_map[, 3, 1] * prop_full

  vax_map
}



##' @name run_onevax_xpvwrh
##' @title run model with a two-dose vaccine for input parameter sets, either
##' from initialisation or from equilibrium, those with waned vaccines are
##' eligible for revaccination (R), and return to the R stratum, those with
##' waned partial vaccines return to the unvaccinated stratum (X) and considered
##' immunologically naive. A user defined proportion of the population is
##' vaccine hesitant and is never vaccinated. Full acciantion (V) with 2 doses
##' gives maximum protection whereas partial vaccinaiton with 1 dose (P) gives
##' less.
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
##'  'gono_params' giving proportion of population accepting the first vaccine
##'  dose only
##' @param r1r2 scalar or numeric vector with same length as
##'  'gono_params' giving proportion of population accepting the both vaccine
##'  doses
##' @param booster_uptake scalar or numeric vector with same length as
##'  'gono_params' giving proportion of population undertaking booster
##'  vaccination after primary vaccination protection has waned.
##'   Defaults to supplied value of `r1r2`
##' @inheritParams run_onevax_xvwv
##' @return A list of transformed model outputs
##' @export
run_onevax_xpvwrh <- function(tt, gono_params, init_params = NULL,
                             dur_v = 1e3, dur_p = dur_v,
                             vea = 0, vei = 0, ved = 0, ves = 0,
                             dur_revax = dur_v,
                             vea_revax = vea, vei_revax = vei,
                             ved_revax = ved, ves_revax = ves,
                             vea_p = vea, vei_p = vei, ved_p = ved, ves_p = ves,
                             vbe = 0,
                             r1 = 0, r1r2 = 0,
                             booster_uptake = r1r2, strategy = NULL,
                             t_stop = 99, hes = 0) {

  stopifnot(all(lengths(list(booster_uptake, r1r2, r1, vea, vei,
                             ved, ves, vea_revax, vei_revax, ved_revax,
                             vea_p, vei_p, ves_p, ved_p,
                             dur_v, dur_p,
                             ves_revax, dur_revax)) %in%
                  c(1, length(gono_params))))

  vax_params <- Map(vax_params_xpvwrh, r1 = r1, r1r2 = r1r2,
                    booster_uptake = booster_uptake, dur_v = dur_v,
                    vea = vea, vei = vei, ved = ved, ves = ves,
                    vea_p = vea_p, vei_p = vei_p, ved_p = ved_p, ves_p = ves_p,
                    dur_revax = dur_revax, dur_p = dur_p,
                    vea_revax = vea_revax, vei_revax = vei_revax,
                    ved_revax = ved_revax, ves_revax = ves_revax, hes = hes,
                    MoreArgs = list(strategy = strategy,
                                    t_stop = t_stop, vbe = vbe))

  if (is.null(init_params)) {
    pars <- lapply(gono_params, model_params)
    init_params <- lapply(pars, initial_params_xpvwrh, hes = hes)
  }

  ret <- Map(run, gono_params = gono_params, vax_params = vax_params,
             init_params = init_params,
             MoreArgs = list(tt = tt))

  # name outputs
  ret <- lapply(ret, name_outputs, c("X", "P", "V", "W", "R", "H"))
  ret
}
