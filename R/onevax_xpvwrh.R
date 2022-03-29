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
##' @return A list parameters in the model input format
vax_params_xvwrh <- function(vea = 0, vei = 0, ved = 0, ves = 0,
                             vea_revax = vea, vei_revax = vei,
                             ved_revax = ved, ves_revax = ves,
                             dur_v = 1e3, dur_p = dur_v / 2, dur_revax = dur_v,
                             r1 = 0, r1r2 = 0,
                             booster_uptake = r1r2, strategy = NULL,
                             vbe = 0, t_stop = 99, hes = 0) {
  
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
  
  i_eligible <- c(1, 4)             #X and W are eligible for vaccination
  i_w <- 4
  i_v <- c(2, 3, 5)                    #P(2), V(3), and R(5) are protected
  
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
                                          r1 = r1, r1r2 = r1r2),
       vod     = create_vax_map_branching(n_vax, p$vod, i_eligible, i_v,
                                          r1 = r1, r1r2 = r1r2),
       vos     = create_vax_map_branching(n_vax, p$vos, i_eligible, i_v,
                                          r1 = r1, r1r2 = r1r2),
       vea     = c(0, vea, 0, vea_revax, 0),
       vei     = c(0, vei, 0, vei_revax, 0),                                #next step work out how to tweak efficacies for P and V
       ved     = c(0, ved, 0, ved_revax, 0),                                # need to be adujustable upstream 
       ves     = c(0, ves, 0, ves_revax, 0),
       w       = create_waning_map(n_vax, i_v, i_w, 1 / c(dur, dur_revax)),              #then waning! 
       vax_t   = c(0, t_stop),
       vax_y   = c(1, 0)
  )
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

create_vax_map_branching <- function(n_vax, v, i_u, i_v, r1, r1r2) {
  
  # ensure vaccine input is of correct length
  n_group <- 2
  stopifnot(length(v) == n_group)
  stopifnot(all(v %in% c(0, 1)))
  stopifnot(max(i_u, i_v) <= n_vax)
  
  # tweak eligibility , repeat over stratum 1 column 1 for ease
  i_u <- c(1, i_u)
  
  # set up vaccination matrix
  vax_map <- array(0, dim = c(n_group, n_vax, n_vax))
  
  for (i in seq_along(i_u)) {
    vax_map[, i_u[i], i_u[i]] <-  v
    vax_map[, i_v[i], i_u[i]] <- -v

  }                                        #maybe split this up otherwise  getting negatives for , , 4 as well! 
  
  #obtain proportions

  tot <- r1r2 + r1
  prop_full <- r1r2 / tot
  prop_part <- r1 / tot             #multiply these proportions through to , , 1
  
  vax_map[, 2, 1] <- vax_map[, 2, 1] * prop_part
  vax_map[, 3, 1] <- vax_map[, 3, 1] * prop_full
  
  vax_map
}