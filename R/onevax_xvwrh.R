##' Create initial conditions for the model
##' @name initial_params_xvwrh
##' @title Initial conditions for the model
##' @param pars A parameter list containing `N0`, `q`, `prev_Asl` and `prev_Ash`
##'   elements.
##' @param coverage  scalar giving initial coverage of vaccination in willing
##'  population
##' @param hes proportion of population vaccine hesitant
##' @return A list of initial conditions
##' @export
initial_params_xvwrh <- function(pars, coverage = 0, hes = 0) {
  assert_scalar_unit_interval(coverage)
  n_vax <- 5

  willing <- 1 - hes
  x_init <- willing * (1 - coverage)
  v_init <- willing * coverage
  cov <- c(x_init, v_init, 0, 0, hes)

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

##' @name vax_params_xvwrh
##' @title create vaccination parameters for use in onevax_xvwrh model
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
##' @param hes proportion of population vaccine hesitant
##' @return A list parameters in the model input format
vax_params_xvwrh <- function(vea = 0, vei = 0, ved = 0, ves = 0,
                             vea_revax = vea, vei_revax = vei,
                             ved_revax = ved, ves_revax = ves,
                             dur = 1e3, dur_revax = dur, primary_uptake = 0,
                             booster_uptake = 0, strategy = "VbE",
                             vbe = 0, t_stop = 99, hes = 0) {
  
  assert_character(strategy)
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
  # a proportion of all 'n' exist only in the hesitant compartment
  # there is no movement between the willing (x,v,w,r) and hesitant (h)
  # 1:x -> 2:v -> 3:w <-> 4:r
  # 5:h
  i_eligible <- c(1, 3)             #X and W are eligible for vaccination 
  i_w <- 3
  i_v <- c(2, 4)                    #V(2) and R(4) are protected 
  
  #number of compartments
  n_vax <- 5
  
  #make uptake vector
  uptake <- c(primary_uptake, booster_uptake)
  
  # ensure duration is not divided by 0
  ved <- min(ved, 1 - 1e-10)
  ved_revax <- min(ved_revax, 1 - 1e-10)
  
  p <- set_strategy_booster(strategy, uptake) 

  list(n_vax = n_vax,
       vbe   = create_vax_map_booster(n_vax, vbe, i_eligible, i_v),
       vod   = create_vax_map_booster(n_vax, p$vod, i_eligible, i_v),
       vos   = create_vax_map_booster(n_vax, p$vos, i_eligible, i_v),
       vea   = c(0, vea, 0, vea_revax, 0),
       vei   = c(0, vei, 0, vei_revax, 0),
       ved   = c(0, ved, 0, ved_revax, 0),
       ves   = c(0, ves, 0, ves_revax, 0),
       willing = c((1 - hes), 0, 0, 0, hes),
       w     = create_waning_map(n_vax, i_v, i_w, 1 / c(dur, dur_revax)),
       vax_t = c(0, t_stop),
       vax_y = c(1, 0)
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
##'  vaccination after primary vaccination protection has waned
##' @inheritParams run_onevax_xvwv
##' @return A list of transformed model outputs
##' @export
run_onevax_xvwrh <- function(tt, gono_params, init_params = NULL,
                             dur = 1e3, vea = 0, vei = 0, ved = 0, ves = 0,
                             dur_revax = dur,
                             vea_revax = vea, vei_revax = vei,
                             ved_revax = ved, ves_revax = ves,
                             vbe = 0, primary_uptake = 0,
                             booster_uptake = 0, strategy = "VbE",
                             t_stop = 99, hes = 0) {
  
  stopifnot(all(lengths(list(booster_uptake, primary_uptake, vea, vei, ved, ves, dur,
                             vea_revax, vei_revax, ved_revax, ves_revax,
                             dur_revax)) %in%
                  c(1, length(gono_params))))
  
  vax_params <- Map(vax_params_xvwrh, primary_uptake = primary_uptake,
                    booster_uptake = booster_uptake, dur = dur,
                    vea = vea, vei = vei, ved = ved, ves = ves,
                    dur_revax = dur_revax,
                    vea_revax = vea_revax, vei_revax = vei_revax,
                    ved_revax = ved_revax, ves_revax = ves_revax, hes = hes,
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

##' @name restart_hes
##' @title uses XVWRH model run in the absence of vaccination or hesitancy.
##' Saves down the number of individuals in each compartment, and moves
##' a given proportion (hes) of them from the X to the H strata to generate
##' new initial conditions in the presence of hesitancy.
##' @param y output of runonevax_xvwrh in absence of vaccination and hes = 0,
##' usually to equilibrium
##' @param n_vax number of stratum in model run y
##' @param hes proportion of population vaccine hesitant
##' @return A list of initial conditions to restart a model with n_vax
##' vaccination levels, and a populated hestitant stratum in the given
##' proportion 'hes'
##' @export

restart_hes <- function(y, n_vax = 5, hes = 0) {

  dim_y <- dim(y[["U"]])

   if (rowSums(y$N[, , n_vax])[dim_y[1]] > 0) {
    stop("Provided model run already contains hesitancy > 0")
   }

  if (rowSums(y$N[, , 2])[dim_y[1]] > 0) {
    stop("Provided model run has vaccination, baseline run should have all V
          = 0")
  }

  i_t <- dim_y[1]                     # number of timepoints
  n_vax <- n_vax %||% dim_y[3]        # number of strata

  n_vax_input <- dim_y[3]
  i_vax <- seq_len(min(n_vax,  n_vax_input))

  #create blank array, 2activity groups by number of strata
  U0 <- I0 <- A0 <- S0 <- T0 <- array(0, c(2, n_vax))

  # set compartments new initial conditions based on final position of y
  U0[, i_vax] <- y$U[i_t, , i_vax]
  I0[, i_vax] <- y$I[i_t, , i_vax]
  A0[, i_vax] <- y$A[i_t, , i_vax]
  S0[, i_vax] <- y$S[i_t, , i_vax]
  T0[, i_vax] <- y$T[i_t, , i_vax]

  # move correct number from equilibrium X to H
  U0[, n_vax] <- h <- U0[, 1] * hes
  U0[, 1] <- U0[, 1] - h

  I0[, n_vax] <- h <- I0[, 1] * hes
  I0[, 1] <- I0[, 1] - h

  A0[, n_vax] <- h <- A0[, 1] * hes
  A0[, 1] <- A0[, 1] - h

  S0[, n_vax] <- h <- S0[, 1] * hes
  S0[, 1] <- S0[, 1] - h

  T0[, n_vax] <- h <- T0[, 1] * hes
  T0[, 1] <- T0[, 1] - h

  list(U0 = U0, I0 = I0, A0 = A0, S0 = S0, T0 = T0, t = y$t[i_t])

}


##' @name create_vax_map_booster
##' @title Create mapping for movement between strata due to vaccination where
##' primary vaccination and booster vaccination can be different %s of uptake
##' to one another
##' @param n_vax Integer in (0, 5) denoting total number of strata
##' @param v vector of length 2 indicating % of population vaccinated according 
##' to the strategy, each value will apply to the eligible stratum separately
##' @param i_u indices of strata eligible for vaccination
##' @param i_v indices of strata being vaccinated
##' @return an array of the mapping

create_vax_map_booster <- function(n_vax, v, i_u, i_v) {           #note to self look at the vax map and see if its right!? 
  
  # ensure vaccine input is of correct length             vbe   = create_vax_map_booster(n_vax, vbe, i_eligible, i_v)
  n_group <- 2
  stopifnot(length(v) %in% c(1, n_group))
  stopifnot(all((v >= 0) & (v <= 1)))
  stopifnot(length(i_v) == length(i_u))
  stopifnot(max(i_u, i_v) <= n_vax)
  
  # set up vaccination matrix
  vax_map <- array(0, dim = c(n_group, n_vax, n_vax))
  
  if (length(v) == 1){
    v <- c(v,v)
  }
  
  for (i in seq_along(i_u)) {
    vax_map[, i_u[i], i_u[i]] <-  v[i]
    vax_map[, i_v[i], i_u[i]] <- -v[i]
  }
  print(vax_map)
  vax_map
}


set_strategy_booster <- function(strategy, uptake) {
  
  if (length(uptake) != 2) {
    stop("uptake must be length 2")
  }

  if (strategy == "VbE") {
    vos <- vod <- c(0,0)
  } else if (strategy == "VoD") {
    vod <- uptake
    vos <- c(0,0)
  } else if (strategy == "VoA") {
    vod <- uptake
    vos <- uptake
  } else if (strategy == "VoD(H)") {
    #vod <- c(0, uptake)
    vod <- uptake
    vos <- c(0,0)
  } else if (strategy == "VoA(H)") {
    #vod <- c(0, uptake)
    #vos <- c(0, uptake) temporary until i come up with a solution 
    vod <- uptake
    vos <- uptake
  } else if (strategy == "VoD(L)+VoA(H)") {
    vod <- uptake
    #vos <- c(0, uptake)
    vos <- uptake
  } else if (strategy == "VoS") {
    vod <- c(0,0)
    vos <- uptake
  } else {
    stop("strategy not recognised")
  }
  
  list(vod = vod, vos = vos)
}

check_gono_params <- function(pars) {
  with(pars, {
    assert_scalar_unit_interval(psi)
    assert_scalar_unit_interval(prev_Asl)
    assert_scalar_unit_interval(prev_Ash)
    assert_scalar_unit_interval(epsilon)
    assert_scalar_positive(sigma)
    assert_scalar_positive(nu)
    assert_scalar_positive(mu)
    assert_scalar_positive(rho)
  })
}