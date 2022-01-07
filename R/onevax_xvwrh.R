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

  initial_params(pars, n_vax, cov)
}

##' @name vax_params_xvwrh
##' @title create vaccination parameters for use in vax_params_xvwv model
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
##' @return A list parameters in the model input format
vax_params_xvwrh <- function(vea = 0, vei = 0, ved = 0, ves = 0,
                            vea_revax = vea, vei_revax = vei,
                            ved_revax = ved, ves_revax = ves,
                            dur = 1e3, dur_revax = dur, uptake = 0,
                            strategy = "VbE",
                            vbe = 0, t_stop = 99) {

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
  assert_scalar_unit_interval(uptake)
  assert_scalar_unit_interval(vbe)
  assert_scalar_positive(t_stop)
  # waned vaccinees move to own stratum, but are eligible for re-vaccination
  # re-vaccination is into a fourth stratum (r)
  # a proportion of all 'n' exist only in the hesitant compartment
  # there is no movement between the willing (x,v,w,r) and hesitant (h)
  # 1:x -> 2:v -> 3:w <-> 4:r
  # 5:h
  i_eligible <- c(1, 3)
  i_w <- 3
  i_v <- c(2, 4)

  #number of compartments
  n_vax <- 5

  # ensure duration is not divided by 0
  ved <- min(ved, 1 - 1e-10)
  ved_revax <- min(ved_revax, 1 - 1e-10)

  p <- set_strategy(strategy, uptake)

  list(n_vax = n_vax,
       vbe   = create_vax_map(n_vax, vbe, i_eligible, i_v),
       vod   = create_vax_map(n_vax, p$vod, i_eligible, i_v),
       vos   = create_vax_map(n_vax, p$vos, i_eligible, i_v),
       vea   = c(0, vea, 0, vea_revax, 0),
       vei   = c(0, vei, 0, vei_revax, 0),
       ved   = c(0, ved, 0, ved_revax, 0),
       ves   = c(0, ves, 0, ves_revax, 0),
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
##' @inheritParams run_onevax_xvwv
##' @return A list of transformed model outputs
##' @export
run_onevax_xvwrh <- function(tt, gono_params, init_params = NULL,
                            dur = 1e3, vea = 0, vei = 0, ved = 0, ves = 0,
                            dur_revax = dur,
                            vea_revax = vea, vei_revax = vei,
                            ved_revax = ved, ves_revax = ves,
                            vbe = 0, uptake = 0, strategy = "VbE",
                            t_stop = 99, hes = 0) {

  stopifnot(all(lengths(list(uptake, vea, vei, ved, ves, dur,
                             vea_revax, vei_revax, ved_revax, ves_revax,
                             dur_revax)) %in%
                  c(1, length(gono_params))))

  vax_params <- Map(vax_params_xvwrh, uptake = uptake, dur = dur,
                    vea = vea, vei = vei, ved = ved, ves = ves,
                    dur_revax = dur_revax,
                    vea_revax = vea_revax, vei_revax = vei_revax,
                    ved_revax = ved_revax, ves_revax = ves_revax,
                    MoreArgs = list(strategy = strategy,
                                    t_stop = t_stop, vbe = vbe))

  if (is.null(init_params)) {

    pars <- lapply(gono_params, model_params)
    init_params <- lapply(pars, initial_params_xvwrh, hes = hes)

    ret <- Map(run, gono_params = gono_params, vax_params = vax_params,
               init_params = init_params,
               MoreArgs = list(tt = tt))

  } else {
    ret <- Map(run, gono_params = gono_params, init_params = init_params,
               vax_params = vax_params,
               MoreArgs = list(tt = tt))
  }

  # name outputs
  ret <- lapply(ret, name_outputs, c("X", "V", "W", "R", "H"))
  ret
}
