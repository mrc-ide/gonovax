##' @name vax_params_xvwv
##' @title create vaccination parameters for use in onevax_xvwv model
##' @param vea scalar indicating efficacy of the vaccine against acquisition
##' (between 0-1)
##' @param vei scalar indicating efficacy of the vaccine against infectiousness
##' (between 0-1)
##' @param ved scalar indicating efficacy of the vaccine against duration
##' (between 0-1)
##' @param ves scalar indicating efficacy of the vaccine against symptoms
##' (between 0-1)
##' @param dur scalar indicating duration of the vaccine (in years)
##' @param vbe scalar indicating pc of population vaccinated before entry
##'  (between 0-1)
##' @param uptake scalar indicating pc of population vaccinated as part
##'  of strategy
##' @param strategy single character string in "VbE", "VoD", "VoD(H)",
##'  "VoA", "VoA(H)", "VoD(L)+VoA(H)"
##' @param t_stop time at which vaccination should stop (years)
##' @return A list parameters in the model input format
vax_params_xvwv <- function(vea = 0, vei = 0, ved = 0, ves = 0,
                            dur = 1e3, uptake = 0, strategy = "VbE",
                            vbe = 0, t_stop = 99) {

  assert_character(strategy)
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
  i_eligible <- c(1, 3)
  i_v <- c(2, 2)
  i_w <- n_vax <- 3

  # compartments to which vaccine efficacy applies
  ve <- c(0, 1, 0)
  ved <- min(ved, 1 - 1e-10) # ensure duration is not divided by 0

  p <- set_strategy(strategy, uptake)

  list(n_vax = n_vax,
       vbe   = create_vax_map(n_vax, vbe, i_eligible, i_v),
       vod   = create_vax_map(n_vax, p$vod, i_eligible, i_v),
       vos   = create_vax_map(n_vax, p$vos, i_eligible, i_v),
       vea   = vea * ve,
       vei   = vei * ve,
       ved   = ved * ve,
       ves   = ves * ve,
       willing = c(1, 0, 0),
       w     = create_waning_map(n_vax, i_v, i_w, 1 / dur),
       vax_t = c(0, t_stop),
       vax_y = c(1, 0)
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
##' @inheritParams run
##' @inheritParams vax_params_xvwv
##' @export
run_onevax_xvwv <- function(tt, gono_params, init_params = NULL, dur = 1e3,
                            vea = 0, vei = 0, ved = 0, ves = 0, vbe = 0,
                            uptake = 0, strategy = "VbE",
                            t_stop = 99) {


  stopifnot(all(lengths(list(uptake, vea, vei, ved, ves, dur)) %in%
                  c(1, length(gono_params))))

  vax_params <- Map(vax_params_xvwv, uptake = uptake, dur = dur,
                    vea = vea, vei = vei, ved = ved, ves = ves,
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
  ret <- lapply(ret, name_outputs, c("X", "V", "W"))

  ret
}
