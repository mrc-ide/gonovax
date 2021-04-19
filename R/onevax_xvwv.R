##' @name vax_params_xvwv
##' @title create vaccination parameters for use in onevax_xvwv model
##' @param eff scalar indicating efficacy of the vaccine (between 0-1)
##' @param dur scalar indicating duration of the vaccine (in years)
##' @param ve scalar indicating pc of population vaccinated before entry
##'  (between 0-1)
##' @param uptake scalar indicating pc of population vaccinated as part
##'  of strategy
##' @param strategy single character string in "VbE", "VoD(all)", "VoD(H)",
##'  "VoA(all)", "VoA(H)", "VoD(L)+VoA(H)"
##' @param t_stop time at which vaccination should stop (years)
##' @return A list parameters in the model input format
vax_params_xvwv <- function(eff = 0, dur = 1e3, uptake = 0, strategy = "VbE",
                        ve = 0, t_stop = 99) {

  assert_character(strategy)
  assert_scalar(eff)
  assert_scalar(dur)
  assert_scalar(uptake)
  assert_scalar(ve)
  assert_scalar(t_stop)
  # waned vaccinees move to own stratum, but are eligible for re-vaccination
  # 1:x -> 2:v <-> 3:w
  i_eligible <- c(1, 3)
  i_v <- c(2, 2)
  i_w <- n_vax <- 3

  p <- set_strategy(strategy, uptake)

  list(n_vax = n_vax,
       ve    = create_vax_map(n_vax, ve, i_eligible, i_v),
       vd    = create_vax_map(n_vax, p$vd, i_eligible, i_v),
       vs    = create_vax_map(n_vax, p$vs, i_eligible, i_v),
       eff   = c(0, eff, 0),
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
##' @param eff scalar or numeric vector with same length as `gono_params` giving
##'  efficacy of the vaccine (between 0-1)
##' @param dur  scalar or numeric vector with same length as `gono_params`
##'  giving duration of the vaccine (in years)
##' @param uptake  scalar or numeric vector with same length as `gono_params`
##'  giving pc of population vaccinated as part of strategy
##' @inheritParams run
##' @inheritParams vax_params_xvwv
##' @export
run_onevax_xvwv <- function(tt, gono_params, init_params = NULL,
                          eff, dur, ve = 0, uptake = 0, strategy = "VbE",
                          t_stop = 99) {


  stopifnot(all(lengths(list(uptake, eff, dur)) %in% c(1, length(gono_params))))

  vax_params <- Map(vax_params_xvwv, uptake = uptake, eff = eff, dur = dur,
                    MoreArgs = list(strategy = strategy, t_stop = t_stop,
                                    ve = ve))

  if (is.null(init_params)) {
    ret <- Map(run, gono_params = gono_params, vax_params = vax_params,
              MoreArgs = list(tt = tt))
  } else {
    ret <- Map(run, gono_params = gono_params, init_params = init_params,
               vax_params = vax_params,
               MoreArgs = list(tt = tt))
  }

  ret
}
