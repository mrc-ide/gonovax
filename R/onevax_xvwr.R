##' @name vax_params_xvwr
##' @title create vaccination parameters for use in onevax_xvwr model
##' @inheritParams vax_params_xvwv
##' @param eff_revax efficacy of revaccination, default to same as primary
##' @param dur_revax duration of protection for revaccination,
##'  default to same as primary
##' @return A list parameters in the model input format
vax_params_xvwr <- function(eff = 0, eff_revax = eff,
                            dur = 1e3, dur_revax = dur, uptake = 0,
                            strategy = "VbE",
                            ve = 0, t_stop = 99) {

  assert_character(strategy)
  assert_scalar(eff)
  assert_scalar(eff_revax)
  assert_scalar(dur)
  assert_scalar(dur_revax)
  assert_scalar(uptake)
  assert_scalar(ve)
  assert_scalar(t_stop)
  # waned vaccinees move to own stratum, but are eligible for re-vaccination
  # re-vaccination is into a fourth stratum (r)
  # 1:x -> 2:v -> 3:w <-> 4:r
  i_eligible <- c(1, 3)
  i_w <- 3
  i_v <- c(2, 4)
  n_vax <- 4

  p <- set_strategy(strategy, uptake)

  list(n_vax = n_vax,
       ve    = create_vax_map(n_vax, ve, i_eligible, i_v),
       vd    = create_vax_map(n_vax, p$vd, i_eligible, i_v),
       vs    = create_vax_map(n_vax, p$vs, i_eligible, i_v),
       eff   = c(0, eff, 0, eff_revax),
       w     = create_waning_map(n_vax, i_v, i_w, 1 / c(dur, dur_revax)),
       vax_t = c(0, t_stop),
       vax_y = c(1, 0)
  )
}

##' @name run_onevax_xvwr
##' @title run model with single vaccine for input parameter sets, either from
##' initialisation or from equilibrium, those with waned vaccines are eligible
##' for revaccination (R), and return to the R stratum
##' @param eff_revax scalar or numeric vector with same length as `gono_params`
##'  giving efficacy of revaccination, default to same as primary
##' @param dur_revax scalar or numeric vector with same length as `gono_params`
##'  giving duration of protection for revaccination, default to same as primary
##' @inheritParams run_onevax_xvwv
##' @return A list of transformed model outputs
##' @export
run_onevax_xvwr <- function(tt, gono_params, init_params = NULL,
                            eff, dur,
                            eff_revax = eff, dur_revax = dur,
                            ve = 0, uptake = 0, strategy = "VbE",
                            t_stop = 99) {

  stopifnot(all(lengths(list(uptake, eff, dur, eff_revax, dur_revax)) %in%
                  c(1, length(gono_params))))

  vax_params <- Map(vax_params_xvwr, uptake = uptake, eff = eff, dur = dur,
                    eff_revax = eff_revax, dur_revax = dur_revax,
                    MoreArgs = list(strategy = strategy,
                                    t_stop = t_stop, ve = ve))

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
