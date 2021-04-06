##' @name vax_params_xvwr
##' @title create vaccination parameters for use in onevax_xvwr model
##' @inheritParams vax_params_xvwv
##' @return A list parameters in the model input format
vax_params_xvwr <- function(eff = 0, dur = 1e3,
                        ve = 0, vs = 0, vd = 0, t_stop = 99) {
  # waned vaccinees move to own compartment, but are eligible for re-vaccination
  i_u <- c(1, 3)
  i_w <- n_vax <- 3

  list(n_vax = n_vax,
       ve    = create_vax_map(n_vax, ve, i_u),
       vd    = create_vax_map(n_vax, vd, i_u),
       vs    = create_vax_map(n_vax, vs, i_u),
       eff   = c(0, eff, 0),
       w     = create_waning_map(n_vax, i_w, 1 / dur),
       vax_t = c(0, t_stop),
       vax_y = c(1, 0)
  )
}

##' @name run_onevax_xvwv
##' @title run model with single vaccine for input parameter sets, either from
##' initialisation or from equilibrium, those with waned vaccines are eligible
##' for revaccination, and return to the V compartment
##' @inheritParams run
##' @param eff single numeric indicating efficacy of the vaccine (between 0-1)
##' @param dur single numeric indicating duration of the vaccine (in years)
##' @param ve single numeric indicating % of population vaccinated before entry
##'  (between 0-1)
##' @param vd numeric indicating % of population vaccinated on diagnosis
##' (between 0-1) a scalar will apply to both activity groups, a vector of
##' length 2 will apply to each activity group separately
##' @param vs numeric indicating % of population vaccinated on screening
##' (between 0-1) a scalar will apply to both activity groups, a vector of
##' length 2 will apply to each activity group separately
##' @param t_stop time at which vaccination should stop (years)
##' @return A list of transformed model outputs
##' @export
run_onevax_xvwv <- function(tt, gono_params, init_params = NULL,
                          eff, dur, ve = 0, vd = 0, vs = 0,
                          t_stop = 99) {

  vax_params <- Map(vax_params_xvwv, eff = eff, dur = dur,
                    ve = ve, vs = vs, vd = vd, t_stop = t_stop)

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
