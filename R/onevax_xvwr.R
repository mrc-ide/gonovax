##' @name vax_params_xvwr
##' @title create vaccination parameters for use in onevax_xvwr model
##' @inheritParams vax_params_xvwv
##' @return A list parameters in the model input format
vax_params_xvwr <- function(eff = 0, dur = 1e3,
                        ve = 0, vs = 0, vd = 0, t_stop = 99) {
  # waned vaccinees move to own stratum, but are eligible for re-vaccination
  # re-vaccination is into a fourth stratum (r)
  # 1:x -> 2:v -> 3:w <-> 4:r
  i_eligible <- c(1, 3)
  i_w <- 3
  i_v <- c(2, 4)
  n_vax <- 4

  list(n_vax = n_vax,
       ve    = create_vax_map(n_vax, ve, i_eligible, i_v),
       vd    = create_vax_map(n_vax, vd, i_eligible, i_v),
       vs    = create_vax_map(n_vax, vs, i_eligible, i_v),
       eff   = c(0, eff, 0),
       w     = create_waning_map(n_vax, i_v, i_w, 1 / dur),
       vax_t = c(0, t_stop),
       vax_y = c(1, 0)
  )
}

##' @name run_onevax_xvwr
##' @title run model with single vaccine for input parameter sets, either from
##' initialisation or from equilibrium, those with waned vaccines are eligible
##' for revaccination (R), and return to the R stratum
##' @inheritParams run_onevax_xvwv
##' @return A list of transformed model outputs
##' @export
run_onevax_xvwr <- function(tt, gono_params, init_params = NULL,
                          eff, dur, ve = 0, vd = 0, vs = 0,
                          t_stop = 99) {

  vax_params <- Map(vax_params_xvwr, eff = eff, dur = dur,
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
