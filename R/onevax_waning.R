
##' @name vax_params1w
##' @title create vaccination parameters for use in onevax_waning model
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
##' @return A list parameters in the model input format
vax_params1w <- function(eff = 0, dur = 1e3,
                         ve = 0, vs = 0, vd = 0) {
  # indices of vaccine strata (u: unvaccinated, v: vaccinated, w: waned)
  i_u <- 1
  i_v <- 2
  i_w <- n_vax <- 5

  #' currently waning is in stratum 5 (with strata 3:4 unused) for future
  #' compatibility with dualvax, however may get performance boost by moving to
  #' stratum 3

  list(n_vax = n_vax,
       ve    = c(1 - ve, ve, 0, 0, 0),
       vs    = create_vax_map(n_vax, vs),
       vd    = create_vax_map(n_vax, vd),
       eff   = c(0, eff, 0, 0, 0),
       w     = create_waning_map(n_vax, i_w, 1 / dur)
  )
}

run_onevaxw_int <- function(n, tt, eff, dur, ve, vd, vs, equilib) {

  if (length(n) != 1) stop("length(n) must equal 1")

  vax_params <- vax_params1w(eff = eff, dur = dur,
                             ve = ve, vs = vs, vd = vd)

  if (equilib) {
    init_params <- restart_params(novax_equilib(n), n_vax = vax_params$n_vax)
  } else {
    # set inital params based on gono_params() input
    init_params <- NULL
  }

  pars <- model_params(gono_params = gono_params(n),
                       init_params = init_params,
                       vax_params = vax_params)
  mod <- model(user = pars, unused_user_action = FALSE)
  y <- mod$run(tt)
  mod$transform_variables(y)
}

##' @name run_onevax_waning
##' @title run model with single vaccine for input parameter sets, either from
##' initialisation or from equilibrium, those with waned vaccines are ineligible
##' for revaccination
##' @param n an integer vector (or value) containing the indices of
##' corresponding parameter set (1:982). If `n = NULL` the equilibrium positions
##' for the full parameter set are returned
##' @param tt a numeric vector of times at which the model state is output
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
##' @param equilib a logical indicating whether to run from equilibrium, default
##' is `FALSE`, i.e. run from initial conditions
##' @return A list of transformed model outputs
##' @export
run_onevax_waning <- function(n = NULL, tt, eff, dur,
                       ve = 0, vd = 0, vs = 0,
                       equilib = FALSE) {
  n <- n %||% seq_len(nrow(gono_params()))
  lapply(n, run_onevaxw_int,
         tt = tt, eff = eff, dur = dur,
         ve = ve, vd = vd, vs = vs,
         equilib = equilib)
}
