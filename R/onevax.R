
vax_params1 <- function(dur = 4, ve = 1, vs = 0, vd = 0, eff = 0.31) {
  z <- 1 / dur
  list(n_vax = 2,
       ve    = c(1 - ve, ve),
       vs    = vs,
       vd    = vd,
       eff   = c(0, eff),
       v     = rbind(c(1, 0),
                     c(-1, 0)),
       w     = rbind(c(0,  z),
                     c(0, -z))
  )
}

run_onevax_int <- function(n = NULL, tt, eff, dur, ve, vd, vs, equilib) {
  init_params <- NULL
  if (equilib) {
    if (length(n) > 1) stop("if length(n) > 1, equilib must be FALSE")
    init_params <- restart_params(novax_equilib(n))
  }
  vax_params <- vax_params1(dur = dur, ve = ve, vs = vs, vd = vd, eff = eff)
  pars <- model_params(gono_params = gono_params(n),
                       init_params = init_params,
                       vax_params = vax_params)
  mod <- model(user = pars, unused_user_action = FALSE)
  y <- mod$run(tt)
  mod$transform_variables(y)
}

##' @name run_onevax
##' @title run model without vaccination for input parameter sets, either from
##' initialisation or from equilibrium
##' @param n an integer vector (or value) containing the indices of
##' corresponding parameter set (1:982). If `n = NULL` the equilibrium positions
##' for the full parameter set are returned
##' @param tt a numeric vector of times at which the model state is output
##' @param eff single numeric indicating efficacy of the vaccine (between 0-1)
##' @param dur single numeric indicating duration of the vaccine (in years)
##' @param ve single numeric indicating % of population vaccinated before entry
##'  (between 0-1)
##' @param vd single numeric indicating % of population vaccinated on diagnosis
##' (between 0-1)
##' @param vs single numeric indicating % of population vaccinated on screening
##' (between 0-1)
##' @param equilib a logical indicating whether to run from equilibrium, default
##' is `FALSE`, i.e. run from initial conditions
##' @return A list of transformed model outputs
##' @export
run_onevax <- function(n = NULL, tt, eff, dur,
                       ve = 0, vd = 0, vs = 0,
                       equilib = FALSE) {
  n <- n %||% seq_len(nrow(gono_params()))
  lapply(n, run_onevax_int,
         tt = tt, eff = eff, dur = dur,
         ve = ve, vd = vd, vs = vs,
         equilib = equilib)
}
