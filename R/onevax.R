
vax_params1 <- function(dur = 1e3, ve = 0, vs = 0, vd = 0, eff = 0) {
  z <- 1 / dur
  v <- rbind(c(1, 0),
            c(-1, 0))
  list(n_vax = 2,
       ve    = c(1 - ve, ve),
       vs    = vs * v,
       vd    = vd * v,
       eff   = c(0, eff),
       w     = rbind(c(0,  z),
                     c(0, -z))
  )
}

run_onevax_int <- function(n, tt, eff, dur, ve, vd, vs, equilib) {

  if (length(n) != 1) stop("length(n) must equal 1")

  vax_params <- vax_params1(dur = dur, ve = ve, vs = vs, vd = vd, eff = eff)

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

##' @name run_onevax
##' @title run model with single vaccine for input parameter sets, either from
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
