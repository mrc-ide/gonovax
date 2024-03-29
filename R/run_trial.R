##' @name run_trial
##' @title Run odin model of gonorrhoea vaccine trial with or without
##' vaccination
##' @param tt a numeric vector of times at which the model state is output in
##' years
##' @param gono_params a data frame of parameters
##' @param init_params = NULL
##' @param vax_params = NULL
##' @param transform = TRUE
##' @param stochastic logical indicating if the run should be made with the
##' default deterministic trial model in continuous time or stochastic trial
##' model in discrete time
##' @param N integer to assign the total number of individuals in the trial
##' (split equally across the two arms)
##' @export run_trial


run_trial <- function(tt, gono_params, init_params = NULL, vax_params = NULL,
                      transform = TRUE, stochastic = FALSE, N = 6e5) {

  demographic_params_trial <-  demographic_params_trial(N = N)
  ret <- c(demographic_params_trial, gono_params)
  pars <- c(ret, init_params, vax_params)

  if (stochastic == TRUE) {
    mod <- model_trial_stochastic$new(user = pars, unused_user_action = FALSE)
    tt <- round(tt * 365)

  } else {
    mod <- model_trial$new(user = pars, unused_user_action = FALSE)
  }

  y <- mod$run(tt)

  if (transform) {
    y <- mod$transform_variables(y)
  }

  y
}
