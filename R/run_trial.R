##' @name run_trial
##' @title Run odin model of gonorrhoea vaccine trial with or without
##' vaccination
##' @param tt a numeric vector of times at which the model state is output
##' @param gono_params a data frame of parameters
##' @param init_params = NULL
##' @param vax_params = NULL
##' @param transform = TRUE
##' @export run


run_trial <- function(tt, gono_params, init_params = NULL, vax_params = NULL,
                transform = TRUE) {

  demographic_params_trial <-  demographic_params_trial()
  ret <- c(demographic_params_trial, gono_params)
  pars <- c(ret, init_params, vax_params)

  mod <- model_trial$new(user = pars, unused_user_action = FALSE)
  y <- mod$run(tt)

  if (transform) {
    y <- mod$transform_variables(y)
  }

  y
}
