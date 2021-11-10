
## runs model_trial.R

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
