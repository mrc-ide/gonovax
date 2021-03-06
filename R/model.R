##' @name model
##' @title Model of gonorrhoea with dual vaccines
##'This is an odin model.
##' @export model
NULL



##' @name run
##' @title Run odin model of gonorrhoea with or without vaccination
##' @param tt a numeric vector of times at which the model state is output
##' @param gono_params a data frame of parameters
##' @param init_params = NULL
##' @param vax_params = NULL
##' @param transform = TRUE
##' @export run

run <- function(tt, gono_params, init_params = NULL, vax_params = NULL,
                transform = TRUE) {

  pars <- model_params(gono_params = gono_params,
                                init_params = init_params,
                                vax_params = vax_params)
  mod <- model(user = pars, unused_user_action = FALSE)
  y <- mod$run(tt)

  if (transform) {
    y <- mod$transform_variables(y)
  }

  y
}
