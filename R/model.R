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
##' @param n_diag_rec integer for the number of diagnosis history substrata
##' @param transform = TRUE
##' @export run

run <- function(tt, gono_params, init_params = NULL, vax_params = NULL, n_diag_rec = 1,
                transform = TRUE) {

  pars <- model_params(gono_params = gono_params,
                       init_params = init_params,
                       vax_params = vax_params,
                       n_diag_rec = n_diag_rec)
  mod <- model$new(user = pars, unused_user_action = FALSE)
  y <- mod$run(tt)

  if (transform) {
    y <- mod$transform_variables(y)
  }

  y
}

##' @name run_Diagnosis
##' @title Run odin model of gonorrhoea with or without vaccination
##' @param tt a numeric vector of times at which the model state is output
##' @param gono_params a data frame of parameters
##' @param init_params = NULL
##' @param vax_params = NULL
##' @param n_erlang integer giving the number of transitions that need to be
##'  made through vaccine-protected strata until that protection has waned
##'  @param n_diag_rec integer for the number of diagnosis history substrata
##' @param transform = TRUE
##' @export run

run_xpvwrh <- function(tt, gono_params, init_params = NULL, vax_params = NULL, n_erlang = 1, n_diag_rec = 1,
                                 transform = TRUE) {
  

  
  pars <- model_params_xpvwrh(gono_params = gono_params,
                                        init_params = init_params,
                                        vax_params = vax_params,
                                        n_erlang = n_erlang,
                                        n_diag_rec = n_diag_rec)
  

  
  mod <- model$new(user = pars, unused_user_action = FALSE)
  y <- mod$run(tt)

  if (transform) {
    y <- mod$transform_variables(y)
  }
  
  y
}