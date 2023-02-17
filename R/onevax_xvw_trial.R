
##' Create initial conditions for the model in a vaccine trial
##' @name initial_params_xvw_trial
##' @title Initial conditions for the model in a vaccine trial
##' @param pars A parameter list containing `N0`, and `q` elements.
##' @param p_v scalar giving proportion of the trial cohort vaccinated
##' @param n_erlang integer giving the number of transitions that need to be made
##' through vaccine-protected strata until that protection has waned
##' @return A list of initial conditions.
##' @export

initial_params_xvw_trial <- function(pars, p_v = 0.5, n_erlang = 1) {

  assert_scalar_unit_interval(p_v)

  # XVW n_vax = 3, if n_erlang = 1, this is the same, if n_erlang > 1 this
  # needs to be accounted for with additional strata
  n_vax <- stratum_index_xvw_trial(n_erlang)$n_vax
  cov <- c(1 - p_v, p_v, rep(0, n_erlang - 1), 0)
  initial_params_trial(pars, n_vax, cov)
}

## sets up vaccination efficacies, who experiences the effects of vaccination,
## how waning occurs

##' @name vax_params_xvw_trial
##' @title Create vaccination parameters for use in onevax_xvw_trial model,
##' assign who experiences vaccine effects, and how waning occurs.
##' @param vea scalar indicating efficacy of the vaccine against acquisition
##' (between 0-1)
##' @param vei scalar indicating efficacy of the vaccine against infectiousness
##' (between 0-1)
##' @param ved scalar indicating efficacy of the vaccine against duration
##' (between 0-1)
##' @param ves scalar indicating efficacy of the vaccine against symptoms
##' (between 0-1)
##' @param dur scalar indicating duration of the vaccine (in years)
##' @param n_erlang integer giving the number of transitions that need to be made
##' through vaccine-protected strata until that protection has waned
##' @return A list of parameters in the model input format

vax_params_xvw_trial <- function(vea = 0, vei = 0, ved = 0, ves = 0,
                           dur = 1e3, n_erlang = 1) {

  assert_scalar_unit_interval(vea)
  assert_scalar_unit_interval(vei)
  assert_scalar_unit_interval(ved)
  assert_scalar_unit_interval(ves)
  assert_scalar_positive(dur)

  # generate indicies for all strata and
  idx <- stratum_index_xvw_trial(n_erlang)
  
  # waned vaccinees move through erlang compartments until they reach
  # the final waned compartment with no protection

  #X + W + n_erlang = total number of strata
   n_vax <- idx$n_vax

  # waned vaccinees move to own stratum, and are not eligible for re-vaccination
  # generate i_v and i_w

  # strata that individuals wane from, i.e all 'V' strata
   i_v <- idx$V

  # strata that individuals wane to
   i_w <- idx$V + 1
   
  # compartments to which vaccine efficacy applies
  ve <- c(0, rep(1, n_erlang), 0)
  ved <- min(ved, 1 - 1e-10) # ensure duration is not divided by 0

  list(n_vax = n_vax,
       vea   = vea * ve,
       vei   = vei * ve,
       ved   = ved * ve,
       ves   = ves * ve,
       w     = create_waning_map_trial(n_vax, i_v, i_w, (n_erlang / dur))
  )
}


## runs the trial model when supplied NHPs

##' @name run_onevax_xvw_trial
##' @title Run model with single vaccine for input parameter sets, either from
##' initialisation or from equilibrium, those with waned vaccines are not
##' eligible for re-vaccination.
##' @param gono_params list of gono params for a vaccination trial
##' @param initial_params_trial list of initial conditions for model trial. Set
##' default as NULL.
##' @param vea scalar or numeric vector with same length as `gono_params` giving
##'  efficacy of the vaccine against acquisition (between 0-1)
##' @param vei scalar or numeric vector with same length as `gono_params` giving
##'  efficacy of the vaccine against infectiousness (between 0-1)
##' @param ved scalar or numeric vector with same length as `gono_params` giving
##'  efficacy of the vaccine against duration (between 0-1)
##' @param ves scalar or numeric vector with same length as `gono_params` giving
##'  efficacy of the vaccine against symptoms (between 0-1)
##' @param dur  scalar or numeric vector with same length as `gono_params`
##'  giving duration of the vaccine (in years)
##' @param p_v scalar giving proportion of the trial cohort vaccinated, default
##'  is 0.5.
##' @param n_erlang integer giving the number of transitions that need to be made
##' through vaccine-protected strata until that protection has waned
##' @inheritParams run_trial
##' @inheritParams vax_params_xvw_trial
##' @export


run_onevax_xvw_trial <- function(tt, gono_params, initial_params_trial = NULL,
                           dur = 1e3,
                           vea = 0, vei = 0, ved = 0, ves = 0,
                           p_v = 0.5, n_erlang = 1) {

  stopifnot(all(lengths(list(vea, vei, ved, ves, dur)) %in%
                  c(1, length(gono_params))))
  assert_scalar_unit_interval(p_v)

  vax_params <- Map(vax_params_xvw_trial, dur = dur,
                    vea = vea, vei = vei, ved = ved, ves = ves,
                    n_erlang = n_erlang)

  if (is.null(initial_params_trial)) {
    pars <- lapply(gono_params, model_params_trial)
    init_params_trial <- Map(initial_params_xvw_trial, pars = pars,
                             p_v = p_v, n_erlang = n_erlang)
  }

  ret <- Map(run_trial, gono_params = gono_params,
             init_params = init_params_trial,
             vax_params = vax_params,
             MoreArgs = list(tt = tt))

  ret
}


##' @name stratum_index_xvw_trial
##' @title Generate the indices of all xvw trial strata
##' @param n_erlang integer giving the number of transitions that need to be 
##' made through vaccine-protected strata until that protection has waned
##' @return A list of strata with their indicies
##' @export

stratum_index_xvw_trial <- function(n_erlang) {
  
  ret <- list(X = 1)
    
    ret$V <- max(ret$X) + seq_len(n_erlang)
    ret$W <- max(ret$V) + 1
    ret$n_vax <- ret$W
  
  ret
  
}

