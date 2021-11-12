## sets up the starting conditions for the state variables in the presence of
## vaccination

initial_params_xvw_trial <- function(pars, coverage = 0) {
  assert_scalar_unit_interval(coverage)
  n_vax <- 3
  cov <- c(1 - coverage, coverage, 0)

 initial_params_trial(pars, n_vax, cov)
}


## sets up vaccination efficacies, who experiences the effects of vaccination,
## how waning occurs

vax_params_xvw_trial <- function(vea = 0, vei = 0, ved = 0, ves = 0,
                           dur = 1e3,
                           t_stop = 99) {

  assert_scalar_unit_interval(vea)
  assert_scalar_unit_interval(vei)
  assert_scalar_unit_interval(ved)
  assert_scalar_unit_interval(ves)
  assert_scalar_positive(dur)
  assert_scalar_positive(t_stop)

  # waned vaccinees move to own stratum, and are not eligible for re-vaccination
   i_v <- 2
   i_w <- n_vax <- 3

  # compartments to which vaccine efficacy applies
  ve <- c(0, 1, 0)
  ved <- min(ved, 1 - 1e-10) # ensure duration is not divided by 0

  list(n_vax = n_vax,
       vea   = vea * ve,
       vei   = vei * ve,
       ved   = ved * ve,
       ves   = ves * ve,
       w     = create_waning_map(n_vax, i_v, i_w, 1 / dur),
       vax_t = c(0, t_stop),
       vax_y = c(1, 0)
  )
}


## runs the trial model when supplied NHPs

run_onevax_xvw_trial <- function(tt, gono_params, initial_params_trial = NULL,
                           dur = 1e3,
                           vea = 0, vei = 0, ved = 0, ves = 0,
                           coverage = 0,
                           t_stop = 99) {

  stopifnot(all(lengths(list(vea, vei, ved, ves, dur)) %in%
                  c(1, length(gono_params))))
  assert_scalar_unit_interval(coverage)

  vax_params <- Map(vax_params_xvw_trial, dur = dur,
                    vea = vea, vei = vei, ved = ved, ves = ves,
                    MoreArgs = list(t_stop = t_stop))

  if (is.null(initial_params_trial)) {
    pars <- lapply(gono_params, model_params_trial)
    init_params_trial <- Map(initial_params_xvw_trial, pars = pars,
                             coverage = coverage)
  }

  ret <- Map(run_trial, gono_params = gono_params,
             init_params = init_params_trial,
             vax_params = vax_params,
             MoreArgs = list(tt = tt))

  ret
}
