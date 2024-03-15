
##' Create initial conditions for the model in a vaccine trial
##' @name initial_params_xvw_trial
##' @title Initial conditions for the model in a vaccine trial
##' @param pars A parameter list containing `N0`, and `q` elements.
##' @param p_v scalar giving proportion of the trial cohort vaccinated
##' @param n_erlang integer giving the number of transitions that need to be
##'  made through vaccine-protected strata until that protection has waned
##' @param n_diag_rec integer giving the number of each X, V(erlang), and W
##' stratum, allowing tracking of diagnosis history. e.g for a n_diag_rec = 2
##' and erlang = 1, there will be X.I, X.II, V1.I, V1.II, W.I, W.II strata.
##' Where '.I' corresponds to never-diagnosed individuals and '.II' is for
##' individuals diagnosed at least once.
##' @return A list of initial conditions.
##' @export

initial_params_xvw_trial <- function(pars, p_v = 0.5, n_erlang = 1,
                                     n_diag_rec = 1) {
  assert_scalar_unit_interval(p_v)

  # XVW n_vax = 3, if n_erlang = 1, this is the same, if n_erlang > 1 this
  # needs to be accounted for with additional strata
  n_vax <- stratum_index_xvw_trial(n_erlang, n_diag_rec)$n_vax
  cov <- c(1 - p_v, rep(0, n_diag_rec - 1), p_v,
           rep(0, (n_diag_rec * n_erlang) - 1),
           rep(0, n_diag_rec))
  initial_params_trial(pars, n_vax, cov, n_diag_rec)
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
##' @param n_erlang integer giving the number of transitions that need to be
##'  made
##' @param stochastic logical indicating if the parameters are for the
##' default deterministic trial model in continuous time or stochastic trial
##' model in discrete time
##' through vaccine-protected strata until that protection has waned
##' @param n_diag_rec integer giving the number of each X, V(erlang), and W
##' stratum, allowing tracking of diagnosis history. e.g for a n_diag_rec = 2
##' and erlang = 1, there will be X.I, X.II, V1.I, V1.II, W.I, W.II strata.
##' Where '.I' corresponds to never-diagnosed individuals and '.II' is for
##' individuals diagnosed at least once.
##' @return A list of parameters in the model input format

vax_params_xvw_trial <- function(vea = 0, vei = 0, ved = 0, ves = 0,
                                 dur = 1e3, n_erlang = 1, stochastic = FALSE,
                                 n_diag_rec = 1) {

  assert_scalar_unit_interval(vea)
  assert_scalar_unit_interval(vei)
  assert_scalar_unit_interval(ved)
  assert_scalar_unit_interval(ves)
  assert_scalar_positive(dur)

  # generate indices for all strata and
  idx <- stratum_index_xvw_trial(n_erlang, n_diag_rec)

  # waned vaccinees move through erlang compartments until they reach
  # the final waned compartment with no protection

  #X + W + n_erlang = total number of strata
  n_vax <- idx$n_vax

  # waned vaccinees move to own stratum, and are not eligible for re-vaccination
  # generate i_v and i_w

  # strata that individuals wane from, i.e all 'V' strata
  i_v <- idx$V

  # strata that individuals wane to
  i_w <- idx$V + n_diag_rec

  # diagnosed individuals move to the next diagnosis-history stratum (if
  # n_diag_rec > 1). These history strata are the same in their characteristics
  # as the stratum the individual moved from e.g Va and Vb both experience
  # protection and waning rate the same. This is for downstream calculation
  # of person-years exposed etc.

  i <- seq_len(idx$n_vax)
  # diagnosed from
  i_eligible <- i[i %% n_diag_rec != 0]

  # diagnosed to
  i_p <- i[i %% n_diag_rec != 1]

  # create diagnosis history mapping
  diag_rec <- create_vax_map_branching(idx$n_vax, c(0, 1), i_eligible, i_p,
                                       set_vbe = FALSE, idx)

  # compartments to which vaccine efficacy applies
  ve <- c(rep(0, n_diag_rec), rep(1, n_erlang * n_diag_rec), rep(0, n_diag_rec))
  ved <- min(ved, 1 - 1e-10) # ensure duration is not divided by 0

  # create waning map
  if (stochastic == TRUE) {
    map <- create_waning_map_trial(n_vax, i_v, i_w, (n_erlang / dur))
    D <- diag(map)
    w <- sign(map)
  } else {
    w <- create_waning_map_trial(n_vax, i_v, i_w, (n_erlang / dur))
    D <- diag(w)
  }


  list(n_vax = n_vax,
    vea   = vea * ve,
    vei   = vei * ve,
    ved   = ved * ve,
    ves   = ves * ve,
    w     = w,
    D     = D,
    diag_rec = diag_rec
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
##' @param n_erlang integer giving the number of transitions that need to be
##'  made
##' @param stochastic logical indicating if the run should be made with the
##' default deterministic trial model in continuous time or stochastic trial
##' model in discrete time
##' @param N integer to assign the total number of individuals in the trial
##' (split equally across the two arms)
##' @param n_diag_rec integer giving the number of each X, V(erlang), and W
##' stratum, allowing tracking of diagnosis history. e.g for a n_diag_rec = 2
##' and erlang = 1, there will be X.I, X.II, V1.I, V1.II, W.I, W.II strata.
##' Where '.I' corresponds to never-diagnosed individuals and '.II' is for
##' individuals diagnosed at least once.
##' @inheritParams run_trial
##' @inheritParams vax_params_xvw_trial
##' @export


run_onevax_xvw_trial <- function(tt, gono_params, initial_params_trial = NULL,
                                 dur = 1e3,
                                 vea = 0, vei = 0, ved = 0, ves = 0,
                                 p_v = 0.5, n_erlang = 1,
                                 stochastic = FALSE,
                                 N = 6e05, n_diag_rec = 1) {

  stopifnot(all(lengths(list(vea, vei, ved, ves, dur)) %in%
                  c(1, length(gono_params))))
  assert_scalar_unit_interval(p_v)

  vax_params <- Map(vax_params_xvw_trial, dur = dur,
                    vea = vea, vei = vei, ved = ved, ves = ves,
                    n_erlang = n_erlang, stochastic = stochastic,
                    n_diag_rec = n_diag_rec)

  if (is.null(initial_params_trial)) {
    pars <- lapply(gono_params, model_params_trial, N = N,
                   n_diag_rec = n_diag_rec)
    init_params_trial <- Map(initial_params_xvw_trial, pars = pars,
                             p_v = p_v, n_erlang = n_erlang,
                             n_diag_rec = n_diag_rec)
  } else {
    init_params_trial <- initial_params_trial
  }

  ret <- Map(run_trial, gono_params = gono_params,
             init_params = init_params_trial,
             vax_params = vax_params,
             stochastic = stochastic,
             MoreArgs = list(tt = tt),
             N = N)

  # name outputs
  ret <- lapply(ret, name_outputs_trial, gen_trial_labels(n_erlang, n_diag_rec))
  ret

}


##' @name stratum_index_xvw_trial
##' @title Generate the indices of all xvw trial strata
##' @param n_erlang integer giving the number of transitions that need to be
##' made through vaccine-protected strata until that protection has waned
##' @param n_diag_rec integer giving the number of each X, V(erlang), and W
##' stratum, allowing tracking of diagnosis history. e.g for a n_diag_rec = 2
##' and erlang = 1, there will be X.I, X.II, V1.I, V1.II, W.I, W.II strata.
##' Where '.I' corresponds to never-diagnosed individuals and '.II' is for
##' individuals diagnosed at least once.
##' @return A list of strata with their indices
##' @export

stratum_index_xvw_trial <- function(n_erlang, n_diag_rec = 1) {

  # for an n_erlang of 3, and n_diag_rec of 2, the list of indexes returned
  # will be in the following order, where roman numerals refer to n_diag_rec,
  # and arabic numerals refer to erlang:
  # X.I, X.II, V1.I, V1.II, V2.I, V2.II, V3.I, V3.II, W.I, W.II

  ret <- list(X = seq_len(n_diag_rec))

  ret$V <- max(ret$X) + seq_len(n_erlang * n_diag_rec)
  ret$W <- max(ret$V) + seq_len(n_diag_rec)
  ret$n_vax <- max(ret$W)
  ret$n_erlang <- n_erlang
  ret$n_diag_rec <- n_diag_rec

  ret
}

##' @name gen_trial_labels
##' @title generates the appropriate strata labels for the number of strata
##' in the model, which depends on the value given to n_erlang and diagnosis
##' history levels desired (n_diag_rec)
##' @param n_erlang integer giving the number of transitions that need to be
##'  made through vaccine-protected strata until that protection has waned
##' @param n_diag_rec integer giving the number of levels of diagnosis history
##' for each X, V(*n_erlang), and W stratum
##' @return a character vector of length n_vax containing strata labels
##' @export
##' @importFrom utils as.roman
gen_trial_labels <- function(n_erlang = 1, n_diag_rec = 1) {


  diag_hist <- paste0(".", as.roman(seq_len(n_diag_rec)))

  output <- c(paste0("X", diag_hist),
              paste0("V", rep(seq_len(n_erlang), each = n_diag_rec), diag_hist),
              paste0("W", diag_hist))

  output

}

name_outputs_trial <- function(res, strata_names) {

  group_names <- c("L", "H")
  state_names <- c("U", "I", "A", "S", "T", "N",
                   "cum_incid", "cum_diag_a", "cum_diag_s",
                   "cum_treated", "cum_screened")

  for (nm in state_names) {

    dimnames(res[[nm]]) <- list(NULL, group_names, strata_names)
  }

  res
}
