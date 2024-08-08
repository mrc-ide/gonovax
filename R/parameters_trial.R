##' @name gono_params_trial
##' @title Posterior parameters of gonorrhoea natural history
##' @param n an integer vector (or value) containing the indices of the required
##' parameter sets (1:982). If `n = NULL` the full parameter set is returned
##' @return A list of parameters
##' @export

gono_params_trial <- function(n = NULL) {
  if (is.null(cache$gono_params_trial)) {
    gp <- read_csv(gonovax_file("extdata/gono_params_trial.csv"))
    cache$gono_params_trial <-
      lapply(seq_len(nrow(gp)), function(i) transform0_trial(gp[i, ]))
  }

  pars <- cache$gono_params_trial
  n_pars <- length(pars)
  # if n not supplied, return all parameters
  i  <- n %||% seq_len(n_pars)
  # limit to parameter sets available
  i <- i[(i > 0) & (i <= n_pars)]

  pars[i]

}



## transforms the imported parameters into workable format

transform0_trial <- function(pars) {
  pars <- as.list(pars)              #converts each row of csv into list
  check_gono_params_trial(pars)
  with(pars, {
    assert_scalar_positive(eta)
  })

  pars

}

check_gono_params_trial <- function(pars) {
  with(pars, {
    assert_scalar_unit_interval(psi)
    assert_scalar_unit_interval(epsilon)
    assert_scalar_positive(sigma)
    assert_scalar_positive(nu)
    assert_scalar_positive(mu)
    assert_scalar_positive(rho)
    assert_scalar_positive(lambda)
  })
}


## sets up total trial size and the proportion that are in the high activity
## group

demographic_params_trial <- function(N = 6e5) {
  list(N0 = N,
    q = c(0, 1)
  )
}


##' Create initial conditions for the model trial
##' @name initial_params_trial
##' @title Initial conditions for the model trial where the entire cohort is
##' in the high sexual activity group.
##' @param pars A parameter list containing `N0`, and `q` elements.
##' @param n_vax an integer indicating the number of vaccine compartments
##' @param p_v a vector of length `n_vax` that sums to 1 denoting the
##' proportion in each vaccine stratum
##' @param n_diag_rec integer giving the number of each X, V(erlang), and W
##' stratum, allowing tracking of diagnosis history. e.g for a n_diag_rec = 2
##' and erlang = 1, there will be X.I, X.II, V1.I, V1.II, W.I, W.II strata.
##' Where '.I' corresponds to never-diagnosed individuals and '.II' is for
##' individuals diagnosed at least once.
##' @return A list of initial model states
##' @export

initial_params_trial <- function(pars, n_vax = 1, p_v = 1,
                                 n_diag_rec = n_diag_rec) {

  stopifnot(length(p_v) == n_vax)
  stopifnot(sum(p_v) == 1)

  U0 <- I0 <- A0 <- S0 <- T0 <- array(0, c(2, n_vax))

  # separate into 1:low and 2:high activity groups and by p_v
  N0 <- pars$N0 * outer(pars$q, p_v)

  # put the vaccinated and placebo individuals all to uninfected
  U0 <- round(N0)

  list(U0 = U0, I0 = I0, A0 = A0, S0 = S0, T0 = T0)
}


### generates list of all parameters needed to run the model

##' @name model_params_trial
##' @title Parameters for the vaccination trial model
##' @param gono_params_trial A dataframe of natural history parameters
##' @param demographic_params_trial A dataframe of demographic parameters
##' @param vax_params A vector of vaccination params
##' @param p_v A scalar indicating the percentage of the trial cohort that
##'  is vaccinated
##' @param initial_params_trial A list of starting conditions
##' @param n_erlang scalar giving the number of transitions that need to be made
##' through vaccine-protected strata until that protection has waned
##' @param N integer to assign the total number of individuals in the trial
##' (split equally across the two arms)
##' @param n_diag_rec integer giving the number of each X, V(erlang), and W
##' stratum, allowing tracking of diagnosis history. e.g for a n_diag_rec = 2
##' and erlang = 1, there will be X.I, X.II, V1.I, V1.II, W.I, W.II strata.
##' Where '.I' corresponds to never-diagnosed individuals and '.II' is for
##' individuals diagnosed at least once.
##' @param asymp_recorded logical indicating if the trial screens for and
##' records asymptomatic diagnosis. If FALSE, asymptomatic infected individuals
##' undergoing treatment do not move diagnosis history stratum
##' @return A list of inputs to the model many of which are fixed and
##'   represent data. These correspond largely to `user()` calls
##'   within the odin code, though some are also used in processing
##'   just before the model is run.
##' @export

model_params_trial <- function(gono_params_trial = NULL,
                               demographic_params_trial = NULL,
                               initial_params_trial = NULL,
                               vax_params = NULL, p_v = 0,
                               n_erlang = 1, N = 6e5,
                               n_diag_rec = 1,
                               asymp_recorded = TRUE) {
  gono_params_trial <- gono_params_trial %||% gono_params_trial(1)[[1]]
  demographic_params_trial <-
    demographic_params_trial  %||% demographic_params_trial(N = N)
  ret <- c(demographic_params_trial, gono_params_trial)

  #check n_erlang supplied in model_params_trial() is same as
  #n_erlang supplied to vax_params_xvw_trial()
  #unless vax_params not supplied
  if (is.null(vax_params) == FALSE) {  #evaluates to TRUE if vax_params supplied
    stopifnot(unique(dim(vax_params$w)) ==  (2 + n_erlang) * n_diag_rec)
  } else {

    #also add in diag_rec if vax_params not supplied
    vax_params <- vax_params0(n_diag_rec = n_diag_rec)
    n_vax <- vax_params$n_vax
    i_eligible <-  if (n_diag_rec == 1) {
      i_eligible <- 0
    } else {
      i_eligible <- seq_len(n_vax)[seq_len(n_vax) %% n_diag_rec != 0]
    }

    i_vaccinees <- seq_len(n_vax)[seq_len(n_vax) %% n_diag_rec != 1]
    #diagnosis history
    #symptomatic
    vax_params$diag_rec_s <- create_vax_map(n_vax, c(1, 1), i_eligible,
                                            i_vaccinees)

    #asymptomatic
    if (asymp_recorded) {
      vax_params$diag_rec_a <- vax_params$diag_rec_s
    } else {
      vax_params$diag_rec_a <- create_vax_map(n_vax, c(0, 0), i_eligible,
                                              i_vaccinees)
    }


  }

  #passing initial parameters
  if (p_v == 0) {                             #no vaccination = placebo arm
    cov <- c(1, rep(0, vax_params$n_vax - 1))

    initial_params <-
      initial_params_trial %||% initial_params_trial(ret, vax_params$n_vax, cov,
                                                     n_diag_rec = n_diag_rec)
    #if no vaccination then no need for V and W strata

  } else {                                  #p_v greater than 0 = vaccinated arm
    initial_params <- initial_params_xvw_trial(pars = ret, p_v = p_v,
                                               n_erlang = n_erlang,
                                               n_diag_rec = n_diag_rec)
  }

  c(ret, initial_params, vax_params)
}

##' @name create_waning_map_trial
##' @title Create mapping for movement between strata due to vaccine waning
##' in a vaccine trial with erlang compartments
##' @param n_vax Integer in (0, 5) denoting total number of strata
##' @param i_v indices of strata under vaccination protection
##' @param i_w indicies denoting which stratum receives waned vaccinees
##' @param z Scalar denoting rate of waning
##' @return an array of the mapping

create_waning_map_trial <- function(n_vax, i_v, i_w, z) {

  stopifnot(z >= 0)

  # set up waning map
  w <- array(0, dim = c(n_vax, n_vax))

  for (i in seq_along(i_v)) {
    w[i_w[i], i_v[i]] <- z
  }

  for (i in i_v) {
    idx <- which(i_v == i)

    w[i, i] <- -w[i_w[idx], i]
  }

  w

}
