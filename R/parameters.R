demographic_params <- function() {
  list(N0 = 6e5,
       q = c(0.85, 0.15), # proportion in group L
       p = c(0.6, 15.6), # partner change rate in group L/H
       enr = 12000,
       exr = 1 / 50
  )
}

##' @name gono_params
##' @title Posterior parameters of gonorrhoea natural history
##' @param n an integer vector (or value) containing the indices of the required
##' parameter sets (1:982). If `n = NULL` the full parameter set is returned
##' @return A list of parameters
##' @export

gono_params <- function(n = NULL) {
  if (is.null(cache$gono_params)) {
    gp <- read_csv(gonovax_file("extdata/gono_params.csv"))
    cache$gono_params <-
      lapply(seq_len(nrow(gp)), function(i) transform0(gp[i, ]))
  }

  pars <- cache$gono_params
  n_pars <- length(pars)
  # if n not supplied, return all parameters
  i  <- n %||% seq_len(n_pars)
  # limit to parameter sets available
  i <- i[(i > 0) & (i <= n_pars)]

  pars[i]
}

transform0 <- function(pars) {
    # reformat time-varying pars
    pars <- as.list(pars)

    check_gono_params(pars)
    with(pars, {
      assert_scalar_positive(beta)
      assert_scalar_positive(eta)
    })


    pars$tt <- c(0, 1e3)
    pars$beta_t <-  rep(pars$beta, 2)
    pars$eta_l_t <- rep(pars$eta, 2)
    pars$eta_h_t <- pars$eta_l_t
    pars$beta <- pars$eta <- NULL

    pars
}

##' Transform fitted parameters into gonovax params
##' @name transform
##' @title Transform fitted parameters into gonovax params
##' @param pars list of fitted parameters
##' @param fix_par_t logical indicating whether parameters with inferred trends
##' are held at their 2020 values after that date.
##' @return A list of parameters for use in the model
##' @export
transform <- function(pars, fix_par_t = TRUE) {
  # reformat time-varying pars
  pars <- as.list(pars)
  check_gono_params(pars)
  with(pars, {
    assert_scalar_positive(beta2009)
    assert_scalar_positive(phi_beta)
    assert_scalar_positive(eta_h)
    assert_scalar_positive(phi_eta)
    assert_scalar_unit_interval(omega)
  })

  t0 <- 2009
  t1 <- 2020
  t_max <- 1e2

  t_fix <- ifelse(fix_par_t, t1, t0 + t_max)
  pars$tt <- gonovax::gonovax_year(pmin(seq(t0, t0 + t_max), t_fix))
  pars$beta_t  <- pars$beta2009 * (1 + pars$phi_beta * pars$tt)

  pars$eta_l_t <- pars$eta_h * pars$omega * (1 + pars$phi_eta * pars$tt)
  pars$eta_h_t <- pars$eta_h * (1 + pars$phi_eta * pars$tt)

  pars$tt <- seq(0, t_max)

  pars
}

##'
##' @name transform_fixed
##' @title  Transform fitted parameters into non-time-varying gonovax params
##' @param pars list of fitted parameters
##' @return A list of parameters for use in the model
##' @export
transform_fixed <- function(pars) {
  # reformat time-varying pars
  pars <- as.list(pars)
  check_gono_params(pars)
  with(pars, {
    assert_scalar_positive(beta)
    assert_scalar_positive(eta_l)
    assert_scalar_positive(eta_h)
  })

  pars$tt <- c(0, 1e3)
  pars$beta_t <-  rep(pars$beta, 2)
  pars$eta_l_t <- rep(pars$eta_l, 2)
  pars$eta_h_t <- rep(pars$eta_h, 2)
  pars$beta <- pars$eta_l <- pars$eta_h <- NULL

  pars
}

##' Create initial conditions for the model
##' @name initial_params
##' @title Initial conditions for the model
##' @param pars A parameter list containing `N0`, `q`, `prev_Asl` and `prev_Ash`
##'   elements.
##' @param n_vax an integer indicating the number of vaccine compartments
##' @param coverage a vector of length `n_vax` that sums to 1 denoting the
##' initial proportion in each vaccine stratum
##' @return A list of initial model states
##' @export
initial_params <- function(pars, n_vax = 1, coverage = 1) {

  stopifnot(length(coverage) == n_vax)
  stopifnot(sum(coverage) == 1)

  U0 <- I0 <- A0 <- S0 <- T0 <- array(0, c(2, n_vax))

  # separate into 1:low and 2:high activity groups and by coverage
  N0 <- pars$N0 * outer(pars$q, coverage)
  # set initial asymptomatic prevalence in each group (unvaccinated only)
  A0[, 1] <- round(N0[, 1] * c(pars$prev_Asl, pars$prev_Ash))

  # set initial uninfecteds
  U0 <- round(N0) - A0

  list(U0 = U0, I0 = I0, A0 = A0, S0 = S0, T0 = T0)
}

##' Create initial conditions based on previous model run
##' @name restart_params
##' @title Create initial conditions to start the model from the end of a run
##' @param y a transformed model run output
##' @param n_vax an integer indicating the number of vaccine compartments,
##' consistent with the input
##' @return A list of initial conditions to restart a model with n_vax
##' vaccination levels
##' @export
restart_params <- function(y, n_vax = NULL) {
  dim_y <- dim(y[["U"]])

  i_t <- dim_y[1]
  n_vax <- n_vax %||% dim_y[3]

  n_vax_input <- dim_y[3]
  i_vax <- seq_len(min(n_vax,  n_vax_input))

  U0 <- I0 <- A0 <- S0 <- T0 <- array(0, c(2, n_vax))
  # set compartments in each group
  U0[, i_vax] <- y$U[i_t, , i_vax]
  I0[, i_vax] <- y$I[i_t, , i_vax]
  A0[, i_vax] <- y$A[i_t, , i_vax]
  S0[, i_vax] <- y$S[i_t, , i_vax]
  T0[, i_vax] <- y$T[i_t, , i_vax]

  list(U0 = U0, I0 = I0, A0 = A0, S0 = S0, T0 = T0, t = y$t[i_t])
}

##' @name model_params
##' @title Parameters for the dualvax model
##' @param gono_params A dataframe of natural history parameters
##' @param demographic_params A dataframe of demographic parameters
##' @param vax_params A vector of vaccination params
##' @param init_params A list of starting conditions
##' @return A list of inputs to the model many of which are fixed and
##'   represent data. These correspond largely to `user()` calls
##'   within the odin code, though some are also used in processing
##'   just before the model is run.
##' @export
model_params <- function(gono_params = NULL,
                         demographic_params = NULL,
                         init_params = NULL,
                         vax_params = NULL) {
  gono_params <- gono_params %||% gono_params(1)[[1]]
  demographic_params <- demographic_params %||% demographic_params()
  ret <- c(demographic_params, gono_params)
  vax_params <- vax_params %||% vax_params0()

  cov <- c(1, rep(0, vax_params$n_vax - 1))
  init_params <- init_params %||% initial_params(ret, vax_params$n_vax, cov)
  c(ret, init_params, vax_params)
}

##' @name create_vax_map
##' @title Create mapping for movement between strata due to vaccination
##' @param n_vax Integer in (0, 5) denoting total number of strata
##' @param v a matrix of dimensions 2x2 (vod, vos) or a scalar (vbe) indicating
##' % of population vaccinated according to the strategy. For vod & vos each
##' value will apply to the eligible stratum separately
##' @param i_u indices of strata eligible for vaccination
##' @param i_v indices of strata being vaccinated
##' @return an array of the mapping

create_vax_map <- function(n_vax, v, i_u, i_v) {

  # ensure vaccine input is of correct length
  n_group <- 2
  stopifnot(all((v >= 0) & (v <= 1)))
  stopifnot(length(i_v) == length(i_u))
  stopifnot(max(i_u, i_v) <= n_vax)
  stopifnot(all(dim(v) == c(n_group, length(i_v))))

  # set up vaccination matrix
  vax_map <- array(0, dim = c(n_group, n_vax, n_vax))

  for (i in seq_along(i_u)) {
    vax_map[, i_u[i], i_u[i]] <-  v[i, ]
    vax_map[, i_v[i], i_u[i]] <- -v[i, ]
  }

  vax_map

}

##' @name create_waning_map
##' @title Create mapping for movement between strata due to vaccine waning
##' @param n_vax Integer in (0, 5) denoting total number of strata
##' @param i_v indices of strata being vaccinated
##' @param i_w Integer in (0, 5) denoting which stratum receives waned vaccinees
##' @param z Scalar denoting rate of waning
##' @return an array of the mapping

create_waning_map <- function(n_vax, i_v, i_w, z) {

  stopifnot(z > 0)
  stopifnot(length(z) %in% c(1, length(i_v)))
  stopifnot(length(i_w) == 1)
  # set up waning map
  w <- array(0, dim = c(n_vax, n_vax))

  w[i_w, i_v] <-  z

  for (i in i_v) {
    w[i, i] <- -w[i_w, i]
  }
  w
}


set_strategy <- function(strategy, primary_uptake,
                         booster_uptake = primary_uptake) {

  if (length(primary_uptake) != 1) {
    stop("primary vaccination uptake must be length 1")
  }

  if (length(booster_uptake) != 1) {
    stop("booster vaccination uptake must be length 1")
  }

uptake <- c(primary_uptake, booster_uptake)

  if (strategy == "VbE") {
    vos <- vod <- matrix(0, 2, 2)
  } else if (strategy == "VoD") {
    vod <- cbind(uptake, uptake, deparse.level = 0)
    vos <- matrix(0, 2, 2)
  } else if (strategy == "VoA") {
    vod <- vos <- cbind(uptake, uptake, deparse.level = 0)
  } else if (strategy == "VoD(H)") {
    vod <- cbind(0, uptake, deparse.level = 0)
    vos <- matrix(0, 2, 2)
  } else if (strategy == "VoA(H)") {
    vod <- vos <- cbind(0, uptake, deparse.level = 0)
  } else if (strategy == "VoD(L)+VoA(H)") {
    vod <- cbind(uptake, uptake, deparse.level = 0)
    vos <- cbind(0, uptake, deparse.level = 0)
  } else if (strategy == "VoS") {
    vod <- matrix(0, 2, 2)
    vos <- cbind(uptake, uptake, deparse.level = 0)
  } else {
    stop("strategy not recognised")
  }

  list(vod = vod, vos = vos)
}

check_gono_params <- function(pars) {
  with(pars, {
    assert_scalar_unit_interval(psi)
    assert_scalar_unit_interval(prev_Asl)
    assert_scalar_unit_interval(prev_Ash)
    assert_scalar_unit_interval(epsilon)
    assert_scalar_positive(sigma)
    assert_scalar_positive(nu)
    assert_scalar_positive(mu)
    assert_scalar_positive(rho)
  })
}
