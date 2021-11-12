
## imports the gonorrhoea fitted NHPs and transforms into list format
gono_params_trial <- function(n = NULL) {
  if (is.null(cache$gono_params_trial)) {
    gp <- read_csv(gonovax_file("extdata/gono_params_updated.csv"))
    gp$lambda <- 1.5   #come back to this later!! Your lambda will be different
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
  check_gono_params(pars)
  with(pars, {
    assert_scalar_positive(beta)
    assert_scalar_positive(eta_h)
    assert_scalar_positive(eta_l)
  })


  pars$tt <- c(0, 1e3)

  pars

}


## sets up total trial size and the proportion of this moving to high activity

demographic_params_trial <- function() {
  list(N0 = 6e5,
       move = c(0, 1)
  )
}



## sets up the starting conditions for state variables of each activity group

initial_params_trial <- function(pars, n_vax = 1, coverage = 1) {

  stopifnot(length(coverage) == n_vax)
  stopifnot(sum(coverage) == 1)

  U0 <- I0 <- A0 <- S0 <- T0 <- array(0, c(2, n_vax))

  # separate into 1:low and 2:high activity groups and by coverage
  N0 <- pars$N0 * outer(pars$move, coverage)

  # put the vaccinated and placebo individuals all to uninfected
  U0 <- round(N0)

  list(U0 = U0, I0 = I0, A0 = A0, S0 = S0, T0 = T0)
}



### generates vector of all starting conditions

model_params_trial <- function(gono_params_trial = NULL,
                        demographic_params_trial = NULL,
                        initial_params_trial = NULL,
                        vax_params = NULL, coverage = 0) {
  gono_params_trial <- gono_params_trial %||% gono_params_trial(1)[[1]]
  dpt <- demographic_params_trial  %||% demographic_params_trial()
  demographic_params_trial <- dpt
  ret <- c(demographic_params_trial, gono_params_trial)
  vax_params <- vax_params %||% vax_params0()

  if (coverage == 0) {
    cov <- c(1, rep(0, vax_params$n_vax - 1))
    initial_params <-
      initial_params_trial %||% initial_params_trial(ret, vax_params$n_vax, cov)

  } else {
    initial_params <- initial_params_xvw_trial(pars = ret, coverage = coverage)
  }

  c(ret, initial_params, vax_params)
}
