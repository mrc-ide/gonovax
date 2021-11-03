

## gonoparams , the oG to tweak :
gono_params_trial <- function(n = NULL) {
  if (is.null(cache$gono_params_trial)) {
    gp <- read_csv(gonovax_file("extdata/gono_params_updated.csv"))
    cache$gono_params_trial <-
      lapply(seq_len(nrow(gp)), function(i) transform0_trial(gp[i, ]))    #returns a 1000 element-length list with each element transformed 
  }
  
  pars <- cache$gono_params_trial
  n_pars <- length(pars)
  # if n not supplied, return all parameters
  i  <- n %||% seq_len(n_pars)
  # limit to parameter sets available
  i <- i[(i > 0) & (i <= n_pars)]
  
  pars[i]
  
  }


#edit transform function also

transform0_trial <- function(pars) {          #pars = each individual row of all columns of the csv
  # reformat time-varying pars     #note, no longer time varying
  
  pars <- as.list(pars)          #makes each row of csv a list 
  
  check_gono_params(pars)            #checking they're positive scalars
  with(pars, {
    assert_scalar_positive(beta)
    assert_scalar_positive(eta_h)
    assert_scalar_positive(eta_l)
  })
  
  
  pars$tt <- c(0, 1e3)        #delete? 

  pars
}

### next step set up demographics correctly 

demographic_params_trial <- function() {
  list(N0 = 6e5,
       move = c(0, 1)
  )
}



initial_params_trial <- function(pars, n_vax = 1, coverage = 1) {
  
  stopifnot(length(coverage) == n_vax)
  stopifnot(sum(coverage) == 1)
  
  U0 <- I0 <- A0 <- S0 <- T0 <- array(0, c(2, n_vax))
  # separate into 1:low and 2:high activity groups and by coverage
  N0 <- pars$N0 * outer(pars$move, coverage)
  
#   set initial asymptomatic prevalence in each group (unvaccinated only)  #seeding infections 
# A0[, 1] <- round(N0[, 1] * c(pars$prev_Asl, pars$prev_Ash))                                   #lambda is constant and we don't need to seed infections
  
  # set initial uninfecteds
  U0 <- round(N0) - A0
  
  list(U0 = U0, I0 = I0, A0 = A0, S0 = S0, T0 = T0)
}



####           #changin for init_params_trial or demographic_params_trial won't work! 

model_params_trial <- function(gono_params_trial = NULL,
                        demographic_params_trial = NULL,
                        initial_params_trial = NULL,
                        vax_params = NULL) {
  gono_params_trial <- gono_params_trial %||% gono_params_trial(1)[[1]]
  demographic_params_trial <- demographic_params_trial  %||% demographic_params_trial()
  ret <- c(demographic_params_trial, gono_params_trial)
  vax_params <- vax_params %||% vax_params0()                                               #vax_params0() is the conditions for no vaccination
  
  cov <- c(1, rep(0, vax_params$n_vax - 1))
  initial_params_trial <- initial_params_trial%||% initial_params_trial(ret, vax_params$n_vax, cov)
  c(ret, initial_params_trial, vax_params)
}


#############