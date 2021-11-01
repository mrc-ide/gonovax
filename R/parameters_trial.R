

## gonoparams , the oG to tweak:
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
  
  
  pars$tt <- c(0, 1e3)                 #adds list element of time 
 # pars$beta <-  rep(pars$beta, 2)      #adds a second dimension 
 # pars$eta_l <- rep(pars$eta_l, 2)
 # pars$eta_h <- rep(pars$eta_h, 2)
 #1 pars$beta <- pars$eta_h <- NULL           #not too sure what this bit does     why go to the effort of checking just to make it null? 
                                            #come back to it!
  
  pars
}
