### gonovax population model whereby individuals in X become vaccinated
# under the given vacc strategy with proportion r1 accepting 1 dose --> P
# or proportion r1r2 accepting both doses --> V
# individuals in P wane back to X, and can then go to P or V again after 
# individuals in V wane to W where they can get a booster vacc to R 
# rate of waning for P is 18 months
# this is a modificaiton of the onevax_xvwrh model where indiviuals pass
# linearly through the strata, now there is a branching after X

# also note to self may need to edit 'restart_hes' function in the onevax_xvwrh
# code so it works for this as well! rather than creating whole new function

##' Create initial conditions for the model
##' @name initial_params_xpvwrh
##' @title Initial conditions for the model
##' @inheritParams initial_params
##' @param hes proportion of population vaccine hesitant
##' @return A list of initial conditions
##' @export
initial_params_xpvwrh <- function(pars, coverage_p = 0, coverage_v = 0,
                                  hes = 0) {
  assert_scalar_unit_interval(coverage)
  n_vax <- 6
  
  willing <- 1 - hes
  x_init <- willing * (1 - coverage_p - coverage_v)
  p_init <- willing * coverage_p
  v_init <- willing * coverage_v
  cov <- c(x_init, p_init, v_init, 0, 0, hes)
  
  stopifnot(length(cov) == n_vax)
  stopifnot(sum(cov) == 1)
  
  U0 <- I0 <- A0 <- S0 <- T0 <- array(0, c(2, n_vax))
  
  # separate into 1:low and 2:high activity groups and by coverage
  N0 <- pars$N0 * outer(pars$q, cov)
  # set initial asymptomatic prevalence in each group (X AND H)
  A0[, 1] <- round(N0[, 1] * c(pars$prev_Asl, pars$prev_Ash))
  A0[, n_vax] <- round(N0[, n_vax] * c(pars$prev_Asl, pars$prev_Ash))
  
  # set initial uninfecteds
  U0 <- round(N0) - A0
  
  list(U0 = U0, I0 = I0, A0 = A0, S0 = S0, T0 = T0)
  
}