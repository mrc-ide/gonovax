# i = Activity group 
## 1: Low
## 2: High
# j = Vaccination status
## 1: Unvaccinated
## 2: Bexsero
## 3: Waned

n_group <- 2                
n_vax   <- user(1)

## calibrate time-varying parameters
# tt runs from t0 = 2009, to t10 = 2019
tt[] <- user()
dim(tt) <- user()            #have left tt here but don't think relevant? 
eta[1] <- eta_l                                     
eta[2] <- eta_h

## Core equations for transitions between compartments:

deriv(U[, ]) <-  - n_UI[i, j] + #- exr * U[i, j] +
  n_AU[i, j] + n_TU[i, j]  + sum(wU[i, j, ]) #-
#sum(n_vbe[i, j, ]) - sum(n_vod[i, j, ]) - sum(n_vos[i, j, ]) #+ entrants[i, j]

#deriv(I[, ]) <- n_UI[i, j] - (sigma + exr) * I[i, j] + sum(wI[i, j, ])
deriv(I[, ]) <- n_UI[i, j] - (sigma) * I[i, j] + sum(wI[i, j, ])

deriv(A[, ]) <- (1 - (1 - ves[j]) * psi) * sigma * I[i, j] - n_AT[i, j] -
 # n_AU[i, j] - exr * A[i, j] + sum(wA[i, j, ])
  n_AU[i, j] + sum(wA[i, j, ])
  
deriv(S[, ]) <- (1 - ves[j]) * psi * sigma * I[i, j] - n_ST[i, j]  + sum(wS[i, j, ])
# -  exr * S[i, j]

  
#deriv(T[, ]) <- n_ST[i, j] + n_AT[i, j] - exr * T[i, j] - n_TU[i, j] +
deriv(T[, ]) <- n_ST[i, j] + n_AT[i, j]  - n_TU[i, j] +
  sum(wT[i, j, ])

## Update population size
N[, ] <- U[i, j] + I[i, j] + A[i, j] + S[i, j] + T[i, j]
#entrants[, 1] <- enr  #* q[i]


n_UI[, ]     <- lambda * (1 - vea[j]) * U[i, j]  #changed lambda[i] to constant 
n_AT[, ]     <- eta[i] * A[i, j]
n_AU[, ]     <- nu / (1 - ved[j]) * A[i, j]
n_ST[, ]     <- mu * S[i, j]
n_TU[, ]     <- rho * T[i, j]
screened[, ] <- eta[i] * U[i, j]

# vaccination --> no need for this 
## time-varying switch
#vax_switch <- interpolate(vax_t, vax_y, "constant")
## at screening
#n_vos[, , ] <- vos[i, j, k] * screened[i, k] * vax_switch
## on diagnosis
#n_vod[, , ] <- vod[i, j, k] * n_TU[i, k] * vax_switch
## on entry - no switch as background rate
#n_vbe[, , ] <- vbe[i, j, k] * entrants[i, k]

# waning
wU[, , ] <- w[j, k] * U[i, k]
wI[, , ] <- w[j, k] * I[i, k]
wA[, , ] <- w[j, k] * A[i, k]
wS[, , ] <- w[j, k] * S[i, k]
wT[, , ] <- w[j, k] * T[i, k]

## outputs

deriv(cum_incid[, ])      <- n_UI[i, j]
deriv(cum_diag_a[, ])     <- n_AT[i, j]
deriv(cum_diag_s[, ])     <- n_ST[i, j]
deriv(cum_treated[, ])    <- n_TU[i, j]
deriv(cum_screened[, ])   <- screened[i, j]
#deriv(cum_vaccinated[, ]) <- n_vos[i, j, j] + n_vod[i, j, j] + n_vbe[i, j, j]
#deriv(cum_vbe[, ]) <- n_vbe[i, j, j]

# aggregated time series for fitting mcmc
output(tot_treated) <- sum(cum_treated)
output(tot_attended) <- sum(cum_treated) + sum(cum_screened)

## Set up compartments
## Initial states are all 0 as we will provide a state vbector
initial(U[, ]) <- U0[i, j]
initial(I[, ]) <- I0[i, j]
initial(A[, ]) <- A0[i, j]
initial(S[, ]) <- S0[i, j]
initial(T[, ]) <- T0[i, j]

U0[, ] <- user()
I0[, ] <- user()
A0[, ] <- user()
S0[, ] <- user()
T0[, ] <- user()

initial(cum_incid[, ])      <- 0
initial(cum_diag_a[, ])     <- 0
initial(cum_diag_s[, ])     <- 0
initial(cum_treated[, ])    <- 0
initial(cum_screened[, ])   <- 0
#initial(cum_vaccinated[, ]) <- 0
#initial(cum_vbe[, ]) <- 0

# set up dimensions of compartments
dim(U) <- c(n_group, n_vax)
dim(I) <- c(n_group, n_vax)
dim(A) <- c(n_group, n_vax)
dim(S) <- c(n_group, n_vax)
dim(T) <- c(n_group, n_vax)

dim(U0) <- c(n_group, n_vax)
dim(I0) <- c(n_group, n_vax)
dim(A0) <- c(n_group, n_vax)
dim(S0) <- c(n_group, n_vax)
dim(T0) <- c(n_group, n_vax)

dim(N)  <- c(n_group, n_vax)
#dim(entrants) <- c(n_group, n_vax)

dim(n_UI)     <- c(n_group, n_vax)
dim(n_AT)     <- c(n_group, n_vax)
dim(n_AU)     <- c(n_group, n_vax)
dim(n_ST)     <- c(n_group, n_vax)
dim(n_TU)     <- c(n_group, n_vax)
dim(screened) <- c(n_group, n_vax)

dim(cum_incid)      <- c(n_group, n_vax)
dim(cum_diag_a)     <- c(n_group, n_vax)
dim(cum_diag_s)     <- c(n_group, n_vax)
dim(cum_treated)    <- c(n_group, n_vax)
dim(cum_screened)   <- c(n_group, n_vax)
#dim(cum_vaccinated) <- c(n_group, n_vax)
#dim(cum_vbe)        <- c(n_group, n_vax)

## Parameters
#move[]      <- user() taking out temporarily 

eta_l     <- user()
eta_h     <- user() 

#enr       <- user()
#exr       <- user()

sigma     <- user()
psi       <- user()
nu        <- user()
mu        <- user()
rho       <- user()
lambda    <- user()         

## vaccination parameters
# vaccination routes
#vbe[, , ] <- user()
#vos[, , ] <- user()
#vod[, , ] <- user()

# vaccine effects
vea[] <- user() # efficacy against acquisition
ved[] <- user() # efficacy against duration of infection
ves[] <- user() # efficacy against symptoms

w[, ]    <- user()
vax_t[]  <- user()
vax_y[]  <- user()

## par dimensions
#dim(move) <- n_group   taking out temporarily 

dim(eta)  <- n_group
dim(vea)  <- n_vax
dim(ved)  <- n_vax
dim(ves)  <- n_vax

#dim(vbe)   <- c(n_group, n_vax, n_vax)
#dim(vod)   <- c(n_group, n_vax, n_vax)
#dim(vos)   <- c(n_group, n_vax, n_vax)
dim(w)    <- c(n_vax, n_vax)
dim(vax_t) <- user()
dim(vax_y) <- user()

#dim(n_vbe) <- c(n_group, n_vax, n_vax)
#dim(n_vos) <- c(n_group, n_vax, n_vax)
#dim(n_vod) <- c(n_group, n_vax, n_vax)
dim(wU)   <- c(n_group, n_vax, n_vax)
dim(wI)   <- c(n_group, n_vax, n_vax)
dim(wA)   <- c(n_group, n_vax, n_vax)
dim(wS)   <- c(n_group, n_vax, n_vax)
dim(wT)   <- c(n_group, n_vax, n_vax)

output(N)   <- N
output(lambda) <- lambda
