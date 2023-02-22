# i = Activity group
## 1: Low
## 2: High
# j = Vaccination status
## 1: Unvaccinated
## 2: Bexsero
## 3: Waned

n_group <- 2
n_vax   <- user(1)

## assign low and high activity etas to the correct level
eta[1] <- eta_l
eta[2] <- eta_h

## Core equations for transitions between compartments:

update(U[, ]) <- U[i, j] - n_UI[i, j] +
  n_AU[i, j] + n_TU[i, j]  + sum(wU[i, j, ])

update(I[, ]) <- I[i,j] + n_UI[i, j] - sigma * I[i, j] + sum(wI[i, j, ])

update(A[, ]) <- A[i,j] +(1 - (1 - ves[j]) * psi) * sigma * I[i, j] 
  - n_AT[i, j] - n_AU[i, j] + sum(wA[i, j, ])

update(S[, ]) <- S[i,j] + (1 - ves[j]) * psi * sigma * I[i, j] -
  n_ST[i, j]  + sum(wS[i, j, ])

update(T[, ]) <- T[i,j] + n_ST[i, j] + n_AT[i, j]  - n_TU[i, j] +
  sum(wT[i, j, ])

## Update population size
N[, ] <- U[i, j] + I[i, j] + A[i, j] + S[i, j] + T[i, j]

n_UI[, ]     <- lambda * (1 - vea[j]) * U[i, j]  # force of infection constant
n_AT[, ]     <- eta[i] * A[i, j]                 # in trial model
n_AU[, ]     <- nu / (1 - ved[j]) * A[i, j]
n_ST[, ]     <- mu * S[i, j]
n_TU[, ]     <- rho * T[i, j]
screened[, ] <- eta[i] * U[i, j]

# vaccination -> no vaccination 'strategies' needed

# waning
wU[, , ] <- w[j, k] * U[i, k]
wI[, , ] <- w[j, k] * I[i, k]
wA[, , ] <- w[j, k] * A[i, k]
wS[, , ] <- w[j, k] * S[i, k]
wT[, , ] <- w[j, k] * T[i, k]

## outputs

update(cum_incid[, ])      <- n_UI[i, j]
update(cum_diag_a[, ])     <- n_AT[i, j]
update(cum_diag_s[, ])     <- n_ST[i, j]
update(cum_treated[, ])    <- n_TU[i, j]
update(cum_screened[, ])   <- screened[i, j]

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

## Parameters

eta_l     <- user()
eta_h     <- user()

sigma     <- user()
psi       <- user()
nu        <- user()
mu        <- user()
rho       <- user()
lambda    <- user()

# vaccine effects
vea[] <- user() # efficacy against acquisition
ved[] <- user() # efficacy against duration of infection
ves[] <- user() # efficacy against symptoms

w[, ]    <- user()

## par dimensions

dim(eta)  <- n_group
dim(vea)  <- n_vax
dim(ved)  <- n_vax
dim(ves)  <- n_vax

dim(w)    <- c(n_vax, n_vax)

dim(wU)   <- c(n_group, n_vax, n_vax)
dim(wI)   <- c(n_group, n_vax, n_vax)
dim(wA)   <- c(n_group, n_vax, n_vax)
dim(wS)   <- c(n_group, n_vax, n_vax)
dim(wT)   <- c(n_group, n_vax, n_vax)

output(N)   <- N