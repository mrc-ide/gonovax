# i = Activity group
## 1: Low
## 2: High
# j = Vaccination status
## 1: Unvaccinated
## 2: Bexsero
## 3: Waned

n_group <- 2
n_vax   <- user(1)

## adding timesteps
steps_per_year <- 365
dt <- 1 / steps_per_year
initial(time) <- 0
update(time) <- (step + 1) * dt

## assign low and high activity etas to the correct level
eta[1] <- eta_l
eta[2] <- eta_h

# individual probabilities of transitioning between infection states
p_UI[, ] <- 1 - exp(- (lambda * (1 - vea[j]) * dt))
p_ST[, ] <- 1 - exp(-mu * dt)
p_TU[, ] <- 1 - exp(-rho * dt)

p_A_or_S[, ] <- 1 - exp(-sigma * dt)

r_AT[, ] <- eta[i]
r_AU[, ] <- nu / (1 - ved[j])
p_T_or_U[, ] <- 1 - exp(- (r_AT[i, j] + r_AU[i, j] * dt))

# draws from binomial distributions for numbers changing between compartments
n_UI[, ] <- rbinom(U[i, j], p_UI[i, j])
n_ST[, ] <- rbinom(S[i, j], p_ST[i, j])
n_TU[, ] <- rbinom(T[i, j], p_TU[i, j])

n_IAS[, ] <- rbinom(I[i, j], p_A_or_S[i, j])
n_IA[, ] <- rbinom(n_IAS[i, j], (1 - ves[j]) * psi)
n_IS[, ] <- n_IAS[i, j] - n_IA[i, j]

n_AUT[, ] <- rbinom(A[i, j], p_T_or_U[i, j])
n_AT[, ] <- rbinom(n_AUT[i, j], r_AT[i, j] / (r_AT[i, j] + r_AU[i, j]))
n_AU[, ] <- n_AUT[i, j] - n_AT[i, j]
  
## Core equations for transitions between compartments:
update(U[, ]) <- U[i, j] - n_UI[i, j] +
  n_AU[i, j] + n_TU[i, j]  + sum(wU[i, j, ])

update(I[, ]) <- I[i, j] + n_UI[i, j] - n_IAS[i, j] + sum(wI[i, j, ])

update(A[, ]) <- A[i, j] + n_IA[i, j] -
               n_AUT[i, j] + sum(wA[i, j, ])

update(S[, ]) <- S[i, j] + n_IS[i, j] -
  n_ST[i, j]  + sum(wS[i, j, ])

update(T[, ]) <- T[i, j] + n_ST[i, j] + n_AT[i, j] - n_TU[i, j] +
  sum(wT[i, j, ])

## Update population size
N[, ] <- U[i, j] + I[i, j] + A[i, j] + S[i, j] + T[i, j]

screened[, ] <- eta[i] * U[i, j]
# vaccination -> no vaccination 'strategies' needed

# waning
n_Uw[, ] <- rbinom(U[i, j] - n_UI[i, j], 1 - exp(D[j] * dt))
n_Iw[, ] <- rbinom(I[i, j] - n_IAS[i, j], 1 - exp(D[j] * dt))
n_Aw[, ] <- rbinom(A[i, j] - n_AUT[i, j], 1 - exp(D[j] * dt))
n_Sw[, ] <- rbinom(S[i, j] - n_ST[i, j], 1 - exp(D[j] * dt))
n_Tw[, ] <- rbinom(T[i, j] - n_TU[i, j], 1 - exp(D[j] * dt))

wU[, , ] <- w[j, k] * n_Uw[i, k]
wI[, , ] <- w[j, k] * n_Iw[i, k]
wA[, , ] <- w[j, k] * n_Aw[i, k]
wS[, , ] <- w[j, k] * n_Sw[i, k]
wT[, , ] <- w[j, k] * n_Tw[i, k]

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

dim(n_Uw) <- c(n_group, n_vax)
dim(n_Iw) <- c(n_group, n_vax)
dim(n_Aw) <- c(n_group, n_vax)
dim(n_Sw) <- c(n_group, n_vax)
dim(n_Tw) <- c(n_group, n_vax)

dim(N)  <- c(n_group, n_vax)

dim(n_UI)     <- c(n_group, n_vax)
dim(n_IAS)    <- c(n_group, n_vax)
dim(n_IA)     <- c(n_group, n_vax)
dim(n_IS)     <- c(n_group, n_vax)
dim(n_AUT)    <- c(n_group, n_vax)
dim(n_AT)     <- c(n_group, n_vax)
dim(n_AU)     <- c(n_group, n_vax)
dim(n_ST)     <- c(n_group, n_vax)
dim(n_TU)     <- c(n_group, n_vax)
dim(screened) <- c(n_group, n_vax)

dim(p_UI)     <- c(n_group, n_vax)
dim(p_A_or_S) <- c(n_group, n_vax)
dim(p_ST)     <- c(n_group, n_vax)
dim(p_T_or_U) <- c(n_group, n_vax)
dim(p_TU)     <- c(n_group, n_vax)

dim(r_AT) <- c(n_group, n_vax)
dim(r_AU) <- c(n_group, n_vax)

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
D[] <- user()

## par dimensions

dim(eta)  <- n_group
dim(vea)  <- n_vax
dim(ved)  <- n_vax
dim(ves)  <- n_vax

dim(w)    <- c(n_vax, n_vax)
dim(D)    <- c(n_vax)

dim(wU)   <- c(n_group, n_vax, n_vax)
dim(wI)   <- c(n_group, n_vax, n_vax)
dim(wA)   <- c(n_group, n_vax, n_vax)
dim(wS)   <- c(n_group, n_vax, n_vax)
dim(wT)   <- c(n_group, n_vax, n_vax)

output(N)   <- N
