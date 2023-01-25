# Key to indices
# i = Activity group
## 1: Low
## 2: High
# j = Vaccination status
## 1: Unvaccinated (X)
## Then optional from the below, depending on model structure
## 2: Vaccinated (V)
## 3: Waned (W)
## 4: Revaccinated (R)
## 5: Hesitant (H)
## When a branching model is used such as XPVWRH (2:P, 3:v, 4:W, 5:R, 6:H)

n_group <- 2
n_vax   <- user(1)
n_erlang <- 1         #set to 1 for now! 

## calibrate time-varying parameters
# tt runs from t0 = 2009, to t10 = 2019
tt[] <- user()
dim(tt) <- user()
beta  <- interpolate(tt, beta_t,  "linear")
eta_l <- interpolate(tt, eta_l_t, "linear")
eta_h <- interpolate(tt, eta_h_t, "linear")
eta[1] <- eta_l
eta[2] <- eta_h

## Core equations for transitions between compartments:

deriv(U[, , ]) <- entrants[i, j, k] - n_UI[i, j, k] - exr * U[i, j, k] +
  n_AU[i, j, k] + n_TU[i, j, k]  + sum(wU[i, j, ]) -
  sum(n_vbe[i, j, ]) - sum(n_vod[i, j, ]) - sum(n_vos[i, j, ])

deriv(I[, , ]) <- n_UI[i, j, k] - (sigma + exr) * I[i, j, k] + sum(wI[i, j, ])

deriv(A[, , ]) <- (1 - (1 - ves[j]) * psi) * sigma * I[i, j, k]
  - n_AT[i, j, k] - n_AU[i, j, k] - exr * A[i, j, k] + sum(wA[i, j, ])

deriv(S[, , ]) <- (1 - ves[j]) * psi * sigma * I[i, j, k] - n_ST[i, j, k] -
  exr * S[i, j, k] + sum(wS[i, j, ])

deriv(T[, , ]) <- n_ST[i, j, k] + n_AT[i, j, k] - exr * T[i, j, k]
  - n_TU[i, j, k] + sum(wT[i, j, ])

## Update population size
N[, , ] <- U[i, j, k] + I[i, j, k] + A[i, j, k] + S[i, j, k] + T[i, j, k]

entrants[, , ] <- enr * q[i] * willing[j]

# calculate mixing matrix, probability of infection and force of infection
C[, , ] <- (1 - vei[j]) * (I[i, j, k] + A[i, j, k] + S[i, j, k])
prop_C[] <- sum(C[i, , ]) / sum(N[i, ,])
Np[]    <- sum(N[i, , ]) * p[i]

foi_LH[] <- prop_C[i] * Np[i] / sum(Np[])
lambda[] <- p[i] * beta * (epsilon * prop_C[i] + (1 - epsilon) * sum(foi_LH[]))

n_UI[, , ]     <- lambda[i] * (1 - vea[j]) * U[i, j, k]
n_AT[, , ]     <- eta[i] * A[i, j, k]
n_AU[, , ]     <- nu / (1 - ved[j]) * A[i, j, k]
n_ST[, , ]     <- mu * S[i, j, k]
n_TU[, , ]     <- rho * T[i, j, k]
screened[, ] <- eta[i] * U[i, j]

# vaccination
## time-varying switch
vax_switch <- interpolate(vax_t, vax_y, "constant")

## Number offered / accepting vaccine
## at screening
n_oos[, , ] <- vos[i, j, k] * screened[i, k] * vax_switch
n_vos[, , ] <- n_oos[i, j, k] * u[i, j, k]
## on diagnosis
n_ood[, , ] <- vod[i, j, k] * n_TU[i, k] * vax_switch
n_vod[, , ] <- n_ood[i, j, k] * u[i, j, k]
## on entry - no switch as background rate, adolescent uptake included in vbe
n_obe[, , ] <- vbe[i, j, k] * entrants[i, k]
n_vbe[, , ] <- n_obe[i, j, k] * u_vbe

# waning
wU[, , ] <- w[j, k] * U[i, k]
wI[, , ] <- w[j, k] * I[i, k]
wA[, , ] <- w[j, k] * A[i, k]
wS[, , ] <- w[j, k] * S[i, k]
wT[, , ] <- w[j, k] * T[i, k]

## outputs

deriv(cum_incid[, ])       <- n_UI[i, j]
deriv(cum_diag_a[, ])      <- n_AT[i, j]
deriv(cum_diag_s[, ])      <- n_ST[i, j]
deriv(cum_treated[, ])     <- n_TU[i, j]
deriv(cum_screened[, ])    <- screened[i, j]
deriv(cum_offered[, ])     <- n_oos[i, j, j] + n_ood[i, j, j] + n_obe[i, j, j]
deriv(cum_vaccinated[, ])  <- n_vos[i, j, j] + n_vod[i, j, j] + n_vbe[i, j, j]
deriv(cum_vbe[, ])         <- n_vbe[i, j, j]
deriv(cum_offered_vbe[, ]) <- n_obe[i, j, j]

# aggregated time series for fitting mcmc
output(tot_treated)  <- sum(cum_treated)
output(tot_attended) <- sum(cum_treated) + sum(cum_screened)

# output time-varying params for checking
output(beta) <- beta
output(eta) <- eta

## Set up compartments
## Initial states are all 0 as we will provide a state vbector
initial(U[, , ]) <- U0[i, j, k]
initial(I[, , ]) <- I0[i, j, k]
initial(A[, , ]) <- A0[i, j, k]
initial(S[, , ]) <- S0[i, j, k]
initial(T[, , ]) <- T0[i, j, k]

U0[, , ] <- user()
I0[, , ] <- user()
A0[, , ] <- user()
S0[, , ] <- user()
T0[, , ] <- user()

initial(cum_incid[, ])       <- 0
initial(cum_diag_a[, ])      <- 0
initial(cum_diag_s[, ])      <- 0
initial(cum_treated[, ])     <- 0
initial(cum_screened[, ])    <- 0
initial(cum_offered[, ])     <- 0
initial(cum_vaccinated[, ])  <- 0
initial(cum_vbe[, ])         <- 0
initial(cum_offered_vbe[, ]) <- 0

# set up dimensions of compartments
dim(U) <- c(n_group, n_vax, n_erlang)
dim(I) <- c(n_group, n_vax, n_erlang)
dim(A) <- c(n_group, n_vax, n_erlang)
dim(S) <- c(n_group, n_vax, n_erlang)
dim(T) <- c(n_group, n_vax, n_erlang)

dim(U0) <- c(n_group, n_vax, n_erlang)
dim(I0) <- c(n_group, n_vax, n_erlang)
dim(A0) <- c(n_group, n_vax, n_erlang)
dim(S0) <- c(n_group, n_vax, n_erlang)
dim(T0) <- c(n_group, n_vax, n_erlang)

dim(C)  <- c(n_group, n_vax, n_erlang)
dim(N)  <- c(n_group, n_vax, n_erlang)
dim(entrants) <- c(n_group, n_vax)
dim(Np)     <- n_group
dim(prop_C) <- n_group
dim(foi_LH) <- n_group
dim(lambda) <- n_group

dim(n_UI)     <- c(n_group, n_vax)
dim(n_AT)     <- c(n_group, n_vax)
dim(n_AU)     <- c(n_group, n_vax)
dim(n_ST)     <- c(n_group, n_vax)
dim(n_TU)     <- c(n_group, n_vax)
dim(screened) <- c(n_group, n_vax)

dim(cum_incid)       <- c(n_group, n_vax)
dim(cum_diag_a)      <- c(n_group, n_vax)
dim(cum_diag_s)      <- c(n_group, n_vax)
dim(cum_treated)     <- c(n_group, n_vax)
dim(cum_screened)    <- c(n_group, n_vax)
dim(cum_offered)     <- c(n_group, n_vax)
dim(cum_vaccinated)  <- c(n_group, n_vax)
dim(cum_vbe)         <- c(n_group, n_vax)
dim(cum_offered_vbe) <- c(n_group, n_vax)

## Parameters
p[]     <- user() # Partner change rate in group L/H
q[]     <- user() # Proportion in group L/H

willing[] <- user() # Proportion willing to be vaccinated
enr       <- user() # Rate of entry into population
exr       <- user() # Rate of exit from population
beta_t[]  <- user() # Time-varying rate of transmission
eta_l_t[] <- user() # Time-varying rate of screening in group L
eta_h_t[] <- user() # Time-varying rate of screening in group H
epsilon   <- user() # Level of assortative mixing
sigma     <- user() # Rate of developing symptoms
psi       <- user() # Proportion of infections that develop symptoms
nu        <- user() # Rate of natural recovery from asymptomatic infection
mu        <- user() # Rate of treatment seeking
rho       <- user() # Rate of recovery after treatment

## vaccination parameters
# vaccination routes
vbe[, , ] <- user() # vaccine map before entry
vos[, , ] <- user() # vaccine map on screening
vod[, , ] <- user() # vaccine map on diagnosis

# vaccine effects
vea[] <- user() # efficacy against acquisition
ved[] <- user() # efficacy against duration of infection
ves[] <- user() # efficacy against symptoms
vei[] <- user() # efficacy against infectiousness

u_vbe    <- user() # uptake of VbE
u[, , ]  <- user() # Uptake matrix for VoD/VoD
w[, ]    <- user() # Waning map
vax_t[]  <- user()
vax_y[]  <- user()

## par dimensions
dim(beta_t)  <- length(tt)
dim(eta_l_t) <- length(tt)
dim(eta_h_t) <- length(tt)

dim(p)    <- n_group
dim(q)    <- n_group
dim(eta)  <- n_group
dim(vea)  <- n_vax
dim(ved)  <- n_vax
dim(ves)  <- n_vax
dim(vei)  <- n_vax
dim(willing) <- n_vax
dim(vbe)   <- c(n_group, n_vax, n_vax)
dim(vod)   <- c(n_group, n_vax, n_vax)
dim(vos)   <- c(n_group, n_vax, n_vax)
dim(w)     <- c(n_vax, n_vax)
dim(u)     <- c(n_group, n_vax, n_vax)
dim(vax_t) <- user()
dim(vax_y) <- user()

dim(n_vbe) <- c(n_group, n_vax, n_vax)
dim(n_vos) <- c(n_group, n_vax, n_vax)
dim(n_vod) <- c(n_group, n_vax, n_vax)
dim(n_obe) <- c(n_group, n_vax, n_vax)
dim(n_oos) <- c(n_group, n_vax, n_vax)
dim(n_ood) <- c(n_group, n_vax, n_vax)
dim(wU)   <- c(n_group, n_vax, n_vax)
dim(wI)   <- c(n_group, n_vax, n_vax)
dim(wA)   <- c(n_group, n_vax, n_vax)
dim(wS)   <- c(n_group, n_vax, n_vax)
dim(wT)   <- c(n_group, n_vax, n_vax)

output(N)   <- N
output(lambda) <- lambda
