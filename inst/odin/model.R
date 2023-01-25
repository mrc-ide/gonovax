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
n_erlang <- 1                    #have set to 1 for now!> change this to user(1)

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
  n_AU[i, j, k] + n_TU[i, j, k]  + sum(wU[i, j, , k]) -
  sum(n_vbe[i, j, , k]) - sum(n_vod[i, j, , k]) - sum(n_vos[i, j, , k])

deriv(I[, , ]) <- n_UI[i, j, k] - (sigma + exr) * I[i, j, k] +
  sum(wI[i, j, k, ])

deriv(A[, , ]) <- (1 - (1 - ves[j]) * psi) * sigma * I[i, j, k] -
  n_AT[i, j, k] - n_AU[i, j, k] - exr * A[i, j, k] + sum(wA[i, j, , k])

deriv(S[, , ]) <- (1 - ves[j]) * psi * sigma * I[i, j, k] - n_ST[i, j, k] -
  exr * S[i, j, k] + sum(wS[i, j, , k])

deriv(T[, , ]) <- n_ST[i, j, k] + n_AT[i, j, k] - exr * T[i, j, k] -
  n_TU[i, j, k] + sum(wT[i, j, , k])

## Update popukation size
N[, , ] <- U[i, j, k] + I[i, j, k] + A[i, j, k] + S[i, j, k] + T[i, j, k]

entrants[, , ] <- enr * q[i] * willing[j]

# calculate mixing matrix, probability of infection and force of infection
C[, , ] <- (1 - vei[j]) * (I[i, j, k] + A[i, j, k] + S[i, j, k])
prop_C[] <- sum(C[i, , ]) / sum(N[i, , ])
Np[]    <- sum(N[i, , ]) * p[i]

foi_LH[] <- prop_C[i] * Np[i] / sum(Np[])
lambda[] <- p[i] * beta * (epsilon * prop_C[i] + (1 - epsilon) * sum(foi_LH[]))

n_UI[, , ]     <- lambda[i] * (1 - vea[j]) * U[i, j, k]
n_AT[, , ]     <- eta[i] * A[i, j, k]
n_AU[, , ]     <- nu / (1 - ved[j]) * A[i, j, k]
n_ST[, , ]     <- mu * S[i, j, k]
n_TU[, , ]     <- rho * T[i, j, k]
screened[, , ] <- eta[i] * U[i, j, k]

# vaccination
## time-varying switch
vax_switch <- interpolate(vax_t, vax_y, "constant")

## Number offered / accepting vaccine
## at screening
n_oos[, , , ] <- vos[i, j, k, l] * screened[i, k, l] * vax_switch
n_vos[, , , ] <- n_oos[i, j, k, l] * u[i, j, k, l]
## on diagnosis
n_ood[, , , ] <- vod[i, j, k, l] * n_TU[i, k, l] * vax_switch
n_vod[, , , ] <- n_ood[i, j, k, l] * u[i, j, k, l]
## on entry - no switch as background rate, adolescent uptake included in vbe
n_obe[, , , ] <- vbe[i, j, k, l] * entrants[i, k, l]
n_vbe[, , , ] <- n_obe[i, j, k, l] * u_vbe

# waning
wU[, , , ] <- w[j, k, l] * U[i, k, l]
wI[, , , ] <- w[j, k, l] * I[i, k, l]
wA[, , , ] <- w[j, k, l] * A[i, k, l]           #come back here
wS[, , , ] <- w[j, k, l] * S[i, k, l]
wT[, , , ] <- w[j, k, l] * T[i, k, l]

## outputs

deriv(cum_incid[, , ])       <- n_UI[i, j, k]
deriv(cum_diag_a[, , ])      <- n_AT[i, j, k]
deriv(cum_diag_s[, , ])      <- n_ST[i, j, k]
deriv(cum_treated[, , ])     <- n_TU[i, j, k]
deriv(cum_screened[, , ])    <- screened[i, j, k]
deriv(cum_offered[, , ])     <- n_oos[i, j, j, k] + n_ood[i, j, j, k] +
  n_obe[i, j, j, k]
deriv(cum_vaccinated[, , ])  <- n_vos[i, j, j, k] + n_vod[i, j, j, k] +
  n_vbe[i, j, j, k]
deriv(cum_vbe[, , ])         <- n_vbe[i, j, j, k]
deriv(cum_offered_vbe[, , ]) <- n_obe[i, j, j, k]

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

initial(cum_incid[, , ])       <- 0
initial(cum_diag_a[, , ])      <- 0
initial(cum_diag_s[, , ])      <- 0
initial(cum_treated[, , ])     <- 0
initial(cum_screened[, , ])    <- 0
initial(cum_offered[, , ])     <- 0
initial(cum_vaccinated[, , ])  <- 0
initial(cum_vbe[, , ])         <- 0
initial(cum_offered_vbe[, , ]) <- 0

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

dim(C)  <- c(n_group, n_vax, n_erlang)    #unsure about this one
dim(N)  <- c(n_group, n_vax, n_erlang)
dim(entrants) <- c(n_group, n_vax, n_erlang)
dim(Np)     <- n_group
dim(prop_C) <- n_group
dim(foi_LH) <- n_group
dim(lambda) <- n_group

dim(n_UI)     <- c(n_group, n_vax, n_erlang)
dim(n_AT)     <- c(n_group, n_vax, n_erlang)
dim(n_AU)     <- c(n_group, n_vax, n_erlang)
dim(n_ST)     <- c(n_group, n_vax, n_erlang)
dim(n_TU)     <- c(n_group, n_vax, n_erlang)
dim(screened) <- c(n_group, n_vax, n_erlang)

dim(cum_incid)       <- c(n_group, n_vax, n_erlang)
dim(cum_diag_a)      <- c(n_group, n_vax, n_erlang)
dim(cum_diag_s)      <- c(n_group, n_vax, n_erlang)
dim(cum_treated)     <- c(n_group, n_vax, n_erlang)
dim(cum_screened)    <- c(n_group, n_vax, n_erlang)
dim(cum_offered)     <- c(n_group, n_vax, n_erlang)
dim(cum_vaccinated)  <- c(n_group, n_vax, n_erlang)
dim(cum_vbe)         <- c(n_group, n_vax, n_erlang)
dim(cum_offered_vbe) <- c(n_group, n_vax, n_erlang)

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
vbe[, , , ] <- user() # vaccine map before entry
vos[, , , ] <- user() # vaccine map on screening
vod[, , , ] <- user() # vaccine map on diagnosis

# vaccine effects
vea[] <- user() # efficacy against acquisition
ved[] <- user() # efficacy against duration of infection
ves[] <- user() # efficacy against symptoms
vei[] <- user() # efficacy against infectiousness

u_vbe    <- user() # uptake of VbE
u[, , , ]  <- user() # Uptake matrix for VoD/VoD
w[, , ]    <- user() # Waning map
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
dim(vbe)   <- c(n_group, n_vax, n_vax,  n_erlang)
dim(vod)   <- c(n_group, n_vax, n_vax,  n_erlang)
dim(vos)   <- c(n_group, n_vax, n_vax,  n_erlang)
dim(w)     <- c(n_vax, n_vax, n_erlang)
dim(u)     <- c(n_group, n_vax, n_vax, n_erlang) ## This may not need an erlang dimension
dim(vax_t) <- user()
dim(vax_y) <- user()

dim(n_vbe) <- c(n_group, n_vax, n_vax, n_erlang)
dim(n_vos) <- c(n_group, n_vax, n_vax, n_erlang)
dim(n_vod) <- c(n_group, n_vax, n_vax, n_erlang)
dim(n_obe) <- c(n_group, n_vax, n_vax, n_erlang)
dim(n_oos) <- c(n_group, n_vax, n_vax, n_erlang)
dim(n_ood) <- c(n_group, n_vax, n_vax, n_erlang)
dim(wU)   <- c(n_group, n_vax, n_vax, n_erlang)
dim(wI)   <- c(n_group, n_vax, n_vax, n_erlang)
dim(wA)   <- c(n_group, n_vax, n_vax, n_erlang)
dim(wS)   <- c(n_group, n_vax, n_vax, n_erlang)
dim(wT)   <- c(n_group, n_vax, n_vax, n_erlang)

output(N)   <- N
output(lambda) <- lambda
