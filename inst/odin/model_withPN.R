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

deriv(U[, ]) <- entrants[i, j] - n_UI[i, j] - exr * U[i, j] +
  n_AU[i, j] + n_TU[i, j] - sum(n_diag_rec[i, j, ])  + sum(wU[i, j, ]) -
  sum(n_vbe[i, j, ]) - sum(n_vod[i, j, ]) - sum(n_vos[i, j, ]) -
  sum(n_vopn[i,j,])  + sum(wdU[i, j, ])

deriv(I[, ]) <- n_UI[i, j] - (sigma + exr) * I[i, j] + sum(wI[i, j, ]) +
  sum(wdI[i, j, ])

deriv(A[, ]) <- (1 - (1 - ves[j]) * psi) * sigma * I[i, j] - n_AT[i, j] -
  n_AU[i, j] - exr * A[i, j] + sum(wA[i, j, ]) + sum(wdA[i, j, ])

deriv(S[, ]) <- (1 - ves[j]) * psi * sigma * I[i, j] - n_ST[i, j] -
  exr * S[i, j] + sum(wS[i, j, ]) + sum(wdS[i, j, ])

deriv(T[, ]) <- n_ST[i, j] + n_AT[i, j] - exr * T[i, j] - n_TU[i, j] +
  sum(wT[i, j, ]) +  sum(wdT[i, j, ])

## Update population size
N[, ] <- U[i, j] + I[i, j] + A[i, j] + S[i, j] + T[i, j]

entrants[, ] <- enr * q[i] * willing[j]

# calculate mixing matrix, probability of infection and force of infection
C[, ] <- (1 - vei[j]) * (I[i, j] + A[i, j] + S[i, j])
prop_C[] <- sum(C[i, ]) / sum(N[i, ])
Np[]    <- sum(N[i, ]) * p[i]

foi_LH[] <- prop_C[i] * Np[i] / sum(Np[])
lambda[] <- p[i] * beta * (epsilon * prop_C[i] + (1 - epsilon) * sum(foi_LH[]))


n_UI[, ]     <- lambda[i] * (1 - vea[j]) * U[i, j]
n_AT[, ]     <- eta[i] * A[i, j]
n_AU[, ]     <- nu / (1 - ved[j]) * A[i, j]
n_ST[, ]     <- mu * S[i, j]
n_TU[, ]     <- rho * T[i, j]
screened[, ] <- eta[i] * U[i, j]



#Calculations required for partner notification

Cp[] <- sum(C[i, ])*p[i]

prop_Usubgroup[,] <- U[i,j] / sum(N[i, ]) 
prop_Asubgroup[,] <- A[i,j] / sum(N[i, ]) 
prop_ACsubgroup[,] <- A[i,j] / sum(C[i, ]) 


omega_N[ , ] = if (i == j) epsilon + (1-epsilon)*Np[j] else (1-epsilon)*Np[j] 
omega_C[ , ] = if (i == j) epsilon + (1-epsilon)*Cp[j] else (1-epsilon)*Cp[j] 

D_S <- 1 / sigma + 1 / mu
D_A[] <- 1 / sigma + 1 / (eta[i] + nu) 

lb_S <- 1/52
lb_A <- 3/12



#Probabilities only in terms of riskgroups
P_S_inf[ , ] = nu / (eta[j] + nu) * ( (1 / D_S) * exp( - ( (1 / D_S) + nu + eta[j] )*lb_S  ) - ( (1 / D_S) + nu + eta[j] ) * exp( - (1 / D_S)*lb_S) + eta[j] + nu ) / ((1 / D_S) + nu + eta[j] )
P_A_inf[ , ] = nu / (eta[j] + nu) * ( (1 / D_A[i]) * exp( - ( (1 / D_A[i]) + nu + eta[j] )*lb_A  ) - ( (1 / D_A[i]) + nu + eta[j] ) * exp( - (1 / D_A[i])*lb_A) + eta[j] + nu ) / ((1 / D_A[i]) + nu + eta[j] )

P_S_Acont[ , ] = nu / (eta[j] + nu) * (lb_S - (1 - exp ( - (eta[j] + nu)*lb_S))/(eta[j] + nu))
P_A_Acont[ , ] = nu / (eta[j] + nu) * (lb_A - (1 - exp ( - (eta[j] + nu)*lb_A))/(eta[j] + nu))


P_S_Ucont[ , ] = (1 - beta) * exp( - (1/D_S) * lb_S) * ( (1 - exp( - lambda[j] * lb_S)) /lambda[j]) +
  (1 - beta) * ( (1/D_S)*exp(-(lambda[j] + (1/D_S))*lb_S) - (lambda[j] + (1/D_S))*exp(- lambda[j]*lb_S) + lambda[j] )/(lambda[j] * (lambda[j] + (1/D_S))) +
  (lambda[j]*exp(-(lambda[j] + (1/D_S))*lb_S) - (lambda[j] + (1/D_S))*exp(- lambda[j]*lb_S) + (1/D_S) )/(lambda[j] * (lambda[j] + (1/D_S)))
P_A_Ucont[ , ] = (1 - beta) * exp( - (1/D_S) * lb_S) * ( (1 - exp( - lambda[j] * lb_S)) /lambda[j]) +
  (1 - beta) * ( (1/D_A[i])*exp(-(lambda[j] + (1/D_A[i]))*lb_A) - (lambda[j] + (1/D_A[i]))*exp(- lambda[j]*lb_A) + lambda[j] )/(lambda[j] * (lambda[j] + (1/D_A[i]))) +
  (lambda[j]*exp(-(lambda[j] + (1/D_A[i]))*lb_A) - (lambda[j] + (1/D_A[i]))*exp(- lambda[j]*lb_A) + (1/D_A[i]) )/(lambda[j] * (lambda[j] + (1/D_A[i])))


#Kvalues in terms of risk groups and vaccination groups
K_S_inf[,,,] = omega_C[i,k] * prop_ACsubgroup[k,l] * P_S_inf[i,k]
K_A_inf[,,,] = omega_C[i,k] * prop_ACsubgroup[k,l] * P_A_inf[i,k]
K_S_Acont[,,,] = p[i]*omega_N[i,k]*prop_Asubgroup[k,l]*P_S_Acont[i,k]
K_A_Acont[,,,] = p[i]*omega_N[i,k]*prop_Asubgroup[k,l]*P_A_Acont[i,k]
K_S_Ucont[,,,] = p[i]*omega_N[i,k]*prop_Usubgroup[k,l]*P_S_Ucont[i,k]
K_A_Ucont[,,,] = p[i]*omega_N[i,k]*prop_Usubgroup[k,l]*P_A_Ucont[i,k]

K_S[,,,] = kappa*(K_S_inf[i,j,k,l] + K_S_Acont[i,j,k,l] + K_S_Ucont[i,j,k,l])
K_A[,,,] = kappa*(K_A_inf[i,j,k,l] + K_A_Acont[i,j,k,l] + K_A_Ucont[i,j,k,l])

phi_group[,,,] = mu*S[i,j]*K_S[i,j,k,l] + eta[k]*A[i,j]*K_A[i,j,k,l]


phi[,] = sum(phi_group[,,i,j])


# mechanism to record number of times infected by moving diagnosed
# individuals into stratum with the relevant diagnosis history
n_diag_rec[, , ] <- diag_rec[i, j, k] * n_TU[i, k]


# vaccination
## time-varying switch
vax_switch <- interpolate(vax_t, vax_y, "constant")

## Number offered / accepting vaccine
## at screening
n_oos[, , ] <- vos[i, j, k] * screened[i, k] * vax_switch
n_vos[, , ] <- n_oos[i, j, k] * u_s[i, j, k]
## on diagnosis
n_ood[, , ] <- vod[i, j, k] * n_TU[i, k] * vax_switch
n_vod[, , ] <- n_ood[i, j, k] * u_d[i, j, k]
## on entry - no switch as background rate, adolescent uptake included in vbe
n_obe[, , ] <- vbe[i, j, k] * entrants[i, k]
n_vbe[, , ] <- n_obe[i, j, k] * u_vbe

# on partner notification
n_oopn[, , ] <- vopn[i,j,k] * phi[i, k] * vax_switch
n_vopn[, , ] <- n_oopn[i,j,k] * u_s[i, j, k]


# waning
wU[, , ] <- w[j, k] * U[i, k]
wI[, , ] <- w[j, k] * I[i, k]
wA[, , ] <- w[j, k] * A[i, k]
wS[, , ] <- w[j, k] * S[i, k]
wT[, , ] <- w[j, k] * T[i, k]


# waning diagnosis
wdU[, , ] <- wd[j, k] * U[i, k]
wdI[, , ] <- wd[j, k] * I[i, k]
wdA[, , ] <- wd[j, k] * A[i, k]
wdS[, , ] <- wd[j, k] * S[i, k]
wdT[, , ] <- wd[j, k] * T[i, k]




## outputs

deriv(cum_incid[, ])       <- n_UI[i, j]
deriv(cum_diag_a[, ])      <- n_AT[i, j]
deriv(cum_diag_s[, ])      <- n_ST[i, j]
deriv(cum_treated[, ])     <- n_TU[i, j]
deriv(cum_screened[, ])    <- screened[i, j]
deriv(cum_offered[, ])     <- n_oos[i, j, j] + n_ood[i, j, j] + n_obe[i, j, j]
deriv(cum_vaccinated[, ])  <- n_vos[i, j, j] + n_vod[i, j, j] + n_vbe[i, j, j]
deriv(cum_vaccinated_screen[, ])  <- n_vos[i, j, j]

deriv(cum_vbe[, ])         <- n_vbe[i, j, j]
deriv(cum_offered_vbe[, ]) <- n_obe[i, j, j]
deriv(cum_entrants[, ]) <- entrants[i, j]


# aggregated time series for fitting mcmc
output(tot_treated)  <- sum(cum_treated)
output(tot_attended) <- sum(cum_treated) + sum(cum_screened)

# output time-varying params for checking
output(beta) <- beta
output(eta) <- eta

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

initial(cum_incid[, ])       <- 0
initial(cum_diag_a[, ])      <- 0
initial(cum_diag_s[, ])      <- 0
initial(cum_treated[, ])     <- 0
initial(cum_screened[, ])    <- 0
initial(cum_offered[, ])     <- 0
initial(cum_vaccinated[, ])  <- 0
initial(cum_vaccinated_screen[, ])  <- 0

initial(cum_vbe[, ])         <- 0
initial(cum_offered_vbe[, ]) <- 0

initial(cum_entrants[, ]) <- 0



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

dim(C)  <- c(n_group, n_vax)
dim(N)  <- c(n_group, n_vax)
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


dim(Cp)     <- n_group
dim(omega_N) <- c(n_group, n_group)
dim(omega_C) <- c(n_group, n_group)
dim(prop_Usubgroup) <- c(n_group, n_vax)
dim(prop_Asubgroup) <- c(n_group, n_vax)
dim(prop_ACsubgroup) <- c(n_group, n_vax)

dim(D_A) <- n_group
dim(P_S_inf) <- c(n_group, n_group)
dim(P_A_inf) <- c(n_group, n_group)
dim(P_S_Acont) <- c(n_group, n_group)
dim(P_A_Acont) <- c(n_group, n_group)
dim(P_S_Ucont) <- c(n_group, n_group)
dim(P_A_Ucont) <- c(n_group, n_group)

dim(K_S_inf) <- c(n_group, n_vax, n_group, n_vax)
dim(K_A_inf) <- c(n_group, n_vax, n_group, n_vax)
dim(K_S_Acont) <- c(n_group, n_vax, n_group, n_vax)
dim(K_A_Acont) <- c(n_group, n_vax, n_group, n_vax)
dim(K_S_Ucont) <- c(n_group, n_vax, n_group, n_vax)
dim(K_A_Ucont) <- c(n_group, n_vax, n_group, n_vax)

dim(K_S) <- c(n_group, n_vax, n_group, n_vax)
dim(K_A) <- c(n_group, n_vax, n_group, n_vax)

dim(phi_group) <- c(n_group, n_vax, n_group, n_vax)
dim(phi) <- c(n_group, n_vax)

dim(cum_incid)       <- c(n_group, n_vax)
dim(cum_diag_a)      <- c(n_group, n_vax)
dim(cum_diag_s)      <- c(n_group, n_vax)
dim(cum_treated)     <- c(n_group, n_vax)
dim(cum_screened)    <- c(n_group, n_vax)
dim(cum_offered)     <- c(n_group, n_vax)
dim(cum_vaccinated)  <- c(n_group, n_vax)
dim(cum_vaccinated_screen)  <- c(n_group, n_vax)


dim(cum_vbe)         <- c(n_group, n_vax)
dim(cum_offered_vbe) <- c(n_group, n_vax)

dim(cum_entrants) <- c(n_group, n_vax)


dim(n_diag_rec) <- c(n_group, n_vax, n_vax)
dim(diag_rec)   <- c(n_group, n_vax, n_vax)


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


kappa     <- user() # proportion of eligible contacts for notification actually notified

## vaccination parameters
# vaccination routes
vbe[, , ] <- user() # vaccine map before entry
vos[, , ] <- user() # vaccine map on screening
vod[, , ] <- user() # vaccine map on diagnosis
vopn[, , ] <- user() # vaccine map on partner notification

# vaccine effects
vea[] <- user() # efficacy against acquisition
ved[] <- user() # efficacy against duration of infection
ves[] <- user() # efficacy against symptoms
vei[] <- user() # efficacy against infectiousness

u_vbe    <- user() # uptake of VbE
u_d[, , ]  <- user() # Uptake matrix for diagnosis
u_s[, , ]  <- user() # Uptake matrix for screening

w[, ]    <- user() # Waning map
vax_t[]  <- user()
vax_y[]  <- user()

wd[, ]   <- user() # Waning of diagnosis map

diag_rec[, , ] <- user()  ## recording diagnosis history mapping




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
dim(vopn)   <- c(n_group, n_vax, n_vax)


dim(w)     <- c(n_vax, n_vax)


dim(u_d)     <- c(n_group, n_vax, n_vax)
dim(u_s)     <- c(n_group, n_vax, n_vax)


dim(vax_t) <- user()
dim(vax_y) <- user()

dim(n_vbe) <- c(n_group, n_vax, n_vax)
dim(n_vos) <- c(n_group, n_vax, n_vax)
dim(n_vod) <- c(n_group, n_vax, n_vax)

dim(n_vopn) <- c(n_group, n_vax, n_vax)



dim(n_obe) <- c(n_group, n_vax, n_vax)
dim(n_oos) <- c(n_group, n_vax, n_vax)
dim(n_ood) <- c(n_group, n_vax, n_vax)

dim(n_oopn) <- c(n_group, n_vax, n_vax)



dim(wU)   <- c(n_group, n_vax, n_vax)
dim(wI)   <- c(n_group, n_vax, n_vax)
dim(wA)   <- c(n_group, n_vax, n_vax)
dim(wS)   <- c(n_group, n_vax, n_vax)
dim(wT)   <- c(n_group, n_vax, n_vax)



dim(wd) <- c(n_vax, n_vax)
dim(wdU)   <- c(n_group, n_vax, n_vax)
dim(wdI)   <- c(n_group, n_vax, n_vax)
dim(wdA)   <- c(n_group, n_vax, n_vax)
dim(wdS)   <- c(n_group, n_vax, n_vax)
dim(wdT)   <- c(n_group, n_vax, n_vax)

output(N)   <- N
output(lambda) <- lambda
