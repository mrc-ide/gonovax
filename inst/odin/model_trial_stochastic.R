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

# individual probabilities of transitioning between infection states
r_AT[, ] <- eta
r_AU[, ] <- nu / (1 - ved[j])

# probability of individuals leaving compartments
p_U_ext[, ] <- 1 - exp(-(lambda * (1 - vea[j]) - D[j]) * dt)
p_I_ext[, ] <- 1 - exp(-(sigma - D[j]) * dt)
p_S_ext[, ] <- 1 - exp(-(mu - D[j]) * dt)
p_A_ext[, ] <- 1 - exp(-(r_AT[i, j] + r_AU[i, j] - D[j]) * dt)
p_T_ext[, ] <- 1 - exp(-(rho - D[j]) * dt)

# draws from binomial distribution for numbers exiting compartments
n_U_ext[, ] <- rbinom(U[i, j], p_U_ext[i, j])
n_I_ext[, ] <- rbinom(I[i, j], p_I_ext[i, j])
n_S_ext[, ] <- rbinom(S[i, j], p_S_ext[i, j])
n_A_ext[, ] <- rbinom(A[i, j], p_A_ext[i, j])
n_T_ext[, ] <- rbinom(T[i, j], p_T_ext[i, j])

#Relative probabilities of infection transition vs vaccine transition
Rel_U[, ] <- if ((lambda * (1 - vea[j]) - D[j]) == 0) 0 else
  lambda * (1 - vea[j]) / (lambda * (1 - vea[j]) - D[j])

Rel_I[, ] <- if (sigma - D[j] == 0) 0 else sigma / (sigma - D[j])
Rel_S[, ] <- if (mu - D[j] == 0) 0 else mu / (mu - D[j])
Rel_A[, ] <- if (r_AT[i, j] + r_AU[i, j] - D[j] == 0) 0 else
  (r_AT[i, j] + r_AU[i, j]) / (r_AT[i, j] + r_AU[i, j] - D[j])

Rel_T[, ] <- if (rho - D[j] == 0) 0 else  rho / (rho - D[j])

#draw numbers changing between specific compartments
n_UI[, ] <- rbinom(n_U_ext[i, j], Rel_U[i, j])
n_Uw[, ] <- n_U_ext[i, j] - n_UI[i, j]

n_IAS[, ] <- rbinom(n_I_ext[i, j], Rel_I[i, j])
n_Iw[, ] <- n_I_ext[i, j] - n_IAS[i, j]
n_IA[, ] <- rbinom(n_IAS[i, j], 1 - (1 - ves[j]) * psi)
n_IS[, ] <- n_IAS[i, j] - n_IA[i, j]

n_ST[, ] <- rbinom(n_S_ext[i, j], Rel_S[i, j])
n_Sw[, ] <- n_S_ext[i, j] - n_ST[i, j]

n_AUT[, ] <- rbinom(n_A_ext[i, j], Rel_A[i, j])
n_Aw[, ] <- n_A_ext[i, j] - n_AUT[i, j]
n_AT[, ] <- rbinom(n_AUT[i, j], r_AT[i, j] / (r_AT[i, j] + r_AU[i, j]))
n_AU[, ] <- n_AUT[i, j] - n_AT[i, j]

n_TU[, ] <- rbinom(n_T_ext[i, j], Rel_T[i, j])
n_Tw[, ] <- n_T_ext[i, j] - n_TU[i, j]


# mechanism to record number of times infected by moving diagnosed
# individuals into stratum with the relevant diagnosis history

n_diag_rec[, , ] <- diag_rec[i, j, k] * n_TU[i, k]

#waning
wU[, , ] <- w[j, k] * n_Uw[i, k]
wI[, , ] <- w[j, k] * n_Iw[i, k]
wA[, , ] <- w[j, k] * n_Aw[i, k]
wS[, , ] <- w[j, k] * n_Sw[i, k]
wT[, , ] <- w[j, k] * n_Tw[i, k]

## Core equations for transitions between compartments:
update(U[, ]) <- U[i, j] - n_UI[i, j] +
  n_AU[i, j] + n_TU[i, j]  + sum(wU[i, j, ]) - sum(n_diag_rec[i, j, ])

update(I[, ]) <- I[i, j] + n_UI[i, j] - n_IAS[i, j] + sum(wI[i, j, ])

update(A[, ]) <- A[i, j] + n_IA[i, j] - n_AUT[i, j] + sum(wA[i, j, ])

update(S[, ]) <- S[i, j] + n_IS[i, j] - n_ST[i, j]  + sum(wS[i, j, ])

update(T[, ]) <- T[i, j] + n_ST[i, j] + n_AT[i, j] - n_TU[i, j] +
  sum(wT[i, j, ])

## Update population size
N[, ] <- U[i, j] + I[i, j] + A[i, j] + S[i, j] + T[i, j]

screened[, ] <- rbinom(U[i, j], 1 - exp(-eta *  dt))

pye_trial[, ] <- U[i, j] + I[i, j] + A[i, j] + S[i, j]

## outputs
update(cum_incid[, ])         <- cum_incid[i, j] + n_UI[i, j]
update(cum_diag_a[, ])        <- cum_diag_a[i, j] + n_AT[i, j]
update(cum_diag_s[, ])        <- cum_diag_s[i, j] + n_ST[i, j]
update(cum_treated[, ])       <- cum_treated[i, j] + n_TU[i, j]
update(cum_screened[, ])      <- cum_screened[i, j] + screened[i, j]
update(cum_pye_trial_pov[, ]) <- cum_pye_trial_pov[i, j] + pye_trial[i, j]
update(cum_pye_true[, ])      <- cum_pye_true[i, j] + U[i, j]

# aggregated time series for fitting mcmc
output(tot_treated) <- sum(cum_treated)
output(tot_attended) <- sum(cum_treated) + sum(cum_screened)

## Set up compartments
## Initial states are all 0 as we will provide a state vector
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
initial(cum_pye_trial_pov[, ]) <- 0
initial(cum_pye_true[, ]) <- 0

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

dim(pye_trial) <- c(n_group, n_vax)

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

# Probability of individuals leaving compartments
dim(p_U_ext) <- c(n_group, n_vax)
dim(p_I_ext) <- c(n_group, n_vax)
dim(p_S_ext) <- c(n_group, n_vax)
dim(p_A_ext) <- c(n_group, n_vax)
dim(p_T_ext) <- c(n_group, n_vax)

#Number of individuals leaving compartments
dim(n_U_ext) <- c(n_group, n_vax)
dim(n_I_ext) <- c(n_group, n_vax)
dim(n_S_ext) <- c(n_group, n_vax)
dim(n_A_ext) <- c(n_group, n_vax)
dim(n_T_ext) <- c(n_group, n_vax)

#Relative probabilities of infection transition
# vs vaccine  transition
dim(Rel_U) <- c(n_group, n_vax)
dim(Rel_I) <- c(n_group, n_vax)
dim(Rel_A) <- c(n_group, n_vax)
dim(Rel_S) <- c(n_group, n_vax)
dim(Rel_T) <- c(n_group, n_vax)

dim(r_AT) <- c(n_group, n_vax)
dim(r_AU) <- c(n_group, n_vax)

dim(cum_incid)         <- c(n_group, n_vax)
dim(cum_diag_a)        <- c(n_group, n_vax)
dim(cum_diag_s)        <- c(n_group, n_vax)
dim(cum_treated)       <- c(n_group, n_vax)
dim(cum_screened)      <- c(n_group, n_vax)
dim(cum_pye_trial_pov) <- c(n_group, n_vax)
dim(cum_pye_true)      <- c(n_group, n_vax)

dim(n_diag_rec) <- c(n_group, n_vax, n_vax)
dim(diag_rec)   <- c(n_group, n_vax, n_vax)

## Parameters
eta       <- user()
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

# mapping
w[, ]    <- user()
D[] <- user()

diag_rec[, , ] <- user()

## par dimensions
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

output(N) <- N
