### Key to indices
## i = parameter set
## j for the activity group (low / high)
## k the compartments will then be replicated five times
### : unvaccinated
### B: Bexsero
### X: Vaccine X
### D: Dual
### W: Waned

# Group is either 1: Low, or 2: High activity
n_par   <- length(beta)
n_group <- 2
n_vax   <- user(1)


## Core equations for transitions between compartments:

deriv(U[, , ]) <- enr * q[j] * ve[k] - incid[i, j, k] - exr * U[i, j, k] +
  nu[i] * A[i, j, k] + treated[i, j, k] -
  sum(n_vd[i, j, k, ]) - sum(n_vs[i, j, k, ]) + sum(wU[i, j, k, ])
deriv(I[, , ]) <- incid[i, j, k] - (sigma[i] + exr) * I[i, j, k] +
  sum(wI[i, j, k, ])
deriv(A[, , ]) <- (1 - psi[i]) * sigma[i] * I[i, j, k] -
  (eta[i] + nu[i] + exr) * A[i, j, k] +
  sum(wA[i, j, k, ])
deriv(S[, , ]) <- psi[i] * sigma[i] * I[i, j, k] - (mu[i] + exr) * S[i, j, k] +
  sum(wS[i, j, k, ])
deriv(T[, , ]) <- mu[i] * S[i, j, k] + eta[i] * A[i, j, k] -
  exr * T[i, j, k] - treated[i, j, k] + sum(wT[i, j, k, ])

## Update population size
C[, , ] <- I[i, j, k] + A[i, j, k] + S[i, j, k]
N[, , ] <- U[i, j, k] + C[i, j, k] + T[i, j, k]

# calculate mixing matrix, probability of infection and force of infection
Np[, ]    <- sum(N[i, j, ]) * p[j]
m[, , ]   <- epsilon[i] * (j == k) + (1 - epsilon[i]) * Np[i, k] / sum(Np[i, ])
foi[, , ] <- beta[i] * p[j] * m[i, j, k] * sum(C[i, k, ]) / sum(N[i, k, ])

incid[, , ]    <- sum(foi[i, j, ]) * (1 - eff[k]) * U[i, j, k]
treated[, , ]  <- rho[i] * T[i, j, k]
screened[, , ] <- eta[i] * U[i, j, k]
# vaccination
## input vax_map (v)
## (1, 0,
## -1, 0)
## at screening
n_vs[, , , ] <- vs[k, l] * screened[i, j, l]
## on diagnosis
n_vd[, , , ] <- vd[k, l] * treated[i, j, l]
## on entry
n_ve[, , 2:n_vax] <- enr * q[j] * ve[k]

# waning
## input waning_map
## (0,  w,
##  0, -w)
wU[, , , ] <- w[k, l] * U[i, j, l]
wI[, , , ] <- w[k, l] * I[i, j, l]
wA[, , , ] <- w[k, l] * A[i, j, l]
wS[, , , ] <- w[k, l] * S[i, j, l]
wT[, , , ] <- w[k, l] * T[i, j, l]

## outputs

deriv(cum_incid[, , ])      <- incid[i, j, k]
deriv(cum_treated[, , ])    <- treated[i, j, k]
deriv(cum_screened[, , ])   <- screened[i, j, k]
deriv(cum_vaccinated[, , ]) <-
  n_vs[i, j, k, k] + n_vd[i, j, k, k] + n_ve[i, j, k]
deriv(cum_screened[, , ])   <- screened[i, j, k]

## Set up compartments
## Initial states are all 0 as we will provide a state vector
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

initial(cum_incid[, , ])      <- 0
initial(cum_treated[, , ])    <- 0
initial(cum_screened[, , ])   <- 0
initial(cum_vaccinated[, , ]) <- 0

# set up dimensions of compartments
dim(U) <- c(n_par, n_group, n_vax)
dim(I) <- c(n_par, n_group, n_vax)
dim(A) <- c(n_par, n_group, n_vax)
dim(S) <- c(n_par, n_group, n_vax)
dim(T) <- c(n_par, n_group, n_vax)

dim(U0) <- c(n_par, n_group, n_vax)
dim(I0) <- c(n_par, n_group, n_vax)
dim(A0) <- c(n_par, n_group, n_vax)
dim(S0) <- c(n_par, n_group, n_vax)
dim(T0) <- c(n_par, n_group, n_vax)

dim(C)  <- c(n_par, n_group, n_vax)
dim(N)  <- c(n_par, n_group, n_vax)
dim(Np) <- c(n_par, n_group)

dim(m)   <- c(n_par, n_group, n_group)
dim(foi) <- c(n_par, n_group, n_group)

dim(incid)    <- c(n_par, n_group, n_vax)
dim(treated)  <- c(n_par, n_group, n_vax)
dim(screened) <- c(n_par, n_group, n_vax)

dim(cum_incid)      <- c(n_par, n_group, n_vax)
dim(cum_treated)    <- c(n_par, n_group, n_vax)
dim(cum_screened)   <- c(n_par, n_group, n_vax)
dim(cum_vaccinated) <- c(n_par, n_group, n_vax)

## Parameters
q[]       <- user()
enr       <- user()
exr       <- user()
beta[]    <- user()
p[]       <- user()
epsilon[] <- user()
sigma[]   <- user()
psi[]     <- user()
eta[]     <- user()
nu[]      <- user()
mu[]      <- user()
rho[]     <- user()

## vaccination pars
ve[]      <- user()
vs[, ]    <- user()
vd[, ]    <- user()
eff[]     <- user()
w[, ]     <- user()


## par dimensions
dim(q)       <- n_group
dim(beta)    <- user()
dim(epsilon) <- n_par
dim(sigma)   <- n_par
dim(psi)     <- n_par
dim(eta)     <- n_par
dim(nu)      <- n_par
dim(mu)      <- n_par
dim(rho)     <- n_par

dim(ve)   <- n_vax
dim(eff)  <- n_vax
dim(vd)   <- c(n_vax, n_vax)
dim(vs)   <- c(n_vax, n_vax)
dim(w)    <- c(n_vax, n_vax)
dim(n_ve) <- c(n_par, n_group, n_vax)
dim(n_vs) <- c(n_par, n_group, n_vax, n_vax)
dim(n_vd) <- c(n_par, n_group, n_vax, n_vax)
dim(wU)   <- c(n_par, n_group, n_vax, n_vax)
dim(wI)   <- c(n_par, n_group, n_vax, n_vax)
dim(wA)   <- c(n_par, n_group, n_vax, n_vax)
dim(wS)   <- c(n_par, n_group, n_vax, n_vax)
dim(wT)   <- c(n_par, n_group, n_vax, n_vax)

dim(p)   <- n_group

output(N)   <- N
output(foi) <- foi
