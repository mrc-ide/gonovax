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
n_par <- length(beta)
n_group <- 2
n_vax <- 1

## Core equations for transitions between compartments:

deriv(U[, ]) <- enr * q[j] * (1 - ve) - n_UI[i, j] - exr * U[i, j] +
  nu[i] * A[i, j] + (1 - vd) * rho[i] * T[i, j]
deriv(I[, ]) <- n_UI[i, j] - (sigma[i] + exr) * I[i, j]
deriv(A[, ]) <- (1 - psi[i]) * sigma[i] * I[i, j] -
  (eta[i] + nu[i] + exr) * A[i, j, k]
deriv(S[, ]) <- psi[i] * sigma[i] * I[i, j] - (mu[i] + exr) * S[i, j]
deriv(T[, ]) <- mu[i] * S[i, j] + eta[i] * A[i, j] - (rho[i] + exr) * T[i, j]

deriv(cum_incid[, , ]) <- n_UI[i, j]


## Update population size
C[, ] <- I[i, j] + A[i, j] + S[i, j] 
N[, ] <- U[i, j] + C[i, j] + T[i, j]

# calculate mixing matrix, probability of infection and force of infection
Np[, ] <- sum(N[i, j]) * p[j]
m[, , ] <- epsilon[i] * (j == k) + (1 - epsilon[i]) * Np[i] / sum(Np[i, ])
foi[, , ] <- beta[i] * p[j] * m[i, j, k] * sum(C[i, k]) / sum(N[i, k])
n_UI[, ] <- sum(foi[i, j, ]) * U[i, j]
n_BI[, ] <- sum(foi[i, j, ]) * eff * BU[i, j]


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

initial(cum_incid[, ]) <- 0

# set up dimensions of compartments
dim(U) <- c(n_par, n_group)
dim(I) <- c(n_par, n_group)
dim(A) <- c(n_par, n_group)
dim(S) <- c(n_par, n_group)
dim(T) <- c(n_par, n_group)

dim(U0) <- c(n_par, n_group)
dim(I0) <- c(n_par, n_group)
dim(A0) <- c(n_par, n_group)
dim(S0) <- c(n_par, n_group)
dim(T0) <- c(n_par, n_group)

dim(C)  <- c(n_par, n_group)
dim(N)  <- c(n_par, n_group)
dim(Np) <- c(n_par, n_group)

dim(m)    <- c(n_par, n_group, n_group)
dim(foi)  <- c(n_par, n_group, n_group)
dim(n_UI) <- c(n_par, n_group)
dim(cum_incid) <- c(n_par, n_group)

# Parameters
q[] <- user()
enr <- user()
exr <- user()
beta[] <- user()
p[] <- user()
epsilon[] <- user()
sigma[] <- user()
psi[] <- user()
eta[] <- user()
nu[] <- user()
mu[] <- user()
rho[] <- user()

## vaccination pars
ve <- user()
vd <- user()
va <- user()
eff <- user()
w <- user()


## par dimensions
dim(q) <- n_group
dim(beta) <- user()
dim(epsilon) <- n_par
dim(sigma) <- n_par
dim(psi) <- n_par
dim(eta) <- n_par
dim(nu) <- n_par
dim(mu) <- n_par
dim(rho) <- n_par

dim(p)   <- n_group

output(N) <- N
output(foi) <- foi
