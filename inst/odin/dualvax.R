### Key to indices
## i for the activity group (low / high)
##  the compartments will then be replicated five times
### : unvaccinated
### B: Bexsero
### X: Vaccine X
### D: Dual
### W: Waned

# Group is either 1: Low, or 2: High activity
n_group <- 2 

## Core equations for transitions between compartments:

deriv(U[]) <- enr[i] - incid[i]  - exr*U[i] + nu*A[i] + rho*T[i]
deriv(I[]) <- incid[i]           - (sigma+exr)*I[i] 
deriv(A[]) <- (1-psi)*sigma*I[i] - (eta+nu+exr)*A[i]
deriv(S[]) <- psi*sigma*I[i]     - (mu+exr)*S[i]
deriv(T[]) <- mu*S[i] + eta*A[i] - (rho+exr)*T[i]

deriv(cum_incid[]) <- incid[i]


## Update population size
C[] <- I[i] + A[i] + S[i] 
N[] <- U[i] + C[i] + T[i]

# calculate mixing matrix, probability of infection and force of infection
Np[] <- N[i] * p[i]
m[, ] <- epsilon * (i==j) + (1-epsilon) * Np[j] / sum(Np)
foi[, ] <- beta * p[i] * m[i, j] * C[j] / N[j]
incid[] <- sum(foi[i, ]) * U[i]

## Set up compartments
## Initial states are all 0 as we will provide a state vector
initial(U[]) <- U0[i]
initial(I[]) <- I0[i]
initial(A[]) <- A0[i]
initial(S[]) <- S0[i]
initial(T[]) <- T0[i]

initial(cum_incid[]) <- 0

U0[] <- user()
I0[] <- user()
A0[] <- user()
S0[] <- user()
T0[] <- user()

# set up dimensions of compartments
dim(U) <- n_group
dim(I) <- n_group
dim(A) <- n_group
dim(S) <- n_group
dim(T) <- n_group

dim(U0) <- n_group
dim(I0) <- n_group
dim(A0) <- n_group
dim(S0) <- n_group
dim(T0) <- n_group

dim(C) <- n_group
dim(N) <- n_group
dim(Np) <- n_group

dim(m) <- c(n_group, n_group)
dim(incid) <- n_group
dim(foi) <- c(n_group, n_group)

dim(cum_incid) <- n_group

## Parameters
enr[] <- user()
exr <- user()
beta <- user()
p[] <- user()
epsilon <- user()
sigma <- user()
psi <- user()
eta <- user()
nu <- user()
mu <- user()
rho <- user()

dim(enr) <- n_group
dim(p) <- n_group

output(N) <- N
output(foi) <- foi

