### Key to indices
## i for the parameter vector
## j for the activity group (low / high)
##  the compartments will then be replicated five times
### U: unvaccinated
### B: Bexsero
### X: Vaccine X
### D: Dual
### W: Waned

n_group <- 2 


## Core equations for transitions between compartments:

deriv(U[]) <- enr[i]             - (foi[i]+exr)*U[i] + nu*A[i] + rho*T[i]
deriv(I[]) <- foi[i]*U[i]        - (sigma+exr)*I[i] 
deriv(A[]) <- (1-psi)*sigma*I[i] - (eta+nu+exr)*A[i]
deriv(S[]) <- psi*sigma*I[i]     - (mu+exr)*S[i]
deriv(T[]) <- mu*S[i] + eta*A[i] - (rho+exr)*T[i]

## Update population size
C[] <- I[i] + A[i] + S[i] 
N[] <- U[i] + C[i] + T[i]

# calculate mixing matrix, probability of infection and force of infection
Np[] <- N[i] * p[i]
m[, ] <- epsilon * (i==j) + (1-epsilon) * p[j] * N[j] / sum(Np)
lambda[, ] <- m[i, j] * C[j] / N[j]
foi[] <- beta * p[i] * sum(lambda[i, ])


## Set up compartments
## Initial states are all 0 as we will provide a state vector
initial(U[]) <- 0
initial(I[]) <- 0
initial(A[]) <- 0
initial(S[]) <- 0
initial(T[]) <- 0

# set up dimensions of compartments
dim(U) <- n_group
dim(I) <- n_group
dim(A) <- n_group
dim(S) <- n_group
dim(T) <- n_group

dim(C) <- n_group
dim(N) <- n_group
dim(Np) <- n_group

dim(m) <- c(n_group, n_group)
dim(lambda) <- c(n_group, n_group)
dim(foi) <- c(n_group)

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
output(m) <- m
output(foi) <- foi
