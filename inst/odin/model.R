# Key to indices
# i = Activity group
## 1: Low
## 2: High
# j = Vaccination status
## 1: Unvaccinated
## 2: Bexsero
## 3: Vaccine X
## 4: Dual
## 5: Waned

n_group <- 2
n_vax   <- user(1)

## Core equations for transitions between compartments:

deriv(U[, ]) <- enr * q[i] * ve[j] - n_UI[i, j] - exr * U[i, j] +
  nu * A[i, j] + n_TU[i, j] -
  sum(n_vd[i, j, ]) - sum(n_vs[i, j, ]) + sum(wU[i, j, ])

deriv(I[, ]) <- n_UI[i, j] - (sigma + exr) * I[i, j] + sum(wI[i, j, ])

deriv(A[, ]) <- (1 - psi) * sigma * I[i, j] - n_AT[i, j] -
  (nu + exr) * A[i, j] + sum(wA[i, j, ])

deriv(S[, ]) <- psi * sigma * I[i, j] - n_ST[i, j] - exr * S[i, j] +
  sum(wS[i, j, ])

deriv(T[, ]) <- n_ST[i, j] + n_AT[i, j] - exr * T[i, j] - n_TU[i, j] +
  sum(wT[i, j, ])

## Update population size
C[, ] <- I[i, j] + A[i, j] + S[i, j]
N[, ] <- U[i, j] + C[i, j] + T[i, j]

# calculate mixing matrix, probability of infection and force of infection
Np[]    <- sum(N[i, ]) * p[i]
m[, ]   <- epsilon * (i == j) + (1 - epsilon) * Np[j] / sum(Np)
foi[, ] <- beta * p[i] * m[i, j] * sum(C[j, ]) / sum(N[j, ])

n_UI[, ]     <- sum(foi[i, ]) * (1 - eff[j]) * U[i, j]
n_AT[, ]     <- eta * A[i, j]
n_ST[, ]   <- mu * S[i, j]
n_TU[, ]     <- rho * T[i, j]
screened[, ] <- eta * U[i, j]

# vaccination
## time-varying switch
vax_switch <- interpolate(vax_t, vax_y, "constant")
## at screening
n_vs[, , ] <- vs[i, j, k] * screened[i, k] * vax_switch
## on diagnosis
n_vd[, , ] <- vd[i, j, k] * n_TU[i, k] * vax_switch
## on entry
n_ve[, 2:n_vax] <- enr * q[i] * ve[j]

# waning
wU[, , ] <- w[j, k] * U[i, k]
wI[, , ] <- w[j, k] * I[i, k]
wA[, , ] <- w[j, k] * A[i, k]
wS[, , ] <- w[j, k] * S[i, k]
wT[, , ] <- w[j, k] * T[i, k]

## outputs

deriv(cum_incid[, ])      <- n_UI[i, j]
deriv(cum_diag_a[, ])     <- n_AT[i, j]
deriv(cum_diag_s[, ])     <- n_ST[i, j]
deriv(cum_treated[, ])    <- n_TU[i, j]
deriv(cum_screened[, ])   <- screened[i, j]
deriv(cum_vaccinated[, ]) <- n_vs[i, j, j] + n_vd[i, j, j] + n_ve[i, j]

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
initial(cum_vaccinated[, ]) <- 0

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
dim(Np) <- n_group

dim(m)   <- c(n_group, n_group)
dim(foi) <- c(n_group, n_group)

dim(n_UI)     <- c(n_group, n_vax)
dim(n_AT)     <- c(n_group, n_vax)
dim(n_ST)     <- c(n_group, n_vax)
dim(n_TU)     <- c(n_group, n_vax)
dim(screened) <- c(n_group, n_vax)

dim(cum_incid)      <- c(n_group, n_vax)
dim(cum_diag_a)     <- c(n_group, n_vax)
dim(cum_diag_s)     <- c(n_group, n_vax)
dim(cum_treated)    <- c(n_group, n_vax)
dim(cum_screened)   <- c(n_group, n_vax)
dim(cum_vaccinated) <- c(n_group, n_vax)

## Parameters
p[]     <- user()
q[]     <- user()

enr     <- user()
exr     <- user()
beta    <- user()
epsilon <- user()
sigma   <- user()
psi     <- user()
eta     <- user()
nu      <- user()
mu      <- user()
rho     <- user()

## vaccination pars
ve[]   <- user()
vs[, , ] <- user()
vd[, , ] <- user()
eff[]  <- user()
w[, ]  <- user()
vax_t[] <- user()
vax_y[] <- user()

## par dimensions
dim(p)    <- n_group
dim(q)    <- n_group
dim(ve)   <- n_vax
dim(eff)  <- n_vax
dim(vd)   <- c(n_group, n_vax, n_vax)
dim(vs)   <- c(n_group, n_vax, n_vax)
dim(w)    <- c(n_vax, n_vax)
dim(vax_t) <- user()
dim(vax_y) <- user()

dim(n_ve) <- c(n_group, n_vax)
dim(n_vs) <- c(n_group, n_vax, n_vax)
dim(n_vd) <- c(n_group, n_vax, n_vax)
dim(wU)   <- c(n_group, n_vax, n_vax)
dim(wI)   <- c(n_group, n_vax, n_vax)
dim(wA)   <- c(n_group, n_vax, n_vax)
dim(wS)   <- c(n_group, n_vax, n_vax)
dim(wT)   <- c(n_group, n_vax, n_vax)

output(N)   <- N
output(foi) <- foi
