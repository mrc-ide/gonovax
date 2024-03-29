% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameters_trial.R
\name{initial_params_trial}
\alias{initial_params_trial}
\title{Initial conditions for the model trial where the entire cohort is
in the high sexual activity group.}
\usage{
initial_params_trial(pars, n_vax = 1, p_v = 1, n_diag_rec = n_diag_rec)
}
\arguments{
\item{pars}{A parameter list containing `N0`, and `q` elements.}

\item{n_vax}{an integer indicating the number of vaccine compartments}

\item{p_v}{a vector of length `n_vax` that sums to 1 denoting the
proportion in each vaccine stratum}

\item{n_diag_rec}{integer giving the number of each X, V(erlang), and W
stratum, allowing tracking of diagnosis history. e.g for a n_diag_rec = 2
and erlang = 1, there will be X.I, X.II, V1.I, V1.II, W.I, W.II strata.
Where '.I' corresponds to never-diagnosed individuals and '.II' is for
individuals diagnosed at least once.}
}
\value{
A list of initial model states
}
\description{
Create initial conditions for the model trial
}
