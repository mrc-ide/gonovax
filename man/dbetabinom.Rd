% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/support.R
\name{dbetabinom}
\alias{dbetabinom}
\title{pdf of a betabinomial parametrised in terms of probability and
over-dispersion}
\usage{
dbetabinom(x, size, prob, rho, log = FALSE)
}
\arguments{
\item{x}{data}

\item{size}{integer of sample size}

\item{prob}{probability of observing a single success}

\item{rho}{overdispersion parameter, has support [0, 1]
with low values being less overdispersion}

\item{log}{logical indicating whether to return log value}
}
\value{
probability of observing x
}
\description{
pdf of a betabinomial parametrised in terms of probability and
over-dispersion
}
