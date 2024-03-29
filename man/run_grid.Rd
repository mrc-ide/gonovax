% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_grid.R
\name{run_grid}
\alias{run_grid}
\title{run model from equilibrium with single vaccine at the input
efficacy / duration grid locations for n parameter sets}
\usage{
run_grid(
  gono_params,
  init_params,
  cost_params,
  baseline,
  model,
  eff,
  dur,
  vbe = 0,
  strategy = NULL,
  uptake_total = 0,
  uptake_second_dose = uptake_total,
  t_stop = 99,
  full_output = FALSE,
  disc_rate = 0
)
}
\arguments{
\item{gono_params}{gono params}

\item{init_params}{initial state parameters}

\item{cost_params}{cost effectiveness parameters}

\item{baseline}{optional input of a baseline to compare to, must be a
gonovax_grid object if supplied}

\item{model}{a gonovax run function, by default `run_onevax`}

\item{eff}{numeric vector (between 0-1) of vaccine efficacy against
acquisition}

\item{dur}{numeric vector of duration in years of the vaccine}

\item{vbe}{single numeric indicating % of population vaccinated before entry
(between 0-1), default = 0}

\item{strategy, }{default is null, no vaccination
`VoD`: vaccination on diagnosis, `VoA`: vaccination on attendance,
`VoD(L)+VoA(H)`: targeted vaccination (i.e. all diagnosed plus group H on
 screening)}

\item{uptake_total}{numeric (0-1) of strategy uptake, default = 0}

\item{uptake_second_dose}{numeric (0-1) of strategy uptake,
default = uptake_total}

\item{t_stop}{time at which vaccination should stop (years), default = 99}

\item{full_output}{logical indicating whether full results should be output}

\item{disc_rate}{discount rate for cost-effectiveness output
per annum, default = 0}
}
\value{
A `gonovax_grid` object
}
\description{
run model from equilibrium with single vaccine at the input
efficacy / duration grid locations for n parameter sets
}
