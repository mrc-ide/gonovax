% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/onevax_xpvwrh.R
\name{run_onevax_xpvwrh}
\alias{run_onevax_xpvwrh}
\title{run model with a two-dose vaccine for input parameter sets, either
from initialisation or from equilibrium, those with waned vaccines are
eligible for revaccination (R), and return to the R stratum, those with
waned partial vaccines return to the unvaccinated stratum (X) and considered
immunologically naive. A user defined proportion of the population is
vaccine hesitant and is never vaccinated. Full acciantion (V) with 2 doses
gives maximum protection whereas partial vaccination with 1 dose (P) gives
less. Individuals can get 2 doses either by committing in X or upon visiting
a clinic in P.}
\usage{
run_onevax_xpvwrh(
  tt,
  gono_params,
  init_params = NULL,
  dur_v = 1000,
  dur_p = dur_v,
  vea = 0,
  vei = 0,
  ved = 0,
  ves = 0,
  dur_revax = dur_v,
  vea_revax = vea,
  vei_revax = vei,
  ved_revax = ved,
  ves_revax = ves,
  vea_p = vea,
  vei_p = vei,
  ved_p = ved,
  ves_p = ves,
  vbe = 0,
  r1 = 0,
  r2 = 0,
  r2_p = 0,
  booster_uptake = (r1 * r2),
  strategy = NULL,
  t_stop = 99,
  hes = 0,
  n_erlang = 1,
  n_diag_rec = 1,
  years_history = 1
)
}
\arguments{
\item{tt}{a numeric vector of times at which the model state is output}

\item{gono_params}{list of gono params}

\item{init_params}{A list of starting conditions}

\item{dur_v}{duration of time spent in V stratum after completing a round of
primary vaccination (fully vaccinated, accepting first and second dose)}

\item{dur_p}{duration of time spent in the P stratum, partially vaccinated
(accepting only the first dose)}

\item{vea}{scalar or numeric vector with same length as `gono_params` giving
efficacy of the vaccine against acquisition (between 0-1)}

\item{vei}{scalar or numeric vector with same length as `gono_params` giving
efficacy of the vaccine against infectiousness (between 0-1)}

\item{ved}{scalar or numeric vector with same length as `gono_params` giving
efficacy of the vaccine against duration (between 0-1)}

\item{ves}{scalar or numeric vector with same length as `gono_params` giving
efficacy of the vaccine against symptoms (between 0-1)}

\item{dur_revax}{scalar or numeric vector with same length as `gono_params`
giving duration of protection for revaccination, default to same as primary}

\item{vea_revax}{scalar or numeric vector with same length as `gono_params`
giving efficacy of revaccination against acquisition, default to same as
primary}

\item{vei_revax}{scalar or numeric vector with same length as `gono_params`
giving efficacy of revaccination against infectiousness, default to same as
primary}

\item{ved_revax}{scalar or numeric vector with same length as `gono_params`
giving efficacy of revaccination against duration of infection, default to
same as primary}

\item{ves_revax}{scalar or numeric vector with same length as `gono_params`
giving efficacy of revaccination against symptoms, default to same as
primary}

\item{vea_p}{scalar indicating efficacy of partial vaccination against
acquisition (between 0-1)}

\item{vei_p}{scalar indicating efficacy of partial vaccination against
infectiousness (between 0-1)}

\item{ved_p}{scalar indicating efficacy of partial vaccination against
duration (between 0-1)}

\item{ves_p}{scalar indicating efficacy of partial vaccination against
symptoms (between 0-1)}

\item{vbe}{scalar indicating pc of population vaccinated before entry
(between 0-1)}

\item{r1}{scalar or numeric vector with same length as
'gono_params' giving proportion of population offered vaccine only
 accepting the first dose, becoming partially vaccinated}

\item{r2}{scalar or numeric vector with same length as
'gono_params' giving proportion of the population who accepted the first
 dose of the vaccine who go on to accept the second dose, becoming fully
 vaccinated}

\item{r2_p}{scalar or numeric vector with same length as 'gono_params'
giving proportion of partially vaccinated individuals who later receive
a second dose when returning to the clinic due to screening or illness}

\item{booster_uptake}{scalar or numeric vector with same length as
'gono_params' giving proportion of population undertaking booster
vaccination after primary vaccination protection has waned.
 Defaults to supplied value of r1 * r2}

\item{strategy}{single character string in "VoD", "VoD(H)",
"VoA", "VoA(H)", "VoD(L)+VoA(H)". Defaults to NULL i.e. no vaccination}

\item{t_stop}{time at which vaccination should stop (years)}

\item{hes}{Proportion of individuals in the population who are vaccine
hesitant}

\item{n_erlang}{integer giving the number of erlang vaccination transitions
through vaccine-protected strata until that protection has waned}

\item{n_diag_rec}{integer for the number of diagnosis history substrata}

\item{years_history}{number of years that diagnosis history is recorded for}
}
\value{
A list of transformed model outputs
}
\description{
run model with a two-dose vaccine for input parameter sets, either
from initialisation or from equilibrium, those with waned vaccines are
eligible for revaccination (R), and return to the R stratum, those with
waned partial vaccines return to the unvaccinated stratum (X) and considered
immunologically naive. A user defined proportion of the population is
vaccine hesitant and is never vaccinated. Full acciantion (V) with 2 doses
gives maximum protection whereas partial vaccination with 1 dose (P) gives
less. Individuals can get 2 doses either by committing in X or upon visiting
a clinic in P.
}
