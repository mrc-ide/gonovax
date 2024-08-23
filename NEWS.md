# gonovax 0.4.20
* person years exposed from the point of view of a non-screening RCT (U, I, A, S, AND T) added

# gonovax 0.4.20
* person-years exposed both from the point of view of an RCT (U, I, A, S) and in the true sense (U only) are output from the model rather than calculated in post

# gonovax 0.4.19
* partner notification and vaccination on notification added into population-level model

# gonovax 0.4.18
* diagnoses history in a trial is recorded at the point of receiving treatment and whether a trial records asymptomatic infections or not can be switched on or off

# gonovax 0.4.17
* ved acts on duration of asymptomatic infection rather than the rate of asymptomatic natural clearance

# gonovax 0.4.16
* Expose demographic_params in run function to facilitate adaptation to non-UK contexts

# gonovax 0.4.15
* the maximum effective rate of natural clearance in the presence of ved is now capped at the value for the rate of care seeking

# gonovax 0.4.14
* extracting diagnoses and vaccination by stratum type now possible with n_erlang > 1, and symptomatic and
asymptomatic diagnoses by vaccine-protected vs non-vaccine-protected can be extracted against (adjusted) baseline

# gonovax 0.4.13
* model updated to include diagnosis history (i.e. the number of times individuals have been diagnosed in last Y years)
at n_diag_rec = 1, the model is the same as previous versions
tested to work for n_diag_rec > 1 with xpvwrh 
tested to work with n_diag_rec = 1 for all other strategies

# gonovax 0.4.12
* model_trial_stochastic updated so that total number of individuals exiting each state drawn first, then the specific transitions made drawn after

#gonovax 0.4.11
* deterministic version of the xvw_trial model can now have parallel strata that account for diagnosis history

# gonovax 0.4.10
* stochastic version of the xvw_trial model can now have parallel strata that account for diagnosis history

# gonovax 0.4.9 
* package now also includes a stochastic version of the xvw_trial model 


# gonovax 0.4.8
* model XPVWRH and XVW_trial can now have 'n' transitions between entering a vaccine-protected stratum and waning to no protection (erlang compartments)

# gonovax 0.4.7
* model outputs present value of cumulative cases averted in a model run compared with baseline transmission in the absence of 
vaccination

# gonovax 0.4.6
* model outputs level of population vaccine protection (number vaccinated weighted by efficacy). Tests also added to check number of individuals
vaccinated is being performed correctly

# gonovax 0.4.5
* model outputs annual number of individuals and proportion of the population 
experiencing vaccine protection (not just the number vaccinated)

# gonovax 0.4.4
* added incremental costs calculation for £35 and £70 per dose

# gonovax 0.4.3
* new downstream functions generate meaningful output metrics for branching
model

# gonovax 0.4.2

* new model structure where primary vaccination can either be partial or full depending on 1 or 2 doses, with primary vaccination waning to unvaccinated and an immunologically naive state, and full vaccination waning with some immunological memory. Partial vaccination with one dose can now provide an assigned level of protection. Individuals with one dose who visit the clinic
can receive their second dose.

# gonovax 0.4.1

* calculate and return incremental vaccination costs for £9 per dose

# gonovax 0.4.0

* calculate and return annual and cumulative # individuals receiving primary and booster vaccination separately and annual and cumulative number of primary and booster doses administered 

# gonovax 0.3.12

* calculate and return # doses wasted. Outputs incremental costs for intermediate price (£50)

# gonovax 0.3.11

* output the cumulative number offered vaccination

# gonovax 0.3.10

* document and export compare_baseline function

# gonovax 0.3.9

* allows proportion of population to be vaccine hesitant, allows primary and booster vaccination uptakes to be different proportions

# gonovax 0.3.8

* add vaccination of trial cohort

# gonovax 0.3.7

* Output costs and benefits separately to allow ICER calculation

# gonovax 0.3.6

* Assume two doses for adolescent vaccination

# gonovax 0.3.5

* output incremental cumulative vaccine doses given

# gonovax 0.3.4

* allow setting of initial coverage in XVW model

# gonovax 0.3.3

* improve tests and add checks for bad parameter inputs

# gonovax 0.3.2

* add three-stratum waning model with no re-vaccination: X -> V -> W 

# gonovax 0.3.1

* add vignette
* add fixed parameter transform function

# gonovax 0.3.0

* add vaccine effects vea (acquisition), vei (infectiousness),
ved (duration of infection), ves (symptoms)
* rename ve -> vbe to avoid confusion

# gonovax 0.2.9

* reduce input size from run_grid

# gonovax 0.2.8

* rename strategies

# gonovax 0.2.7

* fix bug QALY cost bug

# gonovax 0.2.6

* move transform function into package

# gonovax 0.2.5

* run_grid uses length of baseline to determine time horizon

# gonovax 0.2.4

* Allow efficacy to be varied in run_onevax

# gonovax 0.2.3

* move compare function used in paper into package

# gonovax 0.2.2

* refactor to run based on strategy
* allow XVWR model to have two different vaccines (i.e. eff / dur profile)

# gonovax 0.2.1

* add four-stratum waning model: X -> V -> W <-> R

# gonovax 0.2.0

* move globals into package

# gonovax 0.1.24

* allow model restart from a future point (necessary for time-varying params)

# gonovax 0.1.23

* add GRASP 2020 data

# gonovax 0.1.22

* add varying screening rate by group

# gonovax 0.1.21

* data starts in 2010
* add beta binomial function

# gonovax 0.1.20

* add time-varying beta and eta
* add option to change compare function

# gonovax 0.1.19

* add GRASP % symptomatic data
* export data easily
* add user function for running novax

# gonovax 0.1.18

* add mcmc thinning and sampling tools

# gonovax 0.1.17

* add mcmc fitting method to gumcad data

# gonovax 0.1.16

* bug-fix: extract flows now correctly attributes vaccination / re-vaccination
* re-code model so treatment of VbE is now the same structure as VoD / VoA

# gonovax 0.1.15

* output vaccination and revaccination separately
* use parsimonious numbers of compartments in waning model

# gonovax 0.1.14

* output full time series in run_grid()
* extend baseline to 30 years

# gonovax 0.1.13
* add time-varying vaccination parameter

# gonovax 0.1.12

* export aggregate() 
* export set_strategy()

# gonovax 0.1.11

* output asymptomatic and symptomatic diagnoses

# gonovax 0.1.10

* extend baseline to 20 years
* take plot method out of package
* export format_grid

# gonovax 0.1.9

* allow vaccination of activity groups at differing rates

# gonovax 0.1.8

* add onevax_waning model

# gonovax 0.1.7

* remove n_par index as running multiple parameter sets in vectorised mode is
inefficient in ode solver.

# gonovax 0.1.6

* add option to output full results in grid_search()
* add ability to input baseline in grid_search()

# gonovax 0.1.5

* fix bug in novax_baseline()
* fix model bug that caused people to go missing in VoD
* update plot method to avoid ggplot2 warnings

# gonovax 0.1.4

* cache baseline (novax) incidence parameters
* add run_grid_onevax() which run the model with cached parameters across an
eff x dur vaccination grid
* add plot methods for grid objects

# gonovax 0.1.3

* rename model to separate novax, onevax (and in future, dualvax / dualvax_wane)
* add run_novax() and run_onevax() which run the model with cached parameters
from equilibrium if specified
* add VbE to vaccination output

# gonovax 0.1.2

* add novax_equilib() which returns the cached equilibrium compartments for all
parameters in gono_params()

# gonovax 0.1.1

* export gono_params() and allow selection of specific parameter sets
* export dualvax_initial()
* add dualvax_restart() for restarting the model from a previous position

# gonovax 0.1.0

* initial version of the model tested with Bexsero only
