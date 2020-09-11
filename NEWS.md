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
