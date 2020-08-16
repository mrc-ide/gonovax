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
