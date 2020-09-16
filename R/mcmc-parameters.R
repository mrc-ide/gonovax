
## create function to reflect proposal boundaries at pars_min and pars_max
## this ensures the proposal is symmetrical and we can simplify the MH step
reflect_proposal <- mcstate:::reflect_proposal

reflect_proposal_both <- mcstate:::reflect_proposal_both

reflect_proposal_one <- mcstate:::reflect_proposal_one
