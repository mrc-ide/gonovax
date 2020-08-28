# run no-vaccination model for 20 years from equilibrium
res <- gonovax::run_novax(tt = seq(0, 20), equilib = TRUE)

# extract cumulative incidence (aggregated over activity group and vax status)
# this returns a time x parameter set matrix
cum_incid <- sapply(res, function(x) apply(x[["cum_incid"]], 1, sum))

# extract incidence from cumulative incidence and remove cumulative incidence at
# time zero
novax_baseline <- list(incid = apply(cum_incid, 2, diff),
                       cum_incid = cum_incid[-1, ])

saveRDS(novax_baseline, "inst/extdata/novax_baseline.rds")
