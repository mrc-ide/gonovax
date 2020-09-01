# run no-vaccination model for 20 years from equilibrium
res <- gonovax::run_novax(tt = seq(0, 20), equilib = TRUE)

# extract cumulative incidence (aggregated over activity group and vax status)
# this returns a time x parameter set matrix

cum_incid  <- aggregate(res, "cum_incid")
cum_diag_a <- aggregate(res, "cum_diag_a")
cum_diag_s <- aggregate(res, "cum_diag_s")

# extract incidence from cumulative incidence and remove cumulative incidence at
# time zero
novax_baseline <- list(incid = apply(cum_incid, 2, diff),
                       cum_incid = cum_incid[-1, ],
                       cum_diag_a = cum_diag_a[-1, ],
                       cum_diag_s = cum_diag_s[-1, ])

saveRDS(novax_baseline, "inst/extdata/novax_baseline.rds")
